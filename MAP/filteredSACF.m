function [LP_SACF, BFlist, SACF, ACFboundaryValue] = ...
    filteredSACF(inputSignalMatrix, dt, BFlist, params)
% UTIL_filteredSACF computes within-channel, running autocorrelations (ACFs)
%  and finds the running sum across channels (SACF).
%  The SACF is smoothed to give the 'p function' (LP_SACF).
%  (Balaguer-Ballestera, E. Denham, S.L. and Meddis, R. (2008))
%
%
% INPUT
%  	inputSignalMatrix: a matrix (channel x time) of AN activity
%  	dt:		the signal sampling interval in seconds
%   BFlist: channel BFs
%
% params: 	a list of parmeters to guide the computation:
%   params.lags: an array of lags to be computed (seconds)
%   params.acfTau: time constant (sec) of the running ACF
%     if acfTau>1 it is assumed that Wiegrebe'SACF method
%   	for calculating tau is to be used (see below)
%   params.Lambda: smoothing factor for the SACF
%   params.lagsProcedure: strategies for omitting some lags.
%    (Options: 'useAllLags' or 'omitShortLags')
%   params.usePressnitzer applies lower weights longer lags
%   params.plotACFs (=1) creates movie of ACF matrix (optional)
%
% OUTPUT
%  	LP_SACF:	LP_SACF function 	(time x lags), a low pass filtered version of SACF.
%  method: updated version of input 'method' (to include lags used)
%   SACF:  	(time x lags)
%
% Notes:
% ACFboundaryValue refers to segmented evaluation and is currently not
%  supported. However the code may be useful later when this function
%  is incorporated into MAP1_14.

%%
global savedInputSignal

ACFboundaryValue=[];

[nChannels inputLength]= size(inputSignalMatrix);

% create ACF movie
if isfield(params, 'plotACFs') && params.plotACFs==1
    plotACF=1;
    signalTime=dt:dt:dt*length(savedInputSignal);
else
    plotACF=0;  % default
end
params.plotACFsInterval=round(params.plotACFsInterval/dt);

if isfield(params,'lags')
    lags=params.lags;
else
    lags=params.minLag:params.lagStep:params.maxLag;
end
nLags=length(lags);

% Establish lag weightsaccording to params.lagsProcedure
%  lagWeights allow some lag computations to be ignored
%  lagWeights is a (channel x lag) matrix;
if isfield(params, 'lagsProcedure')
    lagsProcedure=params.lagsProcedure;
else
    lagsProcedure='useAllLags';  % default
end

lagWeights=ones(nChannels,nLags);
switch lagsProcedure
    case 'useAllLags'
        % no action required lagWeights set above
    case 'omitShortLags'
        % remove lags that are short relative to CF
        allLags=repmat(lags,nChannels,1);
        allCFs=repmat(BFlist',1,nLags);
        criterionForOmittingLags=1./(params.criterionForOmittingLags*allCFs);
        idx= allLags < criterionForOmittingLags;	% ignore these lags
        lagWeights(idx)=0;
    otherwise
        error ('ACF: params.lagProcedure not recognised')
end


% Establish matrix of lag time constants
%   these are all the same if tau<1
% if acfTau>1, it is assumed that we are using the Wiegrebe method
%   and a different decay factor is applied to each lag
%   ( i.e., longer lags have longer time constants)
acfTau=params.acfTau;
if acfTau<1          % all taus are the same
    acfTaus=repmat(acfTau, 1, nLags);
    acfDecayFactors=ones(size(lags)).*exp(-dt/acfTau);
else                  % use Wiegrebe method: tau= 2*lag (for example)
    WiegrebeFactor=acfTau;
    acfTaus=WiegrebeFactor*lags;
    idx= acfTaus<0.0025; acfTaus(idx)=0.0025;
    acfDecayFactors= exp(-dt./(acfTaus));
end
% make acfDecayFactors into a (channels x lags) matrix for speedy computation
acfDecayFactors=repmat(acfDecayFactors,nChannels, 1);

% LP_SACF function lowpass filter decay (only one value needed)
pDecayFactor=exp(-dt/params.lambda);

% ACF
% lagPointers is a list of pointers relative to 'time now'
lagPointers=round(lags/dt);
if max(lagPointers)+1>inputLength
    error([' filteredSACF: not enough signal to evaluate ACF. Max(lag)= ' num2str(max(lags))])
end


LP_SACF=zeros(nLags,inputLength+1);   % LP_SACF must match segment length +1
SACF=zeros(nLags,inputLength);
method=[];  % legacy programming
if   ~isfield(method,'segmentNumber') || method.segmentNumber==1
    ACF=zeros(nChannels,nLags);
    % create a runup buffer of signal
    buffer= zeros(nChannels, max(lagPointers));
else
    % ACFboundaryValue picks up from a previous calculation
    ACF=params.ACFboundaryValue{1};
    % NB first value is last value of previous segment
    LP_SACF(: , 1)=params.ACFboundaryValue{2};
    buffer=params.ACFboundaryValue{3};
end
inputSignalMatrix=[buffer inputSignalMatrix];
[nChannels inputLength]= size(inputSignalMatrix);

timeCounter=0; biggestSACF=0;
for timePointer= max(lagPointers)+1:inputLength
    % ACF is a continuously changing channels x lags matrix
    %   Only the current value is stored
    % SACF is the vertical sum of ACFs; all values are kept and returned
    % LP_SACF is the smoothed version of SACF:all values are kept and returned
    % lagWeights emphasise some BF/lag combinations and ignore others
    % NB time now begins at the longest lag.
    % E.g. if max(lags) is .04 then this is when the ACf will begin (?).

    % This is the ACF calculation
    timeCounter=timeCounter+1;
    ACF= (repmat(inputSignalMatrix(:, timePointer), 1, nLags) .* ...
        inputSignalMatrix(:, timePointer-lagPointers)).*...
        lagWeights *dt + ACF.* acfDecayFactors;
    x=(mean(ACF,1)./acfTaus)';
    SACF(:,timeCounter)=x;
    
    % smoothed version of SACF
    LP_SACF(:,timeCounter+1)=SACF(:,timeCounter)* (1-pDecayFactor) + ...
        LP_SACF(:,timeCounter)* pDecayFactor;

    % plot ACF at intervals if requested to do so
    if plotACF && ~mod(timePointer,params.plotACFsInterval) && ...
            timePointer*dt>3*max(lags)
        figure(89), clf
        %           plot ACFs one per channel
        subplot(2,1,1)
        UTIL_cascadePlot(ACF, lags)
        title(['running ACF  at ' num2str(dt*timePointer,'%5.3f') ' s'])
        ylabel('channel BF'), xlabel('period (lag, ms)')
        set(gca,'ytickLabel',[])

        %           plot SACF
        subplot(4,1,3), cla
        plotSACF=SACF(:,timeCounter)-min(SACF(:,timeCounter));
        plot(lags*1000, plotSACF, 'k')
        biggestSACF=max(biggestSACF, max(plotSACF));
        if biggestSACF>0, ylim([0 biggestSACF]), else ylim([0 1]), end
        ylabel('SACF'), set(gca,'ytickLabel',[])
%         xlim([min(lags*1000) max(lags*1000)])

        %           plot signal
        subplot(4,1,4)
        plot(signalTime, savedInputSignal, 'k'), hold on
        xlim([0 max(signalTime)])
        a= ylim;
        %       mark cursor on chart to indicate progress
        time=timePointer*dt;
        plot([time time], [a(1) a(1)+(a(2)-a(1))/4], 'r', 'linewidth', 5)
        pause(params.plotMoviePauses)
    end
end
LP_SACF=LP_SACF(:,1:end-1);  % correction for speed up above

% Pressnitzer weights
if ~isfield(params, 'usePressnitzer'),     params.usePressnitzer=0; end
if params.usePressnitzer
    [a nTimePoints]=size(LP_SACF);
    % higher pitches get higher weights
    %     PressnitzerWeights=repmat(min(lags)./lags,nTimePoints,1);
    % weaker weighting
    PressnitzerWeights=repmat((min(lags)./lags).^0.5, nTimePoints,1);
    LP_SACF=LP_SACF.*PressnitzerWeights';
    SACF=SACF.*PressnitzerWeights';
end

% wrap up
% BoundaryValue is legacy programming used in segmented model (keep)
ACFboundaryValue{1}=ACF(:,end-nLags+1:end);
ACFboundaryValue{2}=LP_SACF(:,end);
% save signal buffer for next segment
ACFboundaryValue{3} = inputSignalMatrix(:, end-max(lagPointers)+1 : end);
