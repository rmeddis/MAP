function [P, BFlist, sacf, boundaryValue] = ...
    filteredSACF(inputSignalMatrix, method, params)
% UTIL_filteredSACF computes within-channel, running autocorrelations (acfs)
%  and finds the sum across channels (sacf).
%  The SACF is smoothed to give the 'p function' (P).
%
% INPUT
%  	inputSignalMatrix: a matrix (channel x time) of AN activity
%  	method.dt:			the signal sampling interval in seconds
%   method.segmentNo:
%   method.nonlinCF
%   
% params: 	a list of parmeters to guide the computation:
%   filteredSACFParams.lags: an array of lags to be computed (seconds)
%   filteredSACFParams.acfTau: time constant (sec) of the running acf calculations
%     if acfTau>1 it is assumed that Wiegrebe'sacf method
%   	for calculating tau is to be used (see below)
%   filteredSACFParams.Lambda: time constant  for smoothing thsacf to make P
%   filteredSACFParams.lagsProcedure identifies a strategy for omitting some lags.
%    Options are: 'useAllLags', 'omitShortLags', or 'useBernsteinLagWeights'
%   filteredSACFParams.usePressnitzer applies lower weights longer lags
%   parafilteredSACFParamsms.plotACFs (=1) creates movie of acf matrix (optional)
%
%
% OUTPUT
%  	P:	P function 	(time x lags), a low pass filtered version of sacf.
%  method: updated version of input method (to include lags used)
%   sacf:  	(time x lags)

%%
boundaryValue=[];
[nChannels inputLength]= size(inputSignalMatrix);
% list of BFs must be repeated is many fiber types used
BFlist=method.nonlinCF;
nfibertypes=nChannels/length(BFlist);
BFlist=repmat(BFlist',2,1)';

dt=method.dt;
% adjust sample rate, if required
if isfield(params,'dt')
    [inputSignalMatrix dt]=UTIL_adjustDT(params.dt, method.dt, inputSignalMatrix);
    method.dt=dt;
end

% create acf movie
if isfield(params, 'plotACFs') && params.plotACFs==1
    plotACF=1;
else
    plotACF=0;  % default
end

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
% disp(['lag procedure= ''' lagsProcedure ''''])
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
    case 'useBernsteinLagWeights'
        lagWeights=BernsteinLags(BFlist, lags)';
    otherwise
        error ('acf: params.lagProcedure not recognised')
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
else                      % use Wiegrebe method: tau= 2*lag (for example)
    WiegrebeFactor=acfTau;
    acfTaus=WiegrebeFactor*lags;
    idx= acfTaus<0.0025; acfTaus(idx)=0.0025;
    acfDecayFactors= exp(-dt./(acfTaus));
end
% make acfDecayFactors into a (channels x lags) matrix for speedy computation
acfDecayFactors=repmat(acfDecayFactors,nChannels, 1);

% P function lowpass filter decay (only one value needed)
pDecayFactor=exp(-dt/params.lambda);

% ACF
% lagPointers is a list of pointers relative to 'time now'
lagPointers=round(lags/dt);
if max(lagPointers)+1>inputLength
    error([' filteredSACF: not enough signal to evaluate ACF. Max(lag)= ' num2str(max(lags))])
end


P=zeros(nLags,inputLength+1);   % P must match segment length +1
sacf=zeros(nLags,inputLength);
if   ~isfield(method,'segmentNumber') || method.segmentNumber==1
    acf=zeros(nChannels,nLags);
    % create a runup buffer of signal
    buffer= zeros(nChannels, max(lagPointers));
else
    % boundaryValue picks up from a previous calculation
    acf=params.boundaryValue{1};
    P(: , 1)=params.boundaryValue{2}; % NB first value is last value of previous segment
    buffer=params.boundaryValue{3};
end
inputSignalMatrix=[buffer inputSignalMatrix];
[nChannels inputLength]= size(inputSignalMatrix);

timeCounter=0; biggestSACF=0;
for timePointer= max(lagPointers)+1:inputLength
    % acf is a continuously changing channels x lags matrix
    %   Only the current value is stored
    % sacf is the vertical summary of acf ( a vector) and all values are kept and returned
    % P is the smoothed version of sacf and all values are kept and returned
    % lagWeights emphasise some BF/lag combinations and ignore others
    % NB time now begins at the longest lag.
    % E.g. if max(lags) is .04 then this is when the ACf will begin.
    %            AN                                       AN delayed                             weights               filtering
    
    % This is the ACF calculation
    timeCounter=timeCounter+1;
    acf= (repmat(inputSignalMatrix(:, timePointer), 1, nLags) .* inputSignalMatrix(:, timePointer-lagPointers)).*lagWeights *dt + acf.* acfDecayFactors;
    x=(mean(acf,1)./acfTaus)';
%     disp(num2str(x'))
    sacf(:,timeCounter)=x;
    P(:,timeCounter+1)=sacf(:,timeCounter)*(1-pDecayFactor)+P(:,timeCounter)*pDecayFactor;
    
    % plot at intervals of 200 points
    if plotACF && ~mod(timePointer,params.plotACFsInterval)
        %       mark cursor on chart to signal progress
        % this assumes that the user has already plotted
        % the signal in subplot(2,1,1) of figure (13)
        figure(13)
        hold on
        subplot(4,1,1)
        time=timePointer*dt;
        a =ylim;
        plot([time time], [a(1) a(1)+(a(2)-a(1))/4]) % current signal point marker
        
        %         plot ACFs one per channel
        subplot(2,1,2), cla
        cascadePlot(acf, lags, BFlist)
        xlim([min(lags) max(lags)])
        %         set(gca,'xscale','log')
        title(num2str(method.dt*timePointer))
        ylabel('BF'), xlabel('period (lag)')
        
        %         plot SACF
        subplot(4,1,2), hold off
        plot(lags,sacf(:,timeCounter)-min(sacf(:,timeCounter)))
        biggestSACF=max(biggestSACF, max(sacf(:,timeCounter)));
        if biggestSACF>0, ylim([0 biggestSACF]), end
        %         set(gca,'xscale','log')
        title('SACF')
        pause(params.plotMoviePauses)
    end
end
P=P(:,1:end-1);  % correction for speed up above

% Pressnitzer weights
if ~isfield(params, 'usePressnitzer'),     params.usePressnitzer=0; end
if params.usePressnitzer
    [a nTimePoints]=size(P);
    % higher pitches get higher weights
    %     PressnitzerWeights=repmat(min(lags)./lags,nTimePoints,1);
    % weaker weighting
    PressnitzerWeights=repmat((min(lags)./lags).^0.5, nTimePoints,1);
    P=P.*PressnitzerWeights';
    sacf=sacf.*PressnitzerWeights';
end

% wrap up
method.acfLags=lags;
method.filteredSACFdt=dt;

boundaryValue{1}=acf(:,end-nLags+1:end);
boundaryValue{2}=P(:,end);
% save signal buffer for next segment
boundaryValue{3} = inputSignalMatrix(:, end-max(lagPointers)+1 : end);

method.displaydt=method.filteredSACFdt;

% if ~isfield(params, 'plotUnflteredSACF'), params.plotUnflteredSACF=0; end
% if method.plotGraphs
%     method.plotUnflteredSACF=params.plotUnflteredSACF;
%     if ~method.plotUnflteredSACF
%         method=filteredSACFPlot(P,method);
%     else
%         method=filteredSACFPlot(SACF,method);
%     end
% end


% ------------------------------------------------  plotting ACFs
function cascadePlot(toPlot, lags, BFs)
% % useful code
[nChannels nLags]=size(toPlot);

% cunning code to represent channels as parallel lines
[nRows nCols]=size(toPlot);
if nChannels>1
    % max(toPlot) defines the spacing between lines
    a=max(max(toPlot))*(0:nRows-1)';
    % a is the height to be added to each channel
    peakGain=10;
    % peakGain emphasises the peak height
    x=peakGain*toPlot+repmat(a,1,nCols);
    x=nRows*x/max(max(x));
else
    x=toPlot;                            % used when only the stimulus is returned
end
plot(lags, x','k')
ylim([0 nRows])

