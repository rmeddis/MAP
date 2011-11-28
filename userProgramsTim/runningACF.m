function sampledacf = runningACF(inputSignalMatrix, sfreq, params)
% runningACF computes within-channel, running autocorrelations (acfs)
%  and samples it at given time steps
%
% INPUT
%  	inputSignalMatrix: a matrix (channel x time) of AN activity
%  	sfreq: sampling frequency
%   
% params: 	a list of parmeters to guide the computation:
%   params.lags: an array of lags to be computed (seconds)
%   params.acfTau: time constant (sec) of the running acf calculations
%     
%
%
% OUTPUT
%  	sampledacf:  	(time x lags)

%%
boundaryValue=[];
[nChannels inputLength]= size(inputSignalMatrix);
% list of BFs must be repeated is many fiber types used
%BFlist=method.nonlinCF;
%nfibertypes=nChannels/length(BFlist);
%BFlist=repmat(BFlist',2,1)';

dt=1/sfreq;

samplelength=0.003; %sample length in ms
gridpoints=round([samplelength:samplelength:inputLength*dt]/dt);
gridpoints=[1 gridpoints];

lags=params.lags;

nLags=length(lags);


% Establish matrix of lag time constants 
%   these are all the same if tau<1

if params.acfTau<1          % all taus are the same
    acfTaus=repmat(params.acfTau, 1, nLags);
    acfDecayFactors=ones(size(lags)).*exp(-dt/params.acfTau);
else                     
    error('acfTau must be <1');
end
% make acfDecayFactors into a (channels x lags) matrix for speedy computation
acfDecayFactors=repmat(acfDecayFactors,nChannels, 1);

% P function lowpass filter decay (only one value needed)
%pDecayFactor=exp(-dt/params.lambda);

% ACF
% lagPointers is a list of pointers relative to 'time now'
lagPointers=round(lags/dt);
lagWeights=ones(nChannels,nLags);
if max(lagPointers)+1>inputLength
    error([' filteredSACF: not enough signal to evaluate ACF. Max(lag)= ' num2str(max(lags))])
end



    acf=zeros(nChannels,nLags);
    sampledacf = zeros(length(gridpoints),nChannels,nLags);
    % create a runup buffer of signal
    buffer= zeros(nChannels, max(lagPointers));

inputSignalMatrix=[buffer inputSignalMatrix];
[nChannels inputLength]= size(inputSignalMatrix);

timeCounter=0; biggestSACF=0;
gridcounter =1;
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
    acf = (repmat(inputSignalMatrix(:, timePointer), 1, nLags) .* inputSignalMatrix(:, timePointer-lagPointers)).*lagWeights *dt + acf.* acfDecayFactors;
    
    %just sample on certain samplepoints
    if ~isempty(find(timeCounter==gridpoints))
        sampledacf(gridcounter,:,:)=acf;
        gridcounter = gridcounter+1;
    end
    
    end
end




