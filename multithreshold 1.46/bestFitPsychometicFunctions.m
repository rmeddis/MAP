% ------------------------------------------------------- bestFitPsychometicFunctions
function [psy, levelsPhaseTwoBinVector, logistic, rareEvent]= bestFitPsychometicFunctions (levelsPhaseTwo, responsesPhaseTwo)
% bestFitPsychometicFunctions computes a psychometric function from a matrix of levels and responses
%output
%  psy is the empirical probabbility of a detection associated with levels in levelsPhaseTwoBinVector
%  logistic is a structure with results for fitting the logistic function
%    the logistic fit depends on whether maximum likelihood is used or least squares
%  rareEvent is a structure with results for fitting the rareEvent function
%   this is always calculated and lest squares is always used

global experiment stimulusParameters binFrequencies

% Generate a psychometic function by binning the levelsPhaseTwo
% x is a vector of [levelsPhaseTwo; responsesPhaseTwo]
% experiment.psyBinWidth defines the width of the bin
[psy, levelsPhaseTwoBinVector, binFrequencies,yesResponses, noResponses]= ...
    psychometricFunction(levelsPhaseTwo, responsesPhaseTwo, experiment.psyBinWidth);
% [levelsPhaseTwoBinVector; binFrequencies; psy];

% undefined slope - return
if isempty(psy)
    % slope is undefined
    psy=[]; levelsPhaseTwoBinVector=[];
    logistic.bestThreshold=NaN; logistic.bestK=NaN;  logistic.predictionsLOG=[]; logistic.predictionLevels=[];
    rareEvent.bestGain=NaN; rareEvent.bestVMin=NaN; rareEvent.predictionsRE=[]; rareEvent.predictionLevels=[];
    rareEvent.thresholddB= NaN; rareEvent.bestPaMindB=NaN;
    return
end

% find best fit rare event function
switch experiment.paradigm
    case 'gapDuration'
        % i.e. don't attempt this but visit the function to set null results
        rareEvent=fitRareEvent([], responsesPhaseTwo, stimulusParameters.targetDuration);
    otherwise
        rareEvent=fitRareEvent(levelsPhaseTwo, responsesPhaseTwo, stimulusParameters.targetDuration);
end

% find best logistic fit
logisticFunc=' 1./(1+exp(-a2.*(x-a1)));';
switch experiment.functionEstMethod
    % least squares estimate
    case {'logisticLS', 'rareEvent','peaksAndTroughs'}
        [a1, a2, Euclid]=fitFunctionUsingLeastSquares(levelsPhaseTwo, responsesPhaseTwo, ...
            logisticFunc, experiment.possLogSlopes, experiment.meanSearchStep);
%         [a1, a2, Euclid]=fitLogisticUsingLS(psy, levelsPhaseTwoBinVector); %! does not give same result
        logistic.bestThreshold=a1;
        logistic.bestK=a2;
        logistic.predictionsLOG= ...
            1./(1+exp(-logistic.bestK.*(experiment.predictionLevels-logistic.bestThreshold)));
        logistic.predictionLevels=experiment.predictionLevels;
        logistic.Euclid = Euclid;
        %         predResponses=1./(1+exp(-logistic.bestK.*(levelsPhaseTwo-logistic.bestThreshold)));
        %         [levelsPhaseTwo' responsesPhaseTwo' predResponses' responsesPhaseTwo'-predResponses' ]

        % maximum likelihood fitting
    case 'logisticML'
        [a1, a2, Euclid]=fitFunctionUsingMaxLikelihood (levelsPhaseTwoBinVector,...
            psy, yesResponses, noResponses, logisticFunc, experiment.possLogSlopes,experiment.meanSearchStep);
        logistic.bestThreshold=a1;
        logistic.bestK=a2;
        logistic.predictionsLOG= ...
            1./(1+exp(-logistic.bestK.*(experiment.predictionLevels-logistic.bestThreshold)));
        logistic.predictionLevels=experiment.predictionLevels;
        logistic.Euclid = Euclid;
end

% disp(num2str([logistic.bestThreshold logistic.bestK logistic.Euclid]))
% disp(num2str([ rareEvent.thresholddB rareEvent.bestGain rareEvent.bestVMin rareEvent.Euclid]))

% --------------------------------------------- fitLogisticUsingLS
function [mean, slope, Euclid]=fitLogisticUsingLS(p,L)

%!! each bin needs to be weighted by its size 

idx1=find(p>0);
p=p(idx1);
L=L(idx1);
idx1=find(p<1);
p=p(idx1);
L=L(idx1);

if length(L)<2
    mean=NaN;
    slope=NaN;
    Euclid=NaN;
    return
end

y=L;
x=log(1./p -1);
a =polyfit(x, y, 1);
mean=a(2);
slope=-1/a(1);

% euclid
predy=polyval(a,x);
Euclid=sum((predy-y).^2);
[mean slope Euclid];

% --------------------------------------------- fitFunctionUsingLeastSquares
function [a1, a2, Euclid]=fitFunctionUsingLeastSquares(levelsPhaseTwo, responsesPhaseTwo, func, possSlopes, meanSearchStep)

% nResponses= yesResponses+noResponses;
% idx=find(not(nResponses==0));
% levelsPhaseTwo=levelsPhaseTwo(idx);
% yesResponses=yesResponses(idx);
% noResponses=noResponses(idx);
% nResponses=nResponses(idx);

% psy=yesResponses./nResponses;

possThresholds=min(levelsPhaseTwo):meanSearchStep:max(levelsPhaseTwo);


% create vectors representing all combinations of a1 and a2
nPossThresh=length(possThresholds);
nPossSlopes=length(possSlopes);
a1=repmat(possThresholds',nPossSlopes,1);
a2=repmat(possSlopes,nPossThresh,1);
a2=reshape(a2,nPossThresh*nPossSlopes,1);

% each response is predicted by every possible threshold and slope
LS=ones(size(a1));
for i=1:length(levelsPhaseTwo)
    x=levelsPhaseTwo(i);
    eval(['y=' func]);
    actual=responsesPhaseTwo(i);
    error= (actual - y).^2;
    LS=LS+error;
end

[Euclid, idx]=min(LS);

a1=a1(idx);
a2=a2(idx);

% x=levelsPhaseTwo;
% eval(['y=' func]);
% Euclid= sum((psy-y).^2);

% plot(levelsPhaseTwo,y,'r')
% hold on
% plot(levelsPhaseTwo,psy,'o')

% --------------------------------------------- fitFunctionUsingMaxLikelihood
function [a1, a2, Euclid]=fitFunctionUsingMaxLikelihood...
    (levelsPhaseTwo,psy, yesResponses, noResponses, func, possible_a2, meanSearchStep)

% fitFunctionUsingMaxLikelihood fits the function in 'func' to binned yes-no data.
% levelsPhaseTwo specifies the bin centers
% yesResponses and noResponses are the niumber of yes and no responsesPhaseTwo in each bin
% 'func' is a function of a1, b1, x; e.g. func='1./(1+exp(-a2.*(x-a1)));';
%  and a2, a1 are parameters to be discovered.
% If bins are empty (i.e. neither yes or no), these are eliminated
%

nResponses=yesResponses+noResponses;
possible_a1=min(levelsPhaseTwo):meanSearchStep:max(levelsPhaseTwo);

% if nargin<6
%     possible_a2=-2:.05:4;
% end

% create vectors representing all combinations of a1 and a2
nPossThresh=length(possible_a1);
nPossSlopes=length(possible_a2);
a1=repmat(possible_a1',nPossSlopes,1);
a2=repmat(possible_a2,nPossThresh,1);
a2=reshape(a2,nPossThresh*nPossSlopes,1);

cumulativeProbability=ones(1, length(a2));
cumulativeProbability=cumulativeProbability';

for i=1:length(levelsPhaseTwo)
    x=levelsPhaseTwo(i);
    eval(['y=' func]);
    funcVal=(1-y).^yesResponses(i);
    cumulativeProbability= cumulativeProbability.*funcVal;
    funcVal=y.^noResponses(i);
    cumulativeProbability= cumulativeProbability.*funcVal;
end

[maxProb idx]=max(cumulativeProbability);
a1=a1(idx);
a2=a2(idx);

% x=levelsPhaseTwo;
eval(['y=' func]);
Euclid= sum(nResponses.*(psy-y).^2);

% figure(1), clf
% plot(levelsPhaseTwo,y)
% hold on
% plot(levelsPhaseTwo,psy,'o')



% --------------------------------------------------- fitRareEvent
function rareEvent=fitRareEvent(stimulusLevels, responsesPhaseTwo, duration, gains, Vmins)
% least squares estimate of *rare event* function
% model is: r = g P – A   ... g=events/s/Pa, A= events/s, P= pressure (Pa)
% psy=1- exp(d(gP –A  ))   ... d=duration
% p(event)=gain*levelmPa -Vmin
% 'responsesPhaseTwo' is a binary vector of subject's decision.
% 'stimulusLevels' are the corresponding signal levesl (values)
% duration is required to compute the expectation of an event occurring
% gains is an optional list of gains to be tried
% Vmins is an optional list of Vmins to be tried

global experiment
if nargin<5
    minVmin=.1; maxVmin=10; nVmins=100;
    Vmins=[0 logspace(log10(minVmin),log10(maxVmin),nVmins)];
    % disp(Vmins)
    gainMin=0.0001; gainMax=1; nGains=100; 
    gains=[1e-5 logspace(log10(gainMin),log10(gainMax),nGains)];
end

rareEvent.bestGain=NaN;
rareEvent.bestVMin=NaN;
rareEvent.thresholddB=0;
rareEvent.bestPaMindB=NaN;
rareEvent.predictionLevels=[];
rareEvent.predictionsRE=[];
rareEvent.Euclid=NaN;

if isempty(stimulusLevels), return, end

% expected slope is negative, rareEvent can not be used
if experiment.psyFunSlope<0
    return
end

% NB calculations in microPascals!
stimulusLevelsAsPressure=28 * 10.^(stimulusLevels/20);

% allGains=reshape(repmat(gains,nVmins,1), 1, nVmins*nGains);
% allVmins=repmat(Vmins, 1, nGains);

    predictions=NaN*zeros(1,length(stimulusLevels));
    gainCount=0;
    Euclid=inf; bestVmin=0; bestGain=0;
    for gain= gains
        gainCount=gainCount+1;
        VminCount=0;
        
        % 1 – exp(-d (g P – A))
        for Vmin=Vmins
            VminCount=VminCount+1;
            % all levelsPhaseTwo are simultaneously assessed
            
            % original linear function
%             gP_Vmin=gain*stimulusLevelsAsPressure-Vmin;
%             idx=(gP_Vmin>0);
            %   predictions(idx)= 1-exp(-duration*(gP_Vmin(idx)));
            %   new log function log(r) =gP-A
            
            % log function
            gP_Vmin=gain*stimulusLevelsAsPressure-Vmin;
            idx=(gP_Vmin>0);
            predictions(idx)= 1-exp(-duration*exp(gP_Vmin(idx)));
            predictions(~idx)=0;
            
%             % square function
%             P_Vmin=stimulusLevelsAsPressure-Vmin;
%             idx=(P_Vmin)>0;
%             predictions(idx)= 1-exp(-duration*gain*(P_Vmin(idx)).^2);
%             predictions(~idx)=0;
            
            
            % NB the error is equally weighted for each stimulus/response pair
            % this is not the same as equal weighting for each distinct level
            error=(predictions - responsesPhaseTwo).^2;
            error=mean(error(~isnan(error)));
            if error<Euclid
                Euclid=error;
                bestVmin=Vmin;
                bestVminCount=VminCount;
                %                 bestGainCount=gainCount;
                bestGain=gain;
            end
        end
    end
    % disp(Vmins)

% if bestGainCount==1 | bestGainCount==nGains
%     disp(['gain estimate ' num2str(gains(bestGainCount)) ' may be out of range'])
% end

[rareEvent.Euclid idx]=min(Euclid);
rareEvent.bestGain=bestGain;
rareEvent.bestVMin=bestVmin;
rareEvent.thresholdPa=(-log(0.5)/duration + rareEvent.bestVMin)/rareEvent.bestGain;
rareEvent.thresholddB=20*log10(rareEvent.thresholdPa/28);
rareEvent.bestPaMindB=20*log10((rareEvent.bestVMin/rareEvent.bestGain)/28);

predictionLevels= -50:1:120;
rareEvent.predictionLevels=predictionLevels;
stimulusLevelsAsPressure=28*10.^(predictionLevels/20);

% gP_Vmin=rareEvent.bestGain*stimulusLevelsAsPressure-rareEvent.bestVMin;
% rareEvent.predictionsRE= 1-exp(-duration*(gP_Vmin));

gP_Vmin=rareEvent.bestGain*stimulusLevelsAsPressure-rareEvent.bestVMin;
rareEvent.predictionsRE=  1-exp(-duration*exp(gP_Vmin));

%             P_Vmin=stimulusLevelsAsPressure-Vmin;
%             predictions= 1-exp(-duration*gain*P_Vmin.^2);
% predictions(predictions<0)=0;
%             rareEvent.predictionsRE=predictions;

% rareEvent

