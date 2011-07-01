function [nextStep, msg]=Levitt2 (decision, currentValue)
global LevittControl withinRuns
% descision is a string: 'hit' or 'miss'
% current value is the value of the variable parameter (e.g. level)
% msg:
%     msg='' indicates trial progressing successfully
%     msg='done' indicates that a threshold estimate has been found
%     msg='maximum level exceeded'
%     msg='maximum no of trials exceeded'

% Initialize by running Levitt2 without arguments
%  later LevittControl will be set in aReadAndCheckParameterBoxes in  expGUI_MT
if nargin==0
    
    LevittControl.peakTroughValues=[];   
    LevittControl.trialRunning=0;
    LevittControl.sequence='****'; 
    LevittControl.LevittValuesUsed=[NaN NaN NaN NaN];
    LevittControl.TurnsToSmallSteps=4; % 2 peaks and 2 troughs
    LevittControl.direction='same';
    LevittControl.prevDirection='easier';
    LevittControl.Nreversals=0;
    
    LevittControl.trialRunning=1;
    LevittControl.meanPeakTrough=NaN;
    LevittControl.sd=NaN;
    msg='';
    return
end

% apply Levitt rules to find next stimulus value
rule=LevittControl.rule;
sequence=LevittControl.sequence;
meanPT=LevittControl.meanPeakTrough;
sdPT=LevittControl.sd;

% response sequence: '+' is a hit and '-' is a miss.
if strcmp(decision,'hit')
    sequence=[sequence '+'];
else
    sequence=[sequence '-'];
end
LevittControl.LevittValuesUsed=[LevittControl.LevittValuesUsed withinRuns.levelList(end)];

switch rule
    case '++'
        peakCriterion='-++';
        troughCriteria={'++ -', '++ +-'};
    case '+++'
        peakCriterion='-+++';
        troughCriteria={'+++ -', '+++ +-', '+++ ++-'};
    otherwise
        error('Levitt:')
end	

troughs=[]; allTroughPtrs=[];
for i=1:length(troughCriteria)              % do all trough criteria
    troughCriterion=char(troughCriteria(i));% one criterion at a time
    % identify the location of a trough
    troughPtrs=findstr(sequence, troughCriterion) + length(troughCriterion)-1;
    % identify the level at which it occurred
    troughLevels=LevittControl.LevittValuesUsed(troughPtrs);
    % archive the list
    withinRuns.troughs=troughLevels;
    troughs=[troughs troughLevels];
    allTroughPtrs=[allTroughPtrs troughPtrs];
end
% only one peak criterion used
withinRuns.troughs=troughs;
peakPtrs=findstr(sequence,peakCriterion)+length(peakCriterion)-1;
peakLevels=LevittControl.LevittValuesUsed(peakPtrs);
withinRuns.peaks=peakLevels;

% almagamate and sort into date order
peakTroughList=[peakLevels troughs];
peakTroughPtrs=[peakPtrs allTroughPtrs];
[peakTroughPtrs idx]=sort(peakTroughPtrs);
% it needs to be sequenced so that the algorithm can take the last peaks
% and troughs
peakTroughList=peakTroughList(idx);

% adjust step size as the trial progresses
% a positive step size indicates the 'harder' direction 
if length(peakTroughList)>=LevittControl.TurnsToSmallSteps
    currentStep=LevittControl.steadyLevittStep;
else
    currentStep=LevittControl.startLevelStep;
end    

% base next stimulus on the basis of the sequence
%  any miss requires an 'easier' stimulus next time.
if strcmp(sequence(end),'-')
    nextStep= -currentStep;
    
    % success requires 2 or more successive hits
elseif strcmp(sequence(end-length(rule)+1:end),rule)
    sequence=[sequence ' '];  % add space to prevent success on successive trials
    LevittControl.LevittValuesUsed=[LevittControl.LevittValuesUsed NaN];
    nextStep=currentStep;
    
    % not enough hits to provoke a change
else	
    nextStep=0;
end

LevittControl.sequence=sequence;

LevittControl.Nreversals=length(peakTroughList);
% compute threshold estimate 
% only if minReversals exceeded and even number of peaks and troughs
if LevittControl.Nreversals>=LevittControl.minReversals && rem(length(peakTroughList),2)==0   
    % use only the peaks and troughs at the end of the sequence
    peaksAndTroughs=peakTroughList(end-LevittControl.useLastNturns+1:end);
    disp(['peak/trough sequence= ' num2str(peaksAndTroughs, '%6.0f') ] )
    
    meanPT=mean(peaksAndTroughs);
    sdPT=std(peaksAndTroughs);
    LevittControl.meanPeakTrough=meanPT;
    LevittControl.sd=sdPT;
    fprintf('Levitt, mean, sd= %6.1f,%6.1f\n', meanPT, sdPT)
    % final check that the sd is low enough
    if sdPT<LevittControl.targetsdPT        
        nextValue=currentValue;
        msg='done';
        return
    end
end

nextValue=currentValue + nextStep;
if nextValue>LevittControl.maxLevittValue
    msg='maximum level exceeded'
    return
end

% required for later use
LevittControl.trialRunning=LevittControl.trialRunning+1;
LevittControl.peakTroughValues=peakTroughList;

if LevittControl.trialRunning>LevittControl.maxTrials
    msg='maximum no of trials exceeded'
    return
end

% Trial continues
msg='';
