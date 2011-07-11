function paradigm_SRT(handles)
global stimulusParameters experiment betweenRuns

stimulusParameters.subjectSampleRate=44100;

%  assessment method
% {'oneIntervalUpDown', 'MaxLikelihood', '2I2AFC++', '2I2AFC+++'}
experiment.threshEstMethod='oneIntervalUpDown';
% {'cued', 'noCue'};
stimulusParameters.includeCue=0;
stimulusParameters.cueTestDifference=10;

experiment.singleIntervalMaxTrials=10;
experiment.maxTrials=10;
experiment.allowCatchTrials= 0;

% {'tone','noise', 'pinkNoise','whiteNoise','OHIO'}
stimulusParameters.WRVname='targetLevel';
stimulusParameters.WRVstartValues=60 ;
stimulusParameters.WRVsteps=[5 2];
stimulusParameters.WRVlimits=[-30 110];

% target variable: slope=1, start going down.
experiment.psyFunSlope=1;
withinRuns.direction='down';

betweenRuns.variableName1='targetFrequency';
betweenRuns.variableList1=1000;
betweenRuns.variableName2='maskerDuration';
betweenRuns.variableList2=0.1 ;

% delay > masker > gap > target
stimulusParameters.stimulusDelay=0.3;

% maskerTypes={'tone','noise', 'pinkNoise','TEN','whiteNoise'};
experiment.maskerInUse=0;

stimulusParameters.gapDuration=0.0;

% targetTypes={'tone','noise', 'pinkNoise','whiteNoise','OHIO'};
stimulusParameters.targetType='digitStrings';
stimulusParameters.targetPhase='sin';
stimulusParameters.targetFrequency=1000;
stimulusParameters.targetDuration=2;
stimulusParameters.targetLevel=stimulusParameters.WRVstartValues(1);

stimulusParameters.rampDuration=0.004;
stimulusParameters.stimulusDelay=1;

% forced choice window interval
stimulusParameters.AFCsilenceDuration=0.5;

% {'none','noise', 'pinkNoise', 'TEN','noiseDich', 'pinkNoiseDich','whiteNoise'}
stimulusParameters.backgroundType='24TalkerBabble'; 
stimulusParameters.backgroundLevel= 60;

% instructions to user
%   single interval up/down no cue
stimulusParameters.instructions{1}= ...
    [{'Type three digits in the box (top left)'}, { }, ...
    {'then hit return'}];

stimulusParameters.numOHIOtones=1;

