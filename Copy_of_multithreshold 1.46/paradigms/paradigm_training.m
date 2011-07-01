function paradigm_training(handles)
global stimulusParameters experiment betweenRuns

stimulusParameters.subjectSampleRate=44100; % compatible with file input
% stimulusParameters.subjectSampleRate=64000; % compatible with file input
% stimulusParameters.subjectSampleRate=128000; % compatible with file input

%  assessment method
% {'oneIntervalUpDown', 'MaxLikelihood', '2I2AFC++', '2I2AFC+++'}
experiment.threshEstMethod='oneIntervalUpDown';
% {'cued', 'noCue'};
stimulusParameters.includeCue=1;
stimulusParameters.cueTestDifference=10;

experiment.singleIntervalMaxTrials=10;
experiment.maxTrials=10;
experiment.allowCatchTrials= 1;

% {'tone','noise', 'pinkNoise','whiteNoise','OHIO'}
stimulusParameters.WRVname='targetLevel';
stimulusParameters.WRVstartValues=20 ;
stimulusParameters.WRVsteps=[10 2];
stimulusParameters.WRVlimits=[-30 110];

% target variable: slope=1, start going down.
experiment.psyFunSlope=1;
withinRuns.direction='down';

betweenRuns.variableName1='targetFrequency';
betweenRuns.variableList1=1000;
betweenRuns.variableName2='targetDuration';
betweenRuns.variableList2=0.1 ;
% 1='fixed sequence', 2='randomize within blocks', 3='randomize across
% blocks'
betweenRuns.randomizeSequence=1; % 'random sequence'

% delay > masker > gap > target

stimulusParameters.stimulusDelay=0.3;

% maskerTypes={'tone','noise', 'pinkNoise','TEN','whiteNoise'};
stimulusParameters.maskerType='tone';
stimulusParameters.maskerPhase='sin';
stimulusParameters.maskerDuration=0.0;
stimulusParameters.maskerLevel= -50;
stimulusParameters.maskerRelativeFrequency= 1 ;  % not used

stimulusParameters.gapDuration=0.0;

% targetTypes={'tone','noise', 'pinkNoise','whiteNoise','OHIO'};
stimulusParameters.targetType='tone';
stimulusParameters.targetPhase='sin';
stimulusParameters.targetFrequency=1000;
stimulusParameters.targetDuration=0.1;
stimulusParameters.targetLevel=stimulusParameters.WRVstartValues(1);

stimulusParameters.rampDuration=0.004;

% forced choice window interval
stimulusParameters.AFCsilenceDuration=0.5;

% {'none','noise', 'pinkNoise', 'TEN','noiseDich', 'pinkNoiseDich','whiteNoise'}
stimulusParameters.backgroundType='none'; 
stimulusParameters.backgroundLevel=-100;

% instructions to user
%   single interval up/down no cue
stimulusParameters.instructions{1}= [{'YES if you hear the tone clearly'}, { }, { 'NO if not (or you are uncertain'}];
%   single interval up/down with cue
stimulusParameters.instructions{2}= [{'count the tones you hear clearly'}, { }, { 'ignore indistinct tones'}];

stimulusParameters.numOHIOtones=1;

