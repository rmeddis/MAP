function paradigm_OHIOtemp(handles)
global stimulusParameters experiment betweenRuns

paradigm_training(handles) % default

% find the  threshold for a tonecomplex  consisting of a sequence of 10-ms tones
%   whose frequencies are chosen at random from a list (OHIOfrequencies)
% All tones are presented at the same level (SL) computed using absolute
%   threshols specified in OHIOthresholds;
% The duration of the complex is increased across runs and the number of tones is
%   controlled by OHIOdurations, (for each 20 ms a further tone is added.
% The frequency of the tones is changed on each trial

% fetch thresholds and frequencies
experiment=OHIOthresholds(experiment);

stimulusParameters.WRVname='targetLevel';
stimulusParameters.WRVstartValues=0 ;
stimulusParameters.WRVsteps=[10 2];
stimulusParameters.WRVlimits=[-30 110];
% target variable: slope=1, start going down.
stimulusParameters.cueTestDifference=10;
experiment.psyFunSlope= 1;
withinRuns.direction='down';

betweenRuns.variableName1='numOHIOtones';
betweenRuns.variableList1= [1 2 4 8 12];
betweenRuns.variableName2='stimulusDelay';
betweenRuns.variableList2=0.1;
betweenRuns.randomizeSequence=2; % not random sequence

stimulusParameters.targetType='OHIO';
stimulusParameters.targetPhase='sin';
stimulusParameters.targetFrequency=experiment.OHIOfrequencies;
stimulusParameters.targetDuration=betweenRuns.variableList2;
stimulusParameters.targetLevel=stimulusParameters.WRVstartValues(1);

stimulusParameters.rampDuration=0.005;

% instructions to user
%   single interval up/down no cue
stimulusParameters.instructions{1}= [{'YES if you hear the tone clearly'}, { }, { 'NO if not (or you are uncertain'}];
%   single interval up/down with cue
stimulusParameters.instructions{2}= [{'count the tones you hear clearly'}, { }, { 'ignore indistinct tones'}];
