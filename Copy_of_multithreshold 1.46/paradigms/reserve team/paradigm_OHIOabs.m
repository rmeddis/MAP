function paradigm_OHIOabs(handles)
global stimulusParameters experiment betweenRuns

paradigm_training(handles) % default

% find the  threshold for a tonecomplex  consisting of a sequence of 10-ms tones
%   whose frequencies are chosen at random from a list (OHIOfrequencies)
% All tones are presented at the same level (SL) computed using absolute
%   threshols specified in OHIOthresholds;
% The duration of the complex is increased across runs and the number of tones is
%   controlled by OHIOdurations, (for each 20 ms a further tone is added.
% The frequency of the tones is changed on each trial

experiment.OHIOfrequencies=[494, 663, 870, 1125, 1442, 1838, 2338, 2957, 3725, 4680, 5866,  7334]; %Hz.
% User must specify abs thresholds (dB SPL) of each tone frequency
% experiment.OHIOthresholds= [18	16	16	19	20	22	24	26	27	30	32	35];

%  assessment method
% {'oneIntervalUpDown', 'MaxLikelihood', '2I2AFC++', '2I2AFC+++'}
experiment.threshEstMethod='oneIntervalUpDown';
% {'cued', 'noCue'};

stimulusParameters.WRVname='targetLevel';
stimulusParameters.WRVstartValues=20 ;
stimulusParameters.WRVsteps=[10 2];
stimulusParameters.WRVlimits=[-30 110];

betweenRuns.variableName1='numOHIOtones';
betweenRuns.variableList1= 1:12; % i.e. the frequency to be used
betweenRuns.variableName2='stimulusDelay';
betweenRuns.variableList2=0.05;
betweenRuns.randomizeSequence=2; % not random sequence

stimulusParameters.targetType='OHIO';
stimulusParameters.targetPhase='sin';
stimulusParameters.targetFrequency=experiment.OHIOfrequencies;
stimulusParameters.targetDuration=0.01; % overruled by OHIO program
stimulusParameters.targetLevel=stimulusParameters.WRVstartValues(1);

stimulusParameters.rampDuration=0.005;

% instructions to user
%   single interval up/down no cue
stimulusParameters.instructions{1}= [{'YES if you hear the tone clearly'}, { }, { 'NO if not (or you are uncertain'}];
%   single interval up/down with cue
stimulusParameters.instructions{2}= [{'count the tones you hear clearly'}, { }, { 'ignore indistinct tones'}];

