function paradigm_psychometric(handles)
global stimulusParameters experiment betweenRuns

paradigmBase(handles) % default

experiment.printTracks=1;
experiment.maxTrials=30;

stimulusParameters.WRVname='targetLevel';
stimulusParameters.WRVstartValues=30 ;
stimulusParameters.WRVsteps=[10 2];
stimulusParameters.WRVlimits=[-30 110];

betweenRuns.variableName1='targetFrequency';
betweenRuns.variableList1=1000.01:0.01:1000.05;
betweenRuns.variableName2='targetDuration';
betweenRuns.variableList2=0.1 ;

experiment.maskerInUse=0;

stimulusParameters.targetType='tone';
stimulusParameters.targetPhase='sin';
stimulusParameters.targetFrequency=1000;
stimulusParameters.targetDuration=0.1;
stimulusParameters.targetLevel=stimulusParameters.WRVstartValues(1);

% instructions to user
%   single interval up/down no cue
stimulusParameters.instructions{1}= [{'YES if you hear the tone clearly'}, { }, { 'NO if not (or you are uncertain'}];
%   single interval up/down with cue
stimulusParameters.instructions{2}= [{'count the tones you hear clearly'}, { }, { 'ignore indistinct tones'}];


