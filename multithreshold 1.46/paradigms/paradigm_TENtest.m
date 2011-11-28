function paradigm_TENtest(handles)
global stimulusParameters experiment betweenRuns

paradigmBase(handles) % default

stimulusParameters.WRVname='targetLevel';
stimulusParameters.WRVstartValues=40;
stimulusParameters.WRVsteps=[10 2];
stimulusParameters.WRVlimits=[-30 110];

betweenRuns.variableName1='targetFrequency';
betweenRuns.variableList1=[250 500 1000 2000 4000 8000 ];
betweenRuns.variableName2='targetDuration';
betweenRuns.variableList2= 0.25;

experiment.maskerInUse=0;

stimulusParameters.targetType='tone';
stimulusParameters.targetPhase='sin';
stimulusParameters.targetFrequency=betweenRuns.variableList1;
stimulusParameters.targetDuration=0.5;
stimulusParameters.targetLevel=stimulusParameters.WRVstartValues(1);

% {'none','noise', 'pinkNoise', 'TEN','noiseDich', 'pinkNoiseDich','whiteNoise'}
stimulusParameters.backgroundType='TEN'; 
stimulusParameters.backgroundLevel=20;

% instructions to user
%   single interval up/down no cue
stimulusParameters.instructions{1}=[{'YES if you hear the added click'}, { }, { 'NO if not (or you are uncertain'}];
%   single interval up/down with cue
stimulusParameters.instructions{2}=[{'count how many distinct clicks you hear'},{'ignore the noise'},{' '},...
    {'The clicks must be **clearly distinct** to count'}];

