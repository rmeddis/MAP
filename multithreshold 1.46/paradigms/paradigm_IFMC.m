function paradigm_IFMC(handles)
global stimulusParameters experiment betweenRuns

paradigm_training(handles) % default

stimulusParameters.WRVname='maskerLevel';
stimulusParameters.WRVstartValues=10;
stimulusParameters.WRVsteps= [-10 -2];
stimulusParameters.WRVlimits=[-30 110];

experiment.psyFunSlope = -1;
withinRuns.direction='up';

betweenRuns.variableName1='maskerRelativeFrequency';
betweenRuns.variableList1=[1       0.5       1.6   .9 .7   1.1 1.3      ];
betweenRuns.variableName2='targetFrequency';
% keep old list of target frequencies
betweenRuns.variableList2=str2num(get(handles.edittargetFrequency,'string'));

experiment.maskerInUse=1;
stimulusParameters.maskerType='tone';
stimulusParameters.maskerPhase='sin';
stimulusParameters.maskerDuration=0.108;
stimulusParameters.maskerLevel=stimulusParameters.WRVstartValues(1);
stimulusParameters.maskerRelativeFrequency=betweenRuns.variableList1;

stimulusParameters.gapDuration=0.01;

stimulusParameters.targetType='tone';
stimulusParameters.targetPhase='sin';
stimulusParameters.targetFrequency=betweenRuns.variableList2(1);
stimulusParameters.targetDuration=0.016;
stimulusParameters.targetLevel=NaN;

stimulusParameters.rampDuration=0.004;

% instructions to user
% single interval up/down no cue
stimulusParameters.instructions{1}=[{'YES if you hear the added click'}, { }, { 'NO if not (or you are uncertain'}];
% single interval up/down with cue
stimulusParameters.instructions{2}=[{'count how many distinct clicks you hear'},{'ignore the tones'},{' '},...
    {'The clicks must be **clearly distinct** to count'}];

