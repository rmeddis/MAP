function paradigm_TMC(handles)
global stimulusParameters experiment betweenRuns

paradigmBase(handles) % default

stimulusParameters.WRVname='maskerLevel';
stimulusParameters.WRVstartValues=50;
stimulusParameters.WRVsteps= [-10 -4];
stimulusParameters.WRVlimits=[-30 110];


stimulusParameters.cueTestDifference = 10;
experiment.psyFunSlope = -1;
withinRuns.direction='up';

betweenRuns.variableName1='gapDuration';
betweenRuns.variableList1=[.05 .09 .01 .07 .03 ];
betweenRuns.variableName2='targetFrequency';
% retain existing targetFrequencies
betweenRuns.variableList2=str2num(get(handles.edittargetFrequency,'string'));

experiment.maskerInUse=1;
stimulusParameters.maskerType='tone';
stimulusParameters.maskerPhase='sin';
stimulusParameters.maskerDuration=0.108;
stimulusParameters.maskerLevel=stimulusParameters.WRVstartValues(1);
stimulusParameters.maskerRelativeFrequency=1;

stimulusParameters.gapDuration=betweenRuns.variableList1;

stimulusParameters.targetType='tone';
stimulusParameters.targetPhase='sin';
stimulusParameters.targetFrequency=betweenRuns.variableList2(1);
stimulusParameters.targetDuration=0.016;
stimulusParameters.targetLevel=NaN;

stimulusParameters.rampDuration=0.004;

% instructions to user
%   single interval up/down no cue
stimulusParameters.instructions{1}=[{'YES if you hear the added click'}, { }, { 'NO if not (or you are uncertain'}];
%   single interval up/down with cue
stimulusParameters.instructions{2}=[{'count how many distinct clicks you hear'},{'ignore the tones'},{' '},...
    {'The clicks must be **clearly distinct** to count'}];

