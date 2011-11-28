function paradigm_threshold_duration(handles)
global stimulusParameters experiment betweenRuns

paradigmBase(handles) % default

stimulusParameters.WRVname='targetLevel';
stimulusParameters.WRVstartValues=40;
stimulusParameters.WRVsteps=[10 2];
stimulusParameters.WRVlimits=[-30 110];

betweenRuns.variableName1='targetDuration';
betweenRuns.variableList1=[ .016 .032 .064 .128 .256 .512];
betweenRuns.variableName2='targetFrequency';
betweenRuns.variableList2=1000;

experiment.maskerInUse=0;

stimulusParameters.targetType='tone';
stimulusParameters.targetPhase='sin';
% 		retain current target frequency
x=str2num(get(handles.edittargetFrequency,'string'));
stimulusParameters.targetFrequency=x(1);
stimulusParameters.targetDuration=betweenRuns.variableList2;
stimulusParameters.targetLevel=stimulusParameters.WRVstartValues(1);

stimulusParameters.rampDuration=0.004;

% instructions to user
%   single interval up/down no cue
stimulusParameters.instructions{1}=[{'YES if you hear the added click'}, { }, { 'NO if not (or you are uncertain'}];
%   single interval up/down with cue
stimulusParameters.instructions{2}=[{'count how many distinct clicks you hear'},{'ignore the tones'},{' '},...
    {'The clicks must be **clearly distinct** to count'}];

