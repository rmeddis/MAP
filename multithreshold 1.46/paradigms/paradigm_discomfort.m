function paradigm_discomfort(handles)
global stimulusParameters experiment betweenRuns

paradigmBase(handles) % default

stimulusParameters.WRVname='targetLevel';
stimulusParameters.WRVstartValues=75 ;
stimulusParameters.WRVsteps=[3 3];
stimulusParameters.WRVlimits=[-30 110];

betweenRuns.variableName1='targetFrequency';
betweenRuns.variableList1=[1000];
betweenRuns.variableName2='targetDuration';
betweenRuns.variableList2=0.5 ;
% 'randomize within blocks', 'fixed sequence', 'randomize across blocks'
betweenRuns.randomizeSequence='fixed sequence'; 

stimulusParameters.stimulusDelay=0;

experiment.maskerInUse=0;

stimulusParameters.targetType='tone';
stimulusParameters.targetPhase='sin';
stimulusParameters.targetFrequency=1000;
stimulusParameters.targetDuration=0.5;
stimulusParameters.targetLevel=stimulusParameters.WRVstartValues(1);

stimulusParameters.instructions{1}= ['Is the tone ''comfortable'', ''loud'' or ''uncomfortable''?'];
%   single interval up/down with cue
stimulusParameters.instructions{2}= [];

