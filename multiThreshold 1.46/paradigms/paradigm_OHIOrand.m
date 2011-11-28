function paradigm_OHIOrand(handles)
global stimulusParameters experiment betweenRuns

paradigmBase(handles) % default

betweenRuns.variableName1='OHIOnTones';
betweenRuns.variableList1=...
[2 4 8 12];
betweenRuns.variableName2='targetDuration';
betweenRuns.variableList2= 0.01;

experiment.maskerInUse=0;

stimulusParameters.OHIOnTones=betweenRuns.variableList1;
stimulusParameters.targetDuration=betweenRuns.variableList2;
stimulusParameters.targetLevel=stimulusParameters.WRVstartValues(1);


stimulusParameters.WRVstartValues=30;
experiment.singleIntervalMaxTrials=20;


% forced choice window interval
stimulusParameters.AFCsilenceDuration=0.5;


% instructions to user
%   single interval up/down no cue
stimulusParameters.instructions{1}=[{'YES if you hear the added click'}, { }, { 'NO if not (or you are uncertain'}];
%   single interval up/down with cue
stimulusParameters.instructions{2}=[{'count how many distinct clicks you hear'},{'ignore the tones'},{' '},...
    {'The clicks must be **clearly distinct** to count'}];

