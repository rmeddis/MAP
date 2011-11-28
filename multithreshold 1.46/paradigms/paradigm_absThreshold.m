function paradigm_absThreshold(handles)
global stimulusParameters experiment betweenRuns

paradigmBase(handles) % default

betweenRuns.variableName1='targetFrequency';
betweenRuns.variableList1=[250 500 1000 2000 4000 8000 ];
betweenRuns.variableName2='targetDuration';
betweenRuns.variableList2= 0.25;

experiment.maskerInUse=0;

stimulusParameters.targetFrequency=betweenRuns.variableList1;
stimulusParameters.targetDuration=betweenRuns.variableList2;
stimulusParameters.targetLevel=stimulusParameters.WRVstartValues(1);

stimulusParameters.WRVstartValues=30;


% forced choice window interval
stimulusParameters.AFCsilenceDuration=0.5;


% instructions to user
%   single interval up/down no cue
stimulusParameters.instructions{1}=[{'YES if you hear the added click'}, { }, { 'NO if not (or you are uncertain'}];
%   single interval up/down with cue
stimulusParameters.instructions{2}=[{'count how many distinct clicks you hear'},{'ignore the tones'},{' '},...
    {'The clicks must be **clearly distinct** to count'}];


