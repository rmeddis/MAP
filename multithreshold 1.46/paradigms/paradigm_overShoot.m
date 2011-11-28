function paradigm_overShoot(handles)
global stimulusParameters  betweenRuns experiment

paradigmBase(handles) % default

stimulusParameters.WRVname='targetLevel';
stimulusParameters.WRVstartValues=50;
stimulusParameters.WRVsteps= [10 2];
stimulusParameters.WRVlimits=[-30 110];

betweenRuns.variableName1='gapDuration';
betweenRuns.variableList1=[-.399 -.2];
% betweenRuns.variableList1=[-.288 -.218 -.2 -.190 -.170 -.130 -.070 -.030 -.010 .005 .020 .080];
% betweenRuns.variableList1=[-.350 -.238 -.213 -.180 -.160 -.100 -.040 -.020 0 .010 .040 .140];
betweenRuns.variableName2='maskerLevel';
betweenRuns.variableList2=50;

% delay > masker > gap > target
stimulusParameters.stimulusDelay=0.3;

experiment.maskerInUse=1;
stimulusParameters.maskerType='tone';
stimulusParameters.maskerPhase='sin';
stimulusParameters.maskerDuration=0.4;
stimulusParameters.maskerLevel=betweenRuns.variableList2;
stimulusParameters.maskerRelativeFrequency=1;

stimulusParameters.gapDuration=betweenRuns.variableList1;

stimulusParameters.targetType='tone';
stimulusParameters.targetPhase='sin';
stimulusParameters.targetFrequency=1000;
stimulusParameters.targetDuration=0.01;
stimulusParameters.targetLevel=-stimulusParameters.WRVstartValues(1);

stimulusParameters.rampDuration=0.005;

% instructions to user
% single interval up/down no cue
stimulusParameters.instructions{1}=[{'YES if you hear the added click'}, { }, { 'NO if not (or you are uncertain'}];
% single interval up/down with cue
stimulusParameters.instructions{2}=[{'count how many distinct clicks you hear'},{'ignore the tones'},{' '},...
    {'The clicks must be **clearly distinct** to count'}];

