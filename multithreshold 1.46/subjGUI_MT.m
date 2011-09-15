function varargout = subjGUI_MT(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @subjGUI_MT_OpeningFcn, ...
    'gui_OutputFcn',  @subjGUI_MT_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before subjGUI_MT is made visible.
function subjGUI_MT_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for subjGUI_MT
handles.output = hObject;
initializeGUI(handles)
guidata(hObject, handles);

function varargout = subjGUI_MT_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

% -----------------------------------------------------initializeGUI
function initializeGUI(handles)
global experiment
global subjectGUIHandles expGUIhandles
addpath (['..' filesep 'MAP'], ['..' filesep 'utilities'], ...
    ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'],...
    ['..' filesep 'testPrograms'])

dbstop if error

% subjectGUI size and location         % [left bottom width height]
scrnsize=get(0,'screensize');
set(0, 'units','pixels')
switch experiment.ear
    % use default size unless...
    case {'MAPmodel',  'MAPmodelMultich', 'MAPmodelSingleCh', ...
            'statsModelLogistic','statsModelRareEvent'}
        % 	subjectGUI not needed for modelling so minimize subject GUI
        set(gcf, 'units','pixels')
        y=[0*scrnsize(3) 0.8*scrnsize(4) 0.1*scrnsize(3) 0.2*scrnsize(4)];
        set(gcf,'position',y, 'color',[.871 .961 .996])

    case 'MAPmodelListen',
        % 	subjectGUI is needed for display purposes. Make it large
        set(gcf, 'units','pixels')
        y=[.665*scrnsize(3) 0.02*scrnsize(4) ...
            0.33*scrnsize(3) 0.5*scrnsize(4)]; % alongside
        set(gcf,'position',y, 'color',[.871 .961 .996])
end

switch experiment.ear
    case{'left', 'right','diotic', 'dichoticLeft','dichoticRight'}
        % Look to see if the button box exists and, if so, initialise it
        buttonBoxIntitialize		% harmless if no button box attached
end

% clear display of previous mean values. This is a new measurement series
axes(expGUIhandles.axes4), cla
reset (expGUIhandles.axes4)

% handles needed in non-callback routines below
subjectGUIHandles=handles;

% start immediately
startNewExperiment(handles, expGUIhandles)
% This is the end of the experiment. Exit here and return to ExpGUI.

% ----------------------------------------------------- startNewExperiment
function startNewExperiment(handles, expGUIhandles)
% An experiment consists of a series of 'runs'.
% Resets all relevant variables at the beginning of a new experiment.
global experiment stimulusParameters betweenRuns

% 'start new experiment' button is the only valid action now
experiment.status='waitingForStart';

switch experiment.threshEstMethod
    % add appropriate labels to subject GUI buttons
    case {'2I2AFC++', '2I2AFC+++'}
        set(handles.pushbutton3,'string','')
        set(handles.pushbutton2,'string','2')
        set(handles.pushbutton1,'string','1')
        set(handles.pushbutton0,'string','0')
    case {'MaxLikelihood', 'oneIntervalUpDown'}
        if stimulusParameters.includeCue
            set(handles.pushbutton3,'string','')
            set(handles.pushbutton2,'string','2')
            set(handles.pushbutton1,'string','1')
            set(handles.pushbutton0,'string','0')
        else
            set(handles.pushbutton3,'string','')
            set(handles.pushbutton2,'string','YES')
            set(handles.pushbutton1,'string','NO')
            set(handles.pushbutton0,'string','')
        end
end

switch experiment.paradigm
    case 'discomfort'
        set(handles.pushbutton3,'string','')
        set(handles.pushbutton2,'string','uncomfortable')
        set(handles.pushbutton1,'string','loud')
        set(handles.pushbutton0,'string','comfortable')
        experiment.allowCatchTrials=0;
end

% experiment.subjGUIfontSize is set on expGUI
set(handles.pushbutton3,'FontSize',experiment.subjGUIfontSize)
set(handles.pushbutton2,'FontSize',experiment.subjGUIfontSize)
set(handles.pushbutton1,'FontSize',experiment.subjGUIfontSize)
set(handles.pushbutton0,'FontSize',experiment.subjGUIfontSize)
set(handles.pushbuttoNotSure,'FontSize',experiment.subjGUIfontSize)
set(handles.pushbuttonGO,'FontSize',experiment.subjGUIfontSize)
set(handles.textMSG,'FontSize',experiment.subjGUIfontSize)

set(handles.pushbutton19,'visible','off') % unused button

% start command window summary of progress
fprintf(' \n ----------- NEW MEASUREMENTS\n')
disp(['paradigm:   ' experiment.paradigm])
cla(expGUIhandles.axes1)
cla(expGUIhandles.axes2)
cla(expGUIhandles.axes4)
cla(expGUIhandles.axes5)

experiment.stop=0;              % status of 'stop' button
experiment.pleaseRepeat=0;      % status of 'repeat' button
experiment.buttonBoxStatus='not busy';

% date and time and replace ':' with '_'
date=datestr(now);idx=findstr(':',date);date(idx)='_';
experiment.date=date;
timeNow=clock; betweenRuns.timeNow= timeNow;
experiment.timeAtStart=[num2str(timeNow(4)) ':' num2str(timeNow(5))];
experiment.minElapsed=0;

% unpack catch trial rates. The rate declines from the start rate
%  to the base rate using a time constant.
stimulusParameters.catchTrialRate=stimulusParameters.catchTrialRates(1);
stimulusParameters.catchTrialBaseRate=...
    stimulusParameters.catchTrialRates(2);
stimulusParameters.catchTrialTimeConstant=...
    stimulusParameters.catchTrialRates(3);
if stimulusParameters.catchTrialBaseRate==0
    stimulusParameters.catchTrialRate=0;
end

% for human measurements only, identify the start catch trial rate
switch experiment.ear
    case{'left', 'right','diotic', 'dichoticLeft','dichoticRight'}
        fprintf('stimulusParameters.catchTrialRate= %6.3f\n', ...
            stimulusParameters.catchTrialRate)
end

% Reset betweenRuns parameters. this occurs only at experiment start
% withinRuns values are reset in 'startNewRun'
% this approach creates a more readable structure summary printout.
betweenRuns.thresholds=[];
betweenRuns.thresholds_mean=[];
betweenRuns.thresholds_median=[];
betweenRuns.forceThresholds=[];
betweenRuns.observationCount=[];
betweenRuns.catchTrials=[];
betweenRuns.timesOfFirstReversals=[];
betweenRuns.bestThresholdTracks=[];
betweenRuns.levelTracks=[];
betweenRuns.responseTracks=[];
betweenRuns.slopeKTracks=[];
betweenRuns.gainTracks=[];
betweenRuns.VminTracks=[];
betweenRuns.bestGain=[];
betweenRuns.bestVMin=[];
betweenRuns.bestPaMin=[];
betweenRuns.bestLogisticM=[];
betweenRuns.bestLogisticK=[];
betweenRuns.resets=0;
betweenRuns.runNumber=0;

% Up to two parameters can be changed between runs
% Find the variable parameters and randomize them
% e.g. 'variableList1 = stimulusParameters.targetFrequency;'
eval(['variableList1=stimulusParameters.' betweenRuns.variableName1 ';']);
eval(['variableList2=stimulusParameters.' betweenRuns.variableName2 ';']);
nVar1=length(variableList1);
nVar2=length(variableList2);

% Create two sequence vectors to represent the sequence of var1 and var2
% values. 'var1' changes most rapidly.
switch betweenRuns.randomizeSequence
    % {'randomize within blocks', 'fixed sequence',...
    %  'randomize across blocks'}
    case 'fixed sequence'
        var1Sequence=repmat(betweenRuns.variableList1, 1,nVar2);
        var2Sequence=reshape(repmat(betweenRuns.variableList2, ...
            nVar1,1),1,nVar1*nVar2);
    case 'randomize within blocks'
        % the blocks are not randomized
        var1Sequence=betweenRuns.variableList1;
        ranNums=rand(1, length(var1Sequence)); [x idx]=sort(ranNums);
        var1Sequence=var1Sequence(idx);
        betweenRuns.variableList1=variableList1(idx);
        var1Sequence=repmat(var1Sequence, 1,nVar2);
        var2Sequence=reshape(repmat(betweenRuns.variableList2, nVar1,1)...
            ,1,nVar1*nVar2);
    case 'randomize across blocks'
        var1Sequence=repmat(betweenRuns.variableList1, 1,nVar2);
        var2Sequence=reshape(repmat(betweenRuns.variableList2, nVar1,1),...
            1,nVar1*nVar2);
        ranNums=rand(1, nVar1*nVar2);
        [x idx]=sort(ranNums);
        var1Sequence=var1Sequence(idx);
        var2Sequence=var2Sequence(idx);
        % there should be one start value for every combination
        %  of var1/ var2. In principle this allows these values to be
        % programmed. Not currently in use.
        stimulusParameters.WRVstartValues=...
            stimulusParameters.WRVstartValues(idx);
end
betweenRuns.var1Sequence=var1Sequence;
betweenRuns.var2Sequence=var2Sequence;

% caught out vector needs to be linked to the length of the whole sequence
betweenRuns.caughtOut=zeros(1,length(var1Sequence));

disp('planned sequence:')
if min(var1Sequence)>1
    % use decidaml places only if necessary
    disp([betweenRuns.variableName1 ': ' num2str(var1Sequence,'%6.0f')  ])
else
    disp([betweenRuns.variableName1 ': ' num2str(var1Sequence,'%8.3f')  ])
end
if min(var1Sequence)>1
    disp([betweenRuns.variableName2 ': ' num2str(var2Sequence,'%6.0f') ])
else
    disp([betweenRuns.variableName2 ': ' num2str(var2Sequence,'%8.3f') ])
end

fprintf('\nvariable1 \t  variable2\t  \n')
fprintf('%s \t  %s\t  Threshold  \n',betweenRuns.variableName1,...
    betweenRuns.variableName2)

% Light up 'GO' on subjGUI and advise.
set(handles.editdigitInput,'visible','off')
switch experiment.ear
    case {'statsModelLogistic', 'statsModelRareEvent',...
            'MAPmodel',  'MAPmodelMultiCh','MAPmodelSingleCh'}
        % no changes required if model used
    otherwise
        set(handles.pushbuttonGO,'backgroundcolor','y')
        set(handles.pushbuttonGO,'visible','on')
        set(handles.frame1,'visible','off')
        set(handles.textMSG,'backgroundcolor', 'w')
        msg=[{'Ready to start new Experiment'}, {' '}, {'Please, click on the GO button'}];
        set(handles.textMSG,'string', msg)

        set(handles.pushbuttoNotSure,'visible','off')
        set(handles.pushbuttonWrongButton,'visible','off')
        set(handles.pushbutton3,'visible','off')
        set(handles.pushbutton2,'visible','off')
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton0,'visible','off')
        pause(.1) % to allow display to be drawn
end

% Selecting the 'GO' button is the only valid operation action now
experiment.status='waitingForGO'; 	% i.e. waiting for new run

% control is now either manual, model (MAP) or randomization
switch experiment.ear
    case {'MAPmodel','MAPmodelMultiCh','MAPmodelSingleCh','MAPmodelListen'}                                     % MAP model is now the subject
        stimulusParameters.calibrationdB=0;             % Pascals required!
        MAPmodelRunsGUI(handles)
        % model is now the subject
    case  {'statsModelLogistic', 'statsModelRareEvent'}
        % no catch trials for the statistical model
        stimulusParameters.catchTrialBaseRate=0;
        stimulusParameters.catchTrialRate=0;
        statsModelRunsGUI(handles)
    otherwise
        %manual operation; wait for user to click on 'GO'
end

% Experiment complete (after MAP or randomization)
% return to  'initializeGUI' and then back to expGUI
% Manual control finds its own way home. Program control assumed when
% the user hits the GO button

% -----------------------------------------------------------------   startNewRun
function startNewRun(handles)
% There are many ways to arrive here.
%  Under manual control this is achieved by hitting the GO button
%   either via the button box or a mouse click
%  MAP and randomization methods call this too

global experiment stimulusParameters betweenRuns withinRuns expGUIhandles
global LevittControl rareEvent errormsg

figure(handles.figure1) % guarantee subject GUI visibility

% ignore call if program is not ready
if ~strcmp(experiment.status,'waitingForGO'), return, end

set(handles.pushbuttonGO,'visible','off')

% append message to expGUI message box to alert experimenter that the user
% is active
addToMsg('Starting new trial',0)

cla(expGUIhandles.axes1),  title(''); % stimulus
cla(expGUIhandles.axes2),  title(''); % WRV track
drawnow

betweenRuns.runNumber=betweenRuns.runNumber + 1;

withinRuns.trialNumber=1;
withinRuns.variableValue=...
    stimulusParameters.WRVstartValues(betweenRuns.runNumber);
% add random jitter to start level
if ~experiment.singleShot
    % SS or single shot allows the user to precisely set the WRV
    withinRuns.variableValue=withinRuns.variableValue +...
        (rand-0.5)*stimulusParameters.jitterStartdB;
end

withinRuns.peaks=[];
withinRuns.troughs=[];
withinRuns.levelList=[];
withinRuns.meanEstTrack=[];
withinRuns.bestSlopeK=[];
withinRuns.bestGain=[];
withinRuns.bestVMin=[];
withinRuns.forceThreshold=NaN;
withinRuns.responseList=[];
withinRuns.caughtOut=0;
withinRuns.wrongButton=0;
withinRuns.catchTrialCount=0;
withinRuns.thresholdEstimateTrack=[];

withinRuns.beginningOfPhase2=0;
withinRuns.nowInPhase2=0;
withinRuns.thisIsRepeatTrial=0;

rareEvent.Euclid=NaN;
rareEvent.bestGain=NaN;
rareEvent.bestVMin=NaN;
rareEvent.thresholddB=0;
rareEvent.bestPaMindB=NaN;
rareEvent.predictionLevels=[];
rareEvent.predictionsRE=[];

LevittControl.sequence=[];

% on-screen count of number of runs still to complete
trialsToGo=length(betweenRuns.var1Sequence)-betweenRuns.runNumber;
set(handles.toGoCounter,'string', trialsToGo);

switch experiment.threshEstMethod
    case {'2I2AFC++', '2I2AFC+++'}
        % For 2I2AFC the buttons need to be on the screen ab initio
        Levitt2      % inititalize Levitt2 procedure
end

switch experiment.ear
    case{'left', 'right','diotic', 'dichoticLeft','dichoticRight'}
        % allow subject time to recover from 'go' press
        pause(experiment.clickToStimulusPause)
end

errormsg=nextStimulus(handles);				% get the show on the road

% terminate if there is any kind of problem
if ~isempty(errormsg)
    % e.g. limits exceeded, clipping
    disp(errormsg)
    runCompleted(handles)
    return
end
% return route is variable (see intro to this function)

% -----------------------------------------------------buttonBox_callback
function buttonBox_callback(obj, info)
% deals with a button press on the button box.

global experiment
global serobj
global subjectGUIHandles

% do not accept callback if one is already in process
if strcmp(experiment.buttonBoxStatus,'busy')
    disp(' ignored button press')
    return % to quiescent state
end
experiment.buttonBoxStatus='busy';

% identify the code of the button pressed
buttonPressedNo = fscanf(serobj,'%c',1);

% This is the map from the button to the Cedrus codes
switch experiment.buttonBoxType
    case 'horizontal'
        pbGo='7'; 		pb0='1';
        pb1='2';		pb2='3';
        pbRepeat='4';	pbWrong='6';	pbBlank='5';
    case 'square'
        pbGo='7';		pb0='1';
        pb1='3';		pb2='4';
        pbRepeat='8';	pbWrong='6';	pbBlank='5';
end

% decide what to do
switch experiment.status
    case {'presentingStimulus', 'waitingForStart', 'trialcompleted', ...
            'endOfExperiment'}
        disp(' ignored button press')

    case 'waitingForGO'
        % i.e. waiting for new run
        if strcmp(buttonPressedNo,pbGo)			% only GO button  accepted
            startNewRun(subjectGUIHandles)
        else
            disp(' ignored button press')
        end

    case 'waitingForResponse'
        % response to stimuli
        switch buttonPressedNo
            case pb0						% button 0 (top left)
                switch experiment.threshEstMethod
                    case {'2I2AFC++', '2I2AFC+++'}
                        disp(' ignored button press')
                    otherwise
                        set(subjectGUIHandles.pushbutton0,...
                            'backgroundcolor','r')
                        pause(.1)
                        set(subjectGUIHandles.pushbutton0,...
                            'backgroundcolor',get(0,...
                            'defaultUicontrolBackgroundColor'))
                        userSelects0or1(subjectGUIHandles)
                end

            case pb1						% button 1 (bottom left)
                switch experiment.threshEstMethod
                    case {'2I2AFC++', '2I2AFC+++'}
                        userSelects0or1(subjectGUIHandles)
                    otherwise
                        set(subjectGUIHandles.pushbutton1,...
                            'backgroundcolor','r')
                        pause(.1)
                        set(subjectGUIHandles.pushbutton1,...
                            'backgroundcolor',get(0,...
                            'defaultUicontrolBackgroundColor'))
                        userSelects0or1(subjectGUIHandles)
                end

            case pb2						% button 2 (bottom right)
                switch experiment.threshEstMethod
                    case {'2I2AFC++', '2I2AFC+++'}
                        userSelects2 (subjectGUIHandles)
                    otherwise
                        set(subjectGUIHandles.pushbutton2,...
                            'backgroundcolor','r')
                        pause(.1)
                        set(subjectGUIHandles.pushbutton2,...
                            'backgroundcolor',get(0,...
                            'defaultUicontrolBackgroundColor'))
                        userSelects2 (subjectGUIHandles)
                end

            case pbRepeat                   % extreme right button
                switch experiment.threshEstMethod
                    case {'2I2AFC++', '2I2AFC+++'}
                        disp(' ignored button press')
                    otherwise

                        set(subjectGUIHandles.pushbuttoNotSure,...
                            'backgroundcolor','r')
                        pause(.1)
                        set(subjectGUIHandles.pushbuttoNotSure,...
                            'backgroundcolor',get(0,...
                            'defaultUicontrolBackgroundColor'))
                        userSelectsPleaseRepeat (subjectGUIHandles)
                end

            case {pbWrong, pbBlank}
                disp(' ignored button press')

            otherwise						% unrecognised button
                disp('ignored button press')
        end									% end (button press number)
    otherwise
        disp('ignored button press')
end											% experiment status

% button box remains 'busy' until after the stimulus has been presented
experiment.buttonBoxStatus='not busy';

% -------------------------------------------------- pushbuttonGO_Callback
function pushbuttonGO_Callback(hObject, eventdata, handles)
% This is a mouse click path
% GO function is also called directly from button box
%  and from MAP model and stats model

set(handles.pushbuttonGO,'visible','off')
startNewRun(handles)

% ---------------------------------------------------pushbutton0_Callback
function pushbutton0_Callback(hObject, eventdata, handles)
global experiment
% This is a mouse click path

% ignore 0 button if 2I2AFC used
if findstr(experiment.threshEstMethod,'2I2AFC')
    return  % to quiescent state
end

% userDoesNotHearTarget(handles)		% only possible interpretation
userDecides(handles, false)

% -------------------------------------------------- pushbutton1_Callback
function pushbutton1_Callback(hObject, eventdata, handles)
userSelects0or1(handles)				% also called from buttonBox

% ---------------------------------------------------pushbutton2_Callback
function pushbutton2_Callback(hObject, eventdata, handles)
userSelects2(handles)					% also called from buttonBox

% --------------------------------------------- pushbuttoNotSure_Callback
function pushbuttoNotSure_Callback(hObject, eventdata, handles)
userSelectsPleaseRepeat(handles)		% also called from buttonBox

% -------------------------------------------------- pushbutton3_Callback
function pushbutton3_Callback(hObject, eventdata, handles)

% ------------------------------------------------- pushbutton19_Callback
function pushbutton19_Callback(hObject, eventdata, handles)
% should be invisible (ignore)

% --------------------------------------- pushbuttonWrongButton_Callback
function pushbuttonWrongButton_Callback(hObject, eventdata, handles)
userSelectsWrongButton(handles)

% --------------------------------------- editdigitInput_Callback
function editdigitInput_Callback(hObject, eventdata, handles)
userSelects0or1(handles)				% after digit string input



% ----------------------------------------------------- userSelects0or1
function userSelects0or1(handles)
global experiment withinRuns

switch experiment.threshEstMethod
    case {'2I2AFC++', '2I2AFC+++'}
        switch withinRuns.stimulusOrder
            case 'targetFirst';
                %                 userHearsTarget(handles)
                userDecides(handles, true)
            otherwise
                %                 userDoesNotHearTarget(handles)
                userDecides(handles, false)
        end
    otherwise
        % single interval
        % 0 or 1 are treated as equivalent (i.e. target is not heard)
        userDecides(handles, false)
end
% return to pushButton1 callback

% ----------------------------------------------------- userSelects2
function userSelects2(handles)
global experiment withinRuns
switch experiment.threshEstMethod
    case {'2I2AFC++', '2I2AFC+++'}
        switch withinRuns.stimulusOrder
            case 'targetSecond';
                %                 userDoesNotHearTarget(handles)
                userDecides(handles, true)
            otherwise
                %                 userHearsTarget(handles)
                userDecides(handles, false)
        end
    otherwise
        % single interval (2 targets heard)
        userDecides(handles, true)
end
% return to pushButton2 callback

% ----------------------------------------------------- ---- userDecides
function userDecides(handles, saidYes)
global experiment stimulusParameters betweenRuns withinRuns
global rareEvent logistic psy levelsBinVector errormsg

if experiment.singleShot
    return  % not clear why this should be here
end

% ignore click if not 'waitingForResponse'
if ~strcmp(experiment.status,'waitingForResponse')
    disp('ignored click')
    return % to userSelects
end

% speech reception threshold
if strcmp(stimulusParameters.targetType,'digitStrings')
    % read triple digits from userGUI
    digitsInput=get(handles.editdigitInput,'string');
    % must be three digits
    if ~(length(digitsInput)==3)
        addToMsg(['error message: Wrong no of digits'], 0, 1)
        set(handles.textMSG,'string', 'Wrong no of digits', ...
            'BackgroundColor','r', 'ForegroundColor', 'w')
        set(handles.editdigitInput,'string','')
        return
    end
    % obtain correct answer from file name
    x=stimulusParameters.digitString;
    idx=find(x=='O'); x(idx)='0'; % replace 'oh' with zero
    disp([x '   ' digitsInput])
    if x==digitsInput
        saidYes=1;  % i.e. correct response
    else
        saidYes=0;  % i.e  wrong response
    end
    set(handles.editdigitInput,'string','')
    set(handles.editdigitInput,'visible','off')
    pause(0.1)
end

% no button presses accepted while processing
experiment.status='processingResponse';

% catch trials. Restart trial if caught
if withinRuns.catchTrial
    if saidYes
        disp('catch trial - caught out')
        withinRuns.caughtOut=withinRuns.caughtOut+1;

        % special: estimate caught out rate by allowing the trial
        %  to continue after catch
        if stimulusParameters.catchTrialBaseRate==0.5
            %  To use this facility, set the catchTrialRate and the
            %   catchTrialBaseRate both to 0.5
            %    update false positive rate
            betweenRuns.caughtOut(betweenRuns.runNumber)=...
                withinRuns.caughtOut;
            plotProgressThisTrial(handles)
            nextStimulus(handles);
            return
        end

        % Punishment: caught out restarts the trial
        set(handles.frame1,'backgroundColor','r')
        set(handles.pushbuttonGO, ...
            'visible','on', 'backgroundcolor','y') % and go again
        msg=[{'Start again: catch trial error'}, {' '},...
            {'Please,click on the GO button'}];
        set(handles.textMSG,'string',msg)
        [y,fs]=wavread('ding.wav');
        if ispc
            wavplay(y/100,fs)
        else
            sound(y/100,fs)
        end
            

        % raise catch trial rate temporarily.
        %  this is normally reduced on each new trial (see GO)
        stimulusParameters.catchTrialRate=...
            stimulusParameters.catchTrialRate+0.1;
        if stimulusParameters.catchTrialRate>0.5
            stimulusParameters.catchTrialRate=0.5;
        end
        fprintf('stimulusParameters.catchTrialRate= %6.3f\n', ...
            stimulusParameters.catchTrialRate)

        betweenRuns.caughtOut(betweenRuns.runNumber)=...
            1+betweenRuns.caughtOut(betweenRuns.runNumber);
        betweenRuns.runNumber=betweenRuns.runNumber-1;
        experiment.status='waitingForGO';
        return % unwind and wait for button press
    else % (said No)
        % user claims not to have heard target.
        % This is good as it was not present.
        % So, repeat the stimulus (possibly with target)
        %  and behave as if the last trial did not occur
        errormsg=nextStimulus(handles);

        % terminate if there is any kind of problem
        if ~isempty(errormsg)
            % e.g. limits exceeded, clipping
            disp(['Error nextStimulus: ' errormsg])
            runCompleted(handles)
            return
        end
        return      % no further action - next trial
    end
end     % of catch trial

% Real target: analyse the response, make tracks and define next stim.

% Define response and update response list
if saidYes
    % target was heard, so response=1;
    withinRuns.responseList=[withinRuns.responseList 1];	% 'heard it!'
else
    % target was not hear heard, so response=0;
    withinRuns.responseList=[withinRuns.responseList 0];
end
withinRuns.levelList=[withinRuns.levelList withinRuns.variableValue];
trialNumber=length(withinRuns.responseList);

% keep track of peaks and troughs;
% identify direction of change during initial period
if saidYes
    % default step size before first reversal
    WRVinitialStep=-stimulusParameters.WRVinitialStep;
    WRVsmallStep=-stimulusParameters.WRVsmallStep;
    % if the previous direction was 'less difficult', this must be a peak
    if strcmp(withinRuns.direction,'less difficult') ...
            && length(withinRuns.levelList)>1
        withinRuns.peaks=[withinRuns.peaks withinRuns.variableValue];
    end
    withinRuns.direction='more difficult';
else
    % said 'no'
    % default step size before first reversal
    WRVinitialStep=stimulusParameters.WRVinitialStep;
    WRVsmallStep=stimulusParameters.WRVsmallStep;

    % if the previous direction was 'up', this must be a peak
    if strcmp(withinRuns.direction,'more difficult') ...
            && length(withinRuns.levelList)>1
        withinRuns.troughs=[withinRuns.troughs withinRuns.variableValue];
    end
    withinRuns.direction='less difficult';
end

% phase 2 is all the levels after and incuding the first reversal
%  plus the level before that
% Look for the end of phase 1
if ~withinRuns.nowInPhase2 && length(withinRuns.peaks)+ ...
        length(withinRuns.troughs)>0
    % define phase 2
    withinRuns.beginningOfPhase2=trialNumber-1;
    withinRuns.nowInPhase2=1;
    WRVsmallStep=WRVinitialStep/2;
end

if withinRuns.nowInPhase2
    % keep a record of all levels and responses in phase 2 only
    withinRuns.levelsPhaseTwo=...
        withinRuns.levelList(withinRuns.beginningOfPhase2:end);
    withinRuns.responsesPhaseTwo=...
        withinRuns.responseList(withinRuns.beginningOfPhase2:end);
else
    withinRuns.levelsPhaseTwo=[];
end

% get (or substitute) threshold estimate
switch experiment.threshEstMethod
    case {'2I2AFC++', '2I2AFC+++'}
        % for plotting psychometric function only
        if withinRuns.beginningOfPhase2>0
            [psy, levelsBinVector, logistic, rareEvent]= ...
                bestFitPsychometicFunctions...
                (withinRuns.levelsPhaseTwo,  withinRuns.responsesPhaseTwo);
        end

        if ~isempty(withinRuns.peaks) && ~isempty(withinRuns.troughs)
            thresholdEstimate= ...
                mean([mean(withinRuns.peaks) mean(withinRuns.troughs)]);
        else
            thresholdEstimate=NaN;
        end

    otherwise
        % single interval methods
        try
            % using the s trial after the first reversal
            [psy, levelsBinVector, logistic, rareEvent]= ...
                bestFitPsychometicFunctions(withinRuns.levelsPhaseTwo,...
                withinRuns.responsesPhaseTwo);
        catch
            logistic.bestThreshold=NaN;
        end
end

if withinRuns.nowInPhase2
    % save tracks of threshold estimates for plotting andprinting
    switch experiment.functionEstMethod
        case {'logisticLS', 'logisticML'}
            if withinRuns.nowInPhase2
                withinRuns.meanEstTrack=...
                    [withinRuns.meanEstTrack ...
                    mean(withinRuns.levelsPhaseTwo)];
                withinRuns.thresholdEstimateTrack=...
                    [withinRuns.thresholdEstimateTrack ...
                    logistic.bestThreshold];
            end
        case 'rareEvent'
            withinRuns.meanEstTrack=...
                [withinRuns.meanEstTrack rareEvent.thresholddB];
            withinRuns.thresholdEstimateTrack=...
                [withinRuns.thresholdEstimateTrack logistic.bestThreshold];
        case 'peaksAndTroughs'
            withinRuns.meanEstTrack=...
                [withinRuns.meanEstTrack thresholdEstimate];
            withinRuns.thresholdEstimateTrack=...
                [withinRuns.thresholdEstimateTrack thresholdEstimate];
    end
end

% special discomfort condition
% run is completed when subject hits '2' button
switch experiment.paradigm
    case 'discomfort'
        if saidYes
            runCompleted(handles)
            return
        end
end

% choose the next level for the stimulus
switch experiment.threshEstMethod
    case {'2I2AFC++', '2I2AFC+++'}
        if saidYes
            [WRVinitialStep, msg]=Levitt2('hit', withinRuns.variableValue);
        else
            [WRVinitialStep, msg]=Levitt2('miss',withinRuns.variableValue);
        end

        % empty message means continue as normal
        if ~isempty(msg)
            runCompleted(handles)
            return
        end
        newWRVvalue=withinRuns.variableValue-WRVinitialStep;

    case {'MaxLikelihood', 'oneIntervalUpDown'}
        % run completed by virtue of number of trials
        % or restart because listener is in trouble
        if length(withinRuns.levelsPhaseTwo)== experiment.maxTrials
            % Use bonomial test to decide if there is an imbalance in the
            % number of 'yes'es and 'no's
            yesCount=sum(withinRuns.responseList);
            noCount=length(withinRuns.responseList)-yesCount;
            z=abs(yesCount-noCount)/(yesCount+noCount)^0.5;
            if z>1.96
                betweenRuns.resets=betweenRuns.resets+1;
                disp([ 'reset / z= ' num2str( z)  ...
                    '   Nresets= ' num2str( betweenRuns.resets) ] )
                withinRuns.peaks=[];
                withinRuns.troughs=[];
                withinRuns.levelList=withinRuns.levelList(end);
                withinRuns.meanEstTrack=withinRuns.meanEstTrack(end);
                withinRuns.forceThreshold=NaN;
                withinRuns.responseList=withinRuns.responseList(end);
                withinRuns.beginningOfPhase2=0;
                withinRuns.nowInPhase2=0;
                withinRuns.thresholdEstimateTrack=...
                    withinRuns.thresholdEstimateTrack(end);
            else
                runCompleted(handles)
                return
            end
        end

        % set new value for WRV
        if withinRuns.nowInPhase2
            % phase 2
            currentMeanEst=withinRuns.thresholdEstimateTrack(end);
            switch experiment.threshEstMethod
                case 'MaxLikelihood'
                    newWRVvalue=currentMeanEst;
                case {'oneIntervalUpDown'}
                    newWRVvalue=withinRuns.variableValue+WRVsmallStep;
            end
        else
            % phase 1
            if withinRuns.variableValue+2*WRVinitialStep>...
                    stimulusParameters.WRVlimits(2)
                % use smaller steps when close to maximum
                WRVinitialStep=WRVinitialStep/2;
            end
            newWRVvalue=withinRuns.variableValue+WRVinitialStep;
        end
    otherwise
        error(  'assessment method not recognised')
end

switch experiment.paradigm
    % prevent unrealistic gap durations 'gapDetection' tasks.
    % Note that the gap begins when the ramp ends not when stimulus ends
    case 'gapDetection'
        if newWRVvalue<-2*stimulusParameters.rampDuration
            newWRVvalue=-2*stimulusParameters.rampDuration;
            addToMsg('gap duration fixed at - 2 * ramp!',1, 1)
        end
end

withinRuns.variableValue=newWRVvalue;
withinRuns.trialNumber=withinRuns.trialNumber+1;

% Trial continues
plotProgressThisTrial(handles)

% next stimulus and so the cycle continues
errormsg=nextStimulus(handles);
% after the stimulus is presented, control returns here and the system
% waits for user action.

% terminate if there is any kind of problem
if ~isempty(errormsg)
    % e.g. limits exceeded, clipping
    disp(['Error nextStimulus: ' errormsg])
    runCompleted(handles)
    return
end

% ------------------------------------------------ userSelectsPleaseRepeat
function userSelectsPleaseRepeat(handles)
global experiment withinRuns
% ignore click if not 'waitingForResponse'
if ~strcmp(experiment.status,'waitingForResponse')
    disp('ignored click')
    return
end
% Take no action other than to make a
%  tally of repeat requests
experiment.pleaseRepeat=experiment.pleaseRepeat+1;
withinRuns.thisIsRepeatTrial=1;
nextStimulus(handles);

% ------------------------------------------------ userSelectsWrongButton
function userSelectsWrongButton(handles)
global withinRuns experiment
% restart is the simplest solution for a 'wrong button' request
withinRuns.wrongButton=withinRuns.wrongButton+1;
set(handles.pushbuttonGO, 'visible','on', 'backgroundcolor','y')
msg=[{'Start again: wrong button pressed'}, {' '},...
    {'Please,click on the GO button'}];
set(handles.textMSG,'string',msg)
experiment.status='waitingForGO';

% ------------------------------------------------- plotProgressThisTrial
function plotProgressThisTrial(handles)
% updates GUI: used for all responses

global experiment stimulusParameters betweenRuns withinRuns expGUIhandles
global  psy levelsBinVector binFrequencies rareEvent logistic statsModel

% plot the levelTrack and the threshold track

% Panel 2
% plot the levelList
axes(expGUIhandles.axes2); cla
plot( withinRuns.levelList,'o','markerFaceColor','k'), hold on
% plot the best threshold estimate tracks
if length(withinRuns.meanEstTrack)>=1
    % The length of the levelList is 2 greater than number of thresholds
    ptr=withinRuns.beginningOfPhase2+1;
    plot(ptr: ptr+length(withinRuns.meanEstTrack)-1, ...
        withinRuns.meanEstTrack, 'r')
    plot( ptr: ptr+length(withinRuns.thresholdEstimateTrack)-1, ...
        withinRuns.thresholdEstimateTrack, 'g')
    hold off
    estThresh=withinRuns.thresholdEstimateTrack(end);
    switch experiment.threshEstMethod
        % add appropriate labels to subject GUI buttons
        case {'2I2AFC++', '2I2AFC+++'}
            title([stimulusParameters.WRVname ' = ' ...
                num2str(withinRuns.variableValue, '%5.1f')])
        otherwise
            title([stimulusParameters.WRVname ' = ' ...
                num2str(withinRuns.variableValue, '%5.1f') ...
                ';    TH= ' num2str(estThresh, '%5.1f')])
    end
end
xlim([0 experiment.maxTrials+withinRuns.beginningOfPhase2]);
ylim(stimulusParameters.WRVlimits)
grid on

% Panel 4: Summary of threshold estimates (not used here)
% Estimates from previous runs are set in 'runCompleted'
% It is only necessary to change title showing runs/trials remaining
axes(expGUIhandles.axes4)
runsToGo=length(betweenRuns.var1Sequence)-betweenRuns.runNumber;
if withinRuns.beginningOfPhase2>0
    trialsToGo= experiment.singleIntervalMaxTrials(1) ...
        + withinRuns.beginningOfPhase2- withinRuns.trialNumber;
    title(['trials remaining = ' num2str(trialsToGo) ...
        ':    runs to go= ' num2str(runsToGo)])
end

% plot psychometric function   - panel 5
axes(expGUIhandles.axes5), cla
plot(withinRuns.levelList, withinRuns.responseList,'b.'), hold on
ylim([0 1])
title('')

switch experiment.threshEstMethod
    case {'MaxLikelihood', 'oneIntervalUpDown'}
        if withinRuns.beginningOfPhase2>0
            % display only when in phase 2.
            withinRuns.levelsPhaseTwo=...
                withinRuns.levelList(withinRuns.beginningOfPhase2:end);
            withinRuns.responsesPhaseTwo=...
                withinRuns.responseList(withinRuns.beginningOfPhase2:end);

            % organise data as psychometric function
            [psy, levelsBinVector, binFrequencies]= ...
                psychometricFunction(withinRuns.levelsPhaseTwo,...
                withinRuns.responsesPhaseTwo, experiment.psyBinWidth);

            % Plot the function
            %   point by point with circles of appropiate weighted size
            hold on,
            for i=1:length(psy)
                plot(levelsBinVector(i), psy(i), 'ro', ...
                    'markersize', 50*binFrequencies(i)/sum(binFrequencies))
            end
            % save info for later
            betweenRuns.psychometicFunction{betweenRuns.runNumber}=...
                [levelsBinVector; psy];

            % fitPsychometric functions is  computed in 'userDecides'
            % plot(rareEvent.predictionLevels, rareEvent.predictionsRE,'k')
            plot(logistic.predictionLevels, logistic.predictionsLOG, 'r')
            plot(rareEvent.predictionLevels, rareEvent.predictionsRE, 'k')
            if ~isnan(logistic.bestThreshold )
                xlim([ 0 100 ])
                title(['k= ' num2str(logistic.bestK, '%6.2f') ' g= '...
                    num2str(rareEvent.bestGain,'%6.3f') '  A=' ...
                    num2str(rareEvent.bestVMin,'%8.1f')])
            else
                title(' ')
            end

            switch experiment.ear
                %plot green line for statsModel a priori model
                case 'statsModelLogistic'
                    % plot proTem logistic (green) used by stats model
                    p= 1./(1+exp(-statsModel.logisticSlope...
                        *(levelsBinVector-logistic.bestThreshold)));
                    if experiment.psyFunSlope<0, p=1-p;end
                    titleText=[ ',  statsModel: logistic'];
                    hold on,    plot(levelsBinVector, p,'g')
                case  'statsModelRareEvent'
                    pressure=28*10.^(levelsBinVector/20);
                    p=1-exp(-stimulusParameters.targetDuration...
                        *(statsModel.rareEvenGain...
                        * pressure-statsModel.rareEventVmin));
                    p(p<0)=0;
                    if experiment.psyFunSlope<0, p=1-p;end
                    hold on,    plot(levelsBinVector, p,'g')
            end %(estMethod)
        end

    otherwise           % 2A2IFC
        message3= ...
            ([ 'peaks='  num2str(withinRuns.peaks) ...
            'troughs='  num2str(withinRuns.troughs)]);
        ylimRM([-0.1 1.1])	% 0=no / 1=yes
        set(gca,'ytick',[0 1], 'yTickLabel', {'no';'yes'})
        ylabel('psychometric function'), xlabel('target level')
        if length(levelsBinVector)>1
            xlim([ min(levelsBinVector) max(levelsBinVector)])
            xlim([ 0 100])
        end
end

% command window summary
% Accumulate things to say in the message window
message1= (['responses:      ' num2str(withinRuns.responseList,'%9.0f')]);
switch experiment.paradigm
    % more decimal places needed on GUI
    case { 'gapDetection', 'frequencyDiscrimination', 'forwardMaskingD'}
        message2= ([stimulusParameters.WRVname  ...
            ':       ' num2str(withinRuns.levelList,'%7.3f')]);
        message3= (['Thresh (logistic mean):   ' ...
            num2str(withinRuns.thresholdEstimateTrack,'%7.3f')]);
    otherwise
        message2= ([stimulusParameters.WRVname ':      ' ...
            num2str(withinRuns.levelList,'%7.1f')]);
        message3= (['Thresh (logistic mean):   ' ...
            num2str(withinRuns.thresholdEstimateTrack,'%7.1f')]);
end

addToMsg(str2mat(message1, message2, message3), 0)

% -----------------------------------------------------runCompleted
function runCompleted(handles)
% Used at the end of each run
global experiment stimulusParameters betweenRuns withinRuns
global rareEvent expGUIhandles
% disp('run completed')

experiment.status='runCompleted';
% quick update after final trial just to make sure
plotProgressThisTrial(handles)

switch experiment.ear
    case {'statsModelLogistic', 'statsModelRareEvent','MAPmodel', ...
            'MAPmodelMultiCh','MAPmodelSingleCh', 'MAPmodelListen'}
        % no changes required if model used
    otherwise
        set(handles.frame1,'visible','off')
        set(handles.pushbuttoNotSure,'visible','off')
        set(handles.pushbuttonWrongButton,'visible','off')
        set(handles.pushbutton3,'visible','off')
        set(handles.pushbutton2,'visible','off')
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton0,'visible','off')
        set(handles.pushbuttonGO,'visible','off')
end

if isnan(withinRuns.forceThreshold)
    % the experiment has been aborted for some reason
    threshold=withinRuns.forceThreshold;
    stdev=NaN;
    logistic.bestK=NaN;
    logistic.bestThreshold=NaN;
    medianThreshold=NaN;
    meanThreshold=NaN;
else
    % use only phase 2 levels and responses for calculating thresholds
    withinRuns.levelsPhaseTwo=...
        withinRuns.levelList(withinRuns.beginningOfPhase2:end);
    withinRuns.responsesPhaseTwo=...
        withinRuns.responseList(withinRuns.beginningOfPhase2:end);
    [psy, levelsPhaseTwoBinVector, logistic, rareEvent]= ...
        bestFitPsychometicFunctions...
        (withinRuns.levelsPhaseTwo, withinRuns.responsesPhaseTwo);

    % plot final psychometric function
    axes(expGUIhandles.axes5),cla
    hold on, plot(rareEvent.predictionLevels, rareEvent.predictionsRE, 'k')
    hold on, plot(logistic.predictionLevels, logistic.predictionsLOG, 'r')
    % organise data as psychometric function
    [psy, levelsBinVector, binFrequencies]= ...
        psychometricFunction(withinRuns.levelsPhaseTwo,...
        withinRuns.responsesPhaseTwo, experiment.psyBinWidth);
    %   point by point with circles of appropiate weighted size
    hold on,
    for i=1:length(psy)
        plot(levelsBinVector(i), psy(i), 'ro', ...
            'markersize', 50*binFrequencies(i)/sum(binFrequencies))
    end

    % experimental
    medianThreshold=median(withinRuns.levelsPhaseTwo);
    warning off
    meanThreshold=mean(withinRuns.levelsPhaseTwo);

    % identify the current threshold estimate
    switch experiment.paradigm
        case 'discomfort'
            % most recent value (not truely a mean value)
            threshold=withinRuns.levelList(end);
            stdev=NaN;
        otherwise
            switch experiment.threshEstMethod
                case {'MaxLikelihood', 'oneIntervalUpDown'}
                    % last value in the list
                    %                     threshold=withinRuns.meanEstTrack(end);
                    threshold=withinRuns.thresholdEstimateTrack(end);
                    stdev=NaN;

                case {'2I2AFC++', '2I2AFC+++'}
                    % use peaks and troughs
                    try		% there may not be enough values to use
                        peaksUsed=experiment.peaksUsed;
                        threshold=...
                            mean(...
                            [withinRuns.peaks(end-peaksUsed+1:end) ...
                            withinRuns.troughs(end-peaksUsed+1:end)]);
                        stdev=...
                            std([withinRuns.peaks(end-peaksUsed +1:end) ...
                            withinRuns.troughs(end-peaksUsed:end)]);
                    catch
                        threshold=NaN;
                        stdev=NaN;
                    end
            end
    end
end

% Store thresholds
betweenRuns.thresholds=[betweenRuns.thresholds threshold];
betweenRuns.thresholds_mean=[betweenRuns.thresholds_mean meanThreshold];
betweenRuns.thresholds_median=...
    [betweenRuns.thresholds_median medianThreshold];
betweenRuns.forceThresholds=...
    [betweenRuns.forceThresholds withinRuns.forceThreshold];

% count observations after the startup phase for record keeping
betweenRuns.observationCount=...
    [betweenRuns.observationCount length(withinRuns.levelList)];
betweenRuns.timesOfFirstReversals=...
    [betweenRuns.timesOfFirstReversals withinRuns.beginningOfPhase2];
betweenRuns.catchTrials=...
    [betweenRuns.catchTrials withinRuns.catchTrialCount];

% add variable length tracks to cell arrays
if withinRuns.beginningOfPhase2>0
    betweenRuns.bestThresholdTracks{length(betweenRuns.thresholds)}=...
        withinRuns.thresholdEstimateTrack;
    betweenRuns.levelTracks{length(betweenRuns.thresholds)}=...
        withinRuns.levelList(withinRuns.beginningOfPhase2:end);
    betweenRuns.responseTracks{length(betweenRuns.thresholds)}=...
        withinRuns.responseList(withinRuns.beginningOfPhase2:end);
else
    betweenRuns.bestThresholdTracks{length(betweenRuns.thresholds)}=[];
    betweenRuns.levelTracks{length(betweenRuns.thresholds)}=[];
    betweenRuns.responseTracks{length(betweenRuns.thresholds)}=[];
end

betweenRuns.bestGain=[betweenRuns.bestGain rareEvent.bestGain];
betweenRuns.bestVMin=[betweenRuns.bestVMin rareEvent.bestVMin];
betweenRuns.bestPaMin=[betweenRuns.bestPaMin rareEvent.bestPaMindB];
betweenRuns.bestLogisticM=...
    [betweenRuns.bestLogisticM logistic.bestThreshold];
betweenRuns.bestLogisticK=[betweenRuns.bestLogisticK logistic.bestK];

resultsSoFar=[betweenRuns.var1Sequence(betweenRuns.runNumber)'...
    betweenRuns.var2Sequence(betweenRuns.runNumber)'...
    betweenRuns.thresholds(betweenRuns.runNumber)'
    ];

fprintf('%10.3f \t%10.3f \t%10.1f  \n', resultsSoFar')

switch experiment.ear
    case {'left', 'right', 'diotic', 'dichoticLeft','dichoticRight'}
        disp(['caught out= ' num2str(betweenRuns.caughtOut)])
end

% plot history of thresholds in panel 4
axes(expGUIhandles.axes4), cla
plotColors='rgbmckywrgbmckyw';
for i=1:length(betweenRuns.thresholds)
    faceColor=plotColors(floor(i/length(betweenRuns.variableList1)-.01)+1);
    switch betweenRuns.variableName1
        case {'targetFrequency', 'maskerRelativeFrequency'}
            if min(betweenRuns.var1Sequence)>0
                semilogx(betweenRuns.var1Sequence(i), ...
                    betweenRuns.thresholds(i),  'o', ...
                    'markerSize', 5,'markerFaceColor',faceColor)
            else
                plot(betweenRuns.var1Sequence(1:betweenRuns.runNumber),  ...
                    betweenRuns.thresholds,  'o', ...
                    'markerSize', 5,'markerFaceColor',faceColor)
                plot(betweenRuns.var1Sequence(i),  ...
                    betweenRuns.thresholds(i),  'o', ...
                    'markerSize', 5,'markerFaceColor',faceColor)
            end
        otherwise
            plot(betweenRuns.var1Sequence(i),  ...
                betweenRuns.thresholds(i),  'o', 'markerSize', 5,...
                'markerFaceColor',faceColor)
    end
    hold on
end
xlimRM([ min(betweenRuns.variableList1) max(betweenRuns.variableList1) ])
ylim(stimulusParameters.WRVlimits)
ylabel('thresholds')
xlabel(betweenRuns.variableName1)
set(gca,'ytick', [0 20 40 60 80 100])
try
    % problems if only one x value
    set(gca,'XTick', sort(betweenRuns.variableList1))
catch
end
grid on, set(gca,'XMinorGrid', 'off')

% final run?
if betweenRuns.runNumber==length(betweenRuns.var1Sequence)
    % yes, end of experiment
    fileName=['savedData/' experiment.name experiment.date ...
        experiment.paradigm];
    % 	save (fileName, 'experiment', 'stimulusParameters', 'betweenRuns', 'withinRuns', 'variableNames', 'paradigmNames', 'LevittControl')
    disp('Experiment completed')

    % update subject GUI to acknowledge end of run
    subjGUImsg=[{'Experiment completed'}, {' '}, {'Thank you!'}];
    set(handles.textMSG,'string', subjGUImsg   )
    % play 'Tada'
    [y,fs,nbits]=wavread('TADA.wav');
    musicGain=10^(stimulusParameters.musicLeveldB/20);
    y=y*musicGain;
    if ispc, wavplay(y/100,fs, 'async'), else sound(y/100,fs, 'async'), end

    % update experimenter GUI
    addToMsg('Experiment completed.',1)

    printReport
    experiment.status='endOfExperiment';
    return
else
    % No, hang on.
    switch experiment.ear
        case {'statsModelLogistic', 'statsModelRareEvent','MAPmodel', ...
                'MAPmodelMultiCh','MAPmodelSingleCh', 'MAPmodelListen'}
            % no changes required if model used
        otherwise
            % decrement catchTrialRate towards baseRate
            stimulusParameters.catchTrialRate=...
                stimulusParameters.catchTrialBaseRate + ...
                (stimulusParameters.catchTrialRate...
                -stimulusParameters.catchTrialBaseRate)...
                *(1-exp(-stimulusParameters.catchTrialTimeConstant));
            fprintf('stimulusParameters.catchTrialRate= %6.3f\n', ...
                stimulusParameters.catchTrialRate)

            % and go again
            set(handles.pushbuttonGO,'backgroundcolor','y')
            set(handles.frame1,'visible','off')
            set(handles.pushbuttonGO,'visible','on')
            msg=[{'Ready to start new trial'}, {' '},...
                {'Please,click on the GO button'}];
            set(handles.textMSG,'string',msg)
    end
    experiment.status='waitingForGO';
    %     fprintf('\n')

    [y,fs,nbits]=wavread('CHIMES.wav');
    musicGain=10^(stimulusParameters.musicLeveldB/20);
    y=y*musicGain;
    if ispc, wavplay(y/100,fs, 'async'), else sound(y/100,fs, 'async'), end

end

% -----------------------------------------------------MAPmodelRunsGUI
% The computer presses the buttons
function 	MAPmodelRunsGUI(handles)
global experiment stimulusParameters method expGUIhandles
global AN_IHCsynapseParams
method=[];

while strcmp(experiment.status,'waitingForGO')
    % no catch trials for MAP model
    experiment.allowCatchTrials=0;

    % initiates run and plays first stimulus and it returns
    %  without waiting for button press
    startNewRun(handles)

    % show sample Rate on GUI; it must be set in MAPparams ##??
    set(expGUIhandles.textsampleRate,'string',...
        num2str(stimulusParameters.sampleRate))

    if experiment.singleShot % ##??
        AN_IHCsynapseParams.plotSynapseContents=1;
    else
        AN_IHCsynapseParams.plotSynapseContents=0;
    end

    % continuous loop until the program stops itself
    while strcmp(experiment.status,'waitingForResponse')
        %  NB at this point the stimulus has been played
        pause(0.1)  % to allow interrupt with CTRL/C

        switch experiment.ear
            case { 'MAPmodelListen'}
                % flash the buttons to show model response
                set(handles.pushbutton1,'backgroundcolor','y','visible','on')
                set(handles.pushbutton2,'backgroundcolor','y','visible','on')
        end

        % Analayse the current stimulus using MAP
        [modelResponse earObject]= MAPmodel;

        if experiment.stop || experiment.singleShot
            % trap for single trial or user interrupt using 'stop' button.
            experiment.status= 'waitingForStart';
%             experiment.stop=0;
            errormsg='manually stopped';
            addToMsg(errormsg,1)
            return
        end

        switch modelResponse
            case 1
                %     userDoesNotHearTarget(handles)
                switch experiment.ear
                    case {'MAPmodelListen'}
                        % illuminate appropriate button
                        set(handles.pushbutton1,...
                            'backgroundcolor','r','visible','on')
                        set(handles.pushbutton2,'backgroundcolor','y')
                end
                userDecides(handles, false)
                if experiment.singleShot, return, end

            case 2
                %   userHearsTarget(handles)
                switch experiment.ear
                    case {'MAPmodelListen'}
                        % illuminate appropriate button (DEMO only)
                        set(handles.pushbutton2,'backgroundcolor',...
                            'r','visible','on')
                        set(handles.pushbutton1,'backgroundcolor','y')
                end

                switch experiment.paradigm
                    case 'discomfort'
                        % always treat discomfort as 'not heard'
                        userDecides(handles, false)
                    otherwise
                        userDecides(handles, true)
                end
            otherwise
                % probably an abort
                return
        end
    end
end

% -------------------------------------------------------MAPmodel
function [modelResponse, MacGregorResponse]=MAPmodel

global experiment stimulusParameters audio withinRuns
% global outerMiddleEarParams DRNLParams AN_IHCsynapseParams
global ICoutput ANdt dt savedBFlist ANprobRateOutput expGUIhandles
global paramChanges

savePath=path;
addpath(['..' filesep 'MAP'], ['..' filesep 'utilities'])
modelResponse=[];
MacGregorResponse=[];

% mono only (column vector)
audio=audio(:,1)';

% if stop button pressed earlier
if experiment.stop
    errormsg='manually stopped';
    addToMsg(errormsg,1)
    return
end

% ---------------------------------------------- run Model
MAPparamsName=experiment.name;
showPlotsAndDetails=experiment.MAPplot;

% important buried constant ##??
AN_spikesOrProbability='spikes';
% AN_spikesOrProbability='probability';

% [response, method]=MAPsequenceSeg(audio, method, 1:8);

if sum(strcmp(experiment.ear,{'MAPmodelMultiCh', 'MAPmodelListen'}))
    % use BFlist specified in MAPparams file
    BFlist= -1;
else
    BFlist=stimulusParameters.targetFrequency;
end
paramChanges=get(expGUIhandles.editparamChanges,'string');
% convert from string to a cell array
eval(paramChanges);

MAP1_14(audio, stimulusParameters.sampleRate, BFlist,...
    MAPparamsName, AN_spikesOrProbability, paramChanges);

if showPlotsAndDetails
    options.printModelParameters=0;
    options.showModelOutput=1;
    options.printFiringRates=1;
    options.showACF=0;
    options.showEfferent=1;
    options.surfProbability=0;
    showMapOptions.surfSpikes=0;
    UTIL_showMAP(options)
end

% No response,  probably caused by hitting 'stop' button
if strcmp(AN_spikesOrProbability,'spikes') && isempty(ICoutput)
    return
end
% ---------------------------------------------------------- end model run

if strcmp(AN_spikesOrProbability,'spikes')
    MacGregorResponse= sum(ICoutput,1);                 % use IC
    dt=ANdt;
    time=dt:dt:dt*length(MacGregorResponse);
else
    % for one channel, ANprobResponse=ANprobRateOutput
    % for multi channel take strongest in any epoch
    nChannels=length(savedBFlist);
    % use only HSR fibers
    ANprobRateOutput=ANprobRateOutput(end-nChannels+1:end,:);
    time=dt:dt:dt*length(ANprobRateOutput);
end

% group delay on unit response - these values are iffy
windowOnsetDelay= 0.004;
windowOffsetDelay= 0.022; % long ringing time

% now find the response of the MacGregor model during the target presentation + group delay
switch experiment.threshEstMethod
    case {'2I2AFC++', '2I2AFC+++'}
        idx= time>stimulusParameters.testTargetBegins+windowOnsetDelay ...
            & time<stimulusParameters.testTargetEnds+windowOffsetDelay;
        nSpikesTrueWindow=sum(MacGregorResponse(:,idx));
        idx=find(time>stimulusParameters.testNonTargetBegins+windowOnsetDelay ...
            & time<stimulusParameters.testNonTargetEnds+windowOffsetDelay);
        if strcmp(AN_spikesOrProbability,'spikes')
            nSpikesFalseWindow=sum(MacGregorResponse(:,idx));
        else
            nSpikesDuringTarget=mean(ANprobRateOutput(end,idx));
            % compare with spontaneous rate
            if nSpikesDuringTarget>ANprobRateOutput(end,1)+10
                nSpikesDuringTarget=1;  % i.e. at leastone MacG spike
            else
                nSpikesDuringTarget=0;
            end
        end
        % nSpikesDuringTarget is +ve when more spikes are found
        %   in the target window
        difference= nSpikesTrueWindow-nSpikesFalseWindow;

        if difference>0
            % hit
            nSpikesDuringTarget=experiment.MacGThreshold+1;
        elseif    difference<0
            % miss (wrong choice)
            nSpikesDuringTarget=experiment.MacGThreshold-1;
        else
            if rand>0.5
                % hit (random choice)
                nSpikesDuringTarget=experiment.MacGThreshold+1;
            else
                % miss (random choice)
                nSpikesDuringTarget=experiment.MacGThreshold-1;
            end
        end
        disp(['level target dummy decision: ' ...
            num2str([withinRuns.variableValue nSpikesTrueWindow ...
            nSpikesFalseWindow  nSpikesDuringTarget], '%4.0f') ] )
    otherwise
        % single interval
        idx=find(time>stimulusParameters.testTargetBegins +windowOnsetDelay...
            & time<stimulusParameters.testTargetEnds+windowOffsetDelay);
        if strcmp(AN_spikesOrProbability,'spikes')
            nSpikesDuringTarget=sum(MacGregorResponse(:,idx));
            timeX=time(idx);
        else
            % probability, use channel with the highest average rate
            nSpikesDuringTarget=max(mean(ANprobRateOutput(:,idx),2));
            if nSpikesDuringTarget>ANprobRateOutput(end,1)+10
                nSpikesDuringTarget=1;
            else
                nSpikesDuringTarget=0;
            end

        end
end

if experiment.MAPplot
    % add vertical lines to indicate target region
    figure(99), subplot(6,1,6)
    hold on
    yL=get(gca,'YLim');
    plot([stimulusParameters.testTargetBegins + windowOnsetDelay ...
        stimulusParameters.testTargetBegins   + windowOnsetDelay],yL,'r')
    plot([stimulusParameters.testTargetEnds   + windowOffsetDelay ...
        stimulusParameters.testTargetEnds     + windowOffsetDelay],yL,'r')
end

% specify unambiguous response
if nSpikesDuringTarget>experiment.MacGThreshold
    modelResponse=2;    % stimulus detected
else
    modelResponse=1;    % nothing heard (default)
end

path(savePath)

% -----------------------------------------------------statsModelRunsGUI
% The computer presses the buttons
function 	statsModelRunsGUI(handles)
% Decision are made at random using a prescribe statistical function
% to set probabilities as a function of signal level.
global experiment

experiment.allowCatchTrials=0;

while strcmp(experiment.status,'waitingForGO')
    % i.e. waiting for new run
    if experiment.stop
        % user has requested an abort
        experiment.status= 'waitingForStart';
        addToMsg('manually stopped',1)
        return
    end

    % initiates run and plays first stimulus and it returns
    %  without waiting for button press
    % NB stimulus is not actually generated (for speed)
    startNewRun(handles)

    while strcmp(experiment.status,'waitingForResponse')
        % create artificial response here
        modelResponse=statsModelGetResponse;
        switch modelResponse
            case 1
                %                 userDoesNotHearTarget(handles)
                userDecides(handles, false)
            case 2
                %                 userHearsTarget(handles)
                userDecides(handles, true)
        end
    end
end

% -----------------------------------------------------statsModelGetResponse
function modelResponse=statsModelGetResponse(handles)
global experiment  withinRuns  statsModel stimulusParameters
% use the generating function to decide if a detection occurs or not

% pause(0.1)  % to allow stopping with CTRL/C but slows things down

% first compute the probability that a detection occurs
switch experiment.ear
    case {'statsModelLogistic'}
        prob= 1./(1+exp(-statsModel.logisticSlope.*(withinRuns.variableValue-statsModel.logisticMean)));
        %         if experiment.psyFunSlope<0,
        %             prob=1-prob;
        %         end

    case 'statsModelRareEvent'
        if experiment.psyFunSlope<0
            addToMsg('statsModelRareEvent cannot be used with negative slope',0)
            error('statsModelRareEvent cannot be used with negative slope')
        end

        % standard formula is prob = 1 – exp(-d (g P – A))
        % here A->A;  To find Pmin use A/gain
        pressure=28*10^(withinRuns.variableValue/20);
        gain=statsModel.rareEvenGain;
        A=statsModel.rareEventVmin;
        d=stimulusParameters.targetDuration;
        gP_Vmin=gain*pressure-A;
        if gP_Vmin>0
            prob=1-exp(-d*(gP_Vmin));
        else
            prob=0;
        end
end

% Use the probability to choose whether or not a detection has occurred
switch experiment.threshEstMethod
    case {'MaxLikelihood', 'oneIntervalUpDown'}
        if rand<prob
            modelResponse=2; %bingo
        else
            modelResponse=1; %nothing heard
        end

    case {'2I2AFC++', '2I2AFC+++'}
        if rand<prob
            modelResponse=2; %bingo
        else %if the stimulus is not audible, take a 50:50 chance of getting it right
            if rand<0.5
                modelResponse=2; %bingo
            else
                modelResponse=1; %nothing heard
            end
        end
end


% ------------------------------------------------------- printTabTable
function printTabTable(M, headers)
% printTabTable prints a matrix as a table with tabs
%headers are optional
%headers=strvcat('firstname', 'secondname')
%  printTabTable([1 2; 3 4],strvcat('a1','a2'));

if nargin>1
    [r c]=size(headers);
    for no=1:r
        fprintf('%s\t',headers(no,:))
    end
    fprintf('\n')
end

[r c]=size(M);

for row=1:r
    for col=1:c
        if row==1 && col==1 && M(1,1)==-1000
            %   Print nothing (tab follows below)
        else
            fprintf('%s',num2str(M(row,col)))
        end
        if col<c
            fprintf('\t')
        end
    end
    fprintf('\n')
end

% ------------------------------------------------------- xlimRM
function xlimRM(x)
try
    xlim([x(1) x(2)])
catch
end

% ------------------------------------------------------- ylimRM
function ylimRM(x)
try
    ylim([x(1) x(2)])
catch
end


function editdigitInput_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -----------------------------------------------------buttonBoxIntitialize
function buttonBoxIntitialize
% initialize button box
global  serobj
try
    fclose(serobj);
catch
end

try
    serobj = serial('COM4') ;           	% Creating serial port object now its connected to COM4 !!! button boxes in booths are connected to COM2
    serobj.Baudrate = 9600;           		% Set the baud rate at the specific value
    set(serobj, 'Parity', 'none') ;     	% Set parity as none
    set(serobj, 'Databits', 8) ;          	% set the number of data bits
    set(serobj, 'StopBits', 1) ;         	% set number of stop bits as 1
    set(serobj, 'Terminator', 'CR') ; 		% set the terminator value to carriage return
    set(serobj, 'InputBufferSize', 512) ;  	% Buffer for read operation, default it is 512
    set(serobj,'timeout',10);           	% 10 sec timeout on button press
    set(serobj, 'ReadAsyncMode', 'continuous')
    set(serobj, 'BytesAvailableFcn', @buttonBox_callback)
    set(serobj, 'BytesAvailableFcnCount', 1)
    set(serobj, 'BytesAvailableFcnMode', 'byte')
    % set(serobj, 'BreakInterruptFcn', '@buttonBox_Calback')

    fopen(serobj);
    buttonBoxStatus=get(serobj,'status');
catch
    disp('** no button box found - use mouse **')
end
