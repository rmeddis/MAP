function varargout = speechmanipulation(varargin)
% SPEECHMANIPULATION MATLAB code for speechmanipulation.fig
%      SPEECHMANIPULATION, by itself, creates a new SPEECHMANIPULATION or raises the existing
%      singleton*.
%
%      H = SPEECHMANIPULATION returns the handle to a new SPEECHMANIPULATION or the handle to
%      the existing singleton*.
%
%      SPEECHMANIPULATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPEECHMANIPULATION.M with the given input arguments.
%
%      SPEECHMANIPULATION('Property','Value',...) creates a new SPEECHMANIPULATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before speechmanipulation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to speechmanipulation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help speechmanipulation

% Last Modified by GUIDE v2.5 05-Oct-2011 11:21:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @speechmanipulation_OpeningFcn, ...
                   'gui_OutputFcn',  @speechmanipulation_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before speechmanipulation is made visible.
function speechmanipulation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to speechmanipulation (see VARARGIN)
global wavfilescell;
global rateaxeshandle; %handle to axis6 is needed, if radiobuttons12 and 13 change


% Choose default command line output for speechmanipulation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes speechmanipulation wait for user response (see UIRESUME)
% uiwait(handles.figure1);
wavfiles = dir(['..' filesep 'wavFileStore' filesep '*.wav']);
for iCounter = 1:length(wavfiles)
    wavfilescell{iCounter} = wavfiles(iCounter).name;
end
set(handles.popupmenu1,'String',wavfilescell);
set(handles.popupmenu2,'String',wavfilescell);
set(handles.uipanel8,'SelectionChangeFcn',@radiobuttonselected);
rateaxeshandle = handles.axes6;

% --- Outputs from this function are returned to the command line.
function varargout = speechmanipulation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function edit1_Callback(hObject, eventdata, handles)
function edit2_Callback(hObject, eventdata, handles)
function edit3_Callback(hObject, eventdata, handles)
function edit4_Callback(hObject, eventdata, handles)
function edit5_Callback(hObject, eventdata, handles)
function edit6_Callback(hObject, eventdata, handles)
function edit7_Callback(hObject, eventdata, handles)
function edit8_Callback(hObject, eventdata, handles)
function edit9_Callback(hObject, eventdata, handles)
function edit10_Callback(hObject, eventdata, handles)
function edit11_Callback(hObject, eventdata, handles)
function edit12_Callback(hObject, eventdata, handles)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
%% CHOOSE AND LOAD WAVEFORM
global wavfilescell
global originalsignal
global actualsignal

%get the number from the popupmenu
selectedfileno = get(handles.popupmenu1,'Value');

%load it
[waveform,sfreq] = wavread(['..' filesep 'wavFileStore' filesep wavfilescell{selectedfileno}]);

%if it is from OLLO then cut it
if regexp(wavfilescell{selectedfileno}, 'S[0-9]*[A-Z]_L[0-9]*_V[0-9]_M[0-9]_N[0-9]_CS0.wav') %pattern for OLLO-files
    waveform = cutsignal(waveform,sfreq,'a_a');
end

%store it in global variables
actualsignal.waveform = waveform;
actualsignal.sfreq = sfreq;
originalsignal.waveform = waveform;
originalsignal.sfreq = sfreq;

%plot it
multipleplot(actualsignal,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% LOAD WAVEFORM
global wavfilescell
global originalsignal
global actualsignal

%get the number from the popupmenu
selectedfileno = get(handles.popupmenu1,'Value');

%load it
[waveform,sfreq] = wavread(['..' filesep 'wavFileStore' filesep wavfilescell{selectedfileno}]);

%if it is from OLLO then cut it
if regexp(wavfilescell{selectedfileno}, 'S[0-9]*[A-Z]_L[0-9]*_V[0-9]_M[0-9]_N[0-9]_CS0.wav') %pattern for OLLO-files
    waveform = cutsignal(waveform,sfreq,'a_a');
end

%store it in global variables
actualsignal.waveform = waveform;
actualsignal.sfreq = sfreq;
originalsignal.waveform = waveform;
originalsignal.sfreq = sfreq;

%plot it
multipleplot(actualsignal,handles)


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% PLAY THE ACTUAL (MANIPULATED) SIGNAL
global actualsignal
global toplay

toplay = audioplayer(actualsignal.waveform,actualsignal.sfreq);
play(toplay);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% MANIPULATE SIGNAL
global wavfilescell
global originalsignal
global actualsignal

actualsignal = originalsignal; %reset

%% speed factor change
speedfactor = str2num(get(handles.edit1,'String'));
newfakesfreq=actualsignal.sfreq*speedfactor;
if speedfactor == 1
else
    actualsignal.waveform = resample(actualsignal.waveform,actualsignal.sfreq,newfakesfreq);
end

%% low and high pass filter (butterworth)
% get the desired filter from radio buttons
radiofilterhandles=get(handles.uipanel2,'Children');
for iCounter = 1:length(radiofilterhandles)
    if get(radiofilterhandles(iCounter),'Value') == 1
        filtertype = get(radiofilterhandles(iCounter),'String');
    end
end
% read parameters
tmp = get(handles.edit3,'String');
cutofffreq=str2num(tmp{1});
tmp = get(handles.edit4,'String');
filterorder=str2num(tmp{1});
%do the filter calculations
switch filtertype
    case 'none'
    case 'High pass'
        if isempty(cutofffreq) || isempty(filterorder)
            error('Please specify cutoff frequency and/or filter order')
        else
            [a,b] = butter(filterorder,cutofffreq*2/actualsignal.sfreq,'high');
            actualsignal.waveform=filter(a,b,actualsignal.waveform);
        end
    case 'Low pass'
        if isempty(cutofffreq)|| isempty(filterorder)
            error('Please specify cutoff frequency')
        else
            [a,b] = butter(filterorder,cutofffreq*2/actualsignal.sfreq,'low');
            actualsignal.waveform=filter(a,b,actualsignal.waveform);
        end
end

%additive noise
% get the desired filter from radio buttons
radionoisehandles=get(handles.uipanel4,'Children');
for iCounter = 1:length(radionoisehandles)
    if get(radionoisehandles(iCounter),'Value') == 1
        noisetype = get(radionoisehandles(iCounter),'String');
    end
end      
% read parameters
tmp = get(handles.edit5,'String');
snr=str2num(tmp{1});
%do the noise calculations
switch noisetype
    case 'none'
    case 'White noise'
        if isempty(snr) 
            error('Please specify SNR')
        else
            whitenoise = rand(length(actualsignal.waveform),1)-0.5;
            %set level
            levelsignal = 20*log10(sqrt(mean(actualsignal.waveform.^2)));
            levelnoise = levelsignal-snr;
            whitenoise = whitenoise./sqrt(mean(whitenoise.^2)).*10^(levelnoise/20);
            %20*log10(sqrt(mean(whitenoise.^2)))
            actualsignal.waveform=actualsignal.waveform+whitenoise;
        end
    case 'From file'
        if isempty(snr) 
            error('Please specify SNR')
        else
            selectedfileno = get(handles.popupmenu2,'Value');

            [noisewaveform,noisesfreq] = wavread(['..' filesep 'wavFileStore' filesep wavfilescell{selectedfileno}]);
            if noisesfreq == actualsignal.sfreq
            else
                noisewaveform = resample(noisewaveform,actualsignal.sfreq,noisesfreq);
                noisesfreq = actualsignal.sfreq;
            end
            if length(noisewaveform) < length(actualsignal.waveform)
                warning('Noise waveform too short. Noise is looped without fading at the endings!');
                noisewaveform = repmat(noisewaveform, ceil(length(actualsignal.waveform)/length(noisewaveform)),1);
            end
            noisewaveform=noisewaveform(1:length(actualsignal.waveform));
            %set level
            levelsignal = 20*log10(sqrt(mean(actualsignal.waveform.^2)));
            levelnoise = levelsignal-snr;
            noisewaveform = noisewaveform./sqrt(mean(noisewaveform.^2)).*10^(levelnoise/20);
            %20*log10(sqrt(mean(noisewaveform.^2)))
            %finally add the two signals
            actualsignal.waveform = actualsignal.waveform + noisewaveform;
        end
end

%plot the resulting manipulated waveform
multipleplot(actualsignal,handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% PLAY THE ORIGINAL WAVEFORM
global originalsignal
global toplay
toplay = audioplayer(originalsignal.waveform,originalsignal.sfreq);
play(toplay)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global toplay
stop(toplay);


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% CALCULATE MAP 1.14 auditory nerve probability output
addpath(['..' filesep 'MAP']);
addpath(['..' filesep 'utilities']);
global actualsignal
global AN_HSRoutput

tmp = get(handles.edit6,'String');
level4MAP=str2num(tmp);
if isempty(level4MAP)
    error('Please specify level');
end

tmp = get(handles.edit7,'String');
parameterfile=tmp;
if isempty(parameterfile)
    error('Please specify parameter file');
end

% identify model parameter changes if any
paramChanges=get(handles.edit8,'string');
if ~strcmp(paramChanges, ';'), paramChanges=[paramChanges ';']; end
eval(paramChanges);

%set level
actualsignal.waveform = actualsignal.waveform./sqrt(mean(actualsignal.waveform.^2)).*10^(-(94-level4MAP)/20);
%20*log10(sqrt(mean(actualsignal.waveform.^2))/20e-6) %reference pressure: 20uPa

MAP1_14(actualsignal.waveform,actualsignal.sfreq,-1,parameterfile, ...
    'probability',paramChanges);
global ANprobRateOutput savedBFlist

%take only the HSR fibers
AN_HSRoutput = ANprobRateOutput(size(ANprobRateOutput)/2+1:end,:);

%plot the auditory nerve firing probability as line plot
start_time = size(AN_HSRoutput,2)/2/actualsignal.sfreq-0.025; %start time to plot in s
end_time = start_time+0.05; %plot 50ms
plotIFRAN(AN_HSRoutput,start_time,end_time,actualsignal.sfreq,savedBFlist,handles.axes4);
set(handles.edit11,'String',num2str(1000*size(AN_HSRoutput,2)/2/actualsignal.sfreq));

%plot the fourierhistogram as image plot
formantpattern = fourierautocorrelationhistogram_direct(AN_HSRoutput,actualsignal.sfreq,handles.axes5);
caxis([0 2000]);
colorbar;

%plot the rate output
plotrateOutput(AN_HSRoutput,actualsignal.sfreq,handles.axes6,[33.7 300]);
set(handles.radiobutton12,'Value',1);

%calculate the IPIH
ipih=track_formants_from_IPI_guy(AN_HSRoutput, actualsignal.sfreq);
%the following code assumes that the bin width is 1/actualsignal.sfreq and 
%the temporal step size is 3ms

ipihfreqaxis=1./[1/actualsignal.sfreq:1/actualsignal.sfreq:size(ipih,1)/actualsignal.sfreq];
ipihtimeaxis=[0:3:size(ipih,2)*3];
set(gcf,'CurrentAxes',handles.axes7);
YTickIdx = 1:floor(numel(ipihfreqaxis)/6):numel(ipihfreqaxis);
XTickIdx = 1:floor(numel(ipihtimeaxis)/6):numel(ipihtimeaxis);
imagesc(ipih);
set(gca, 'YTick', YTickIdx);
set(gca, 'YTickLabel', num2str(  ipihfreqaxis(YTickIdx)', '%0.0f' ));
ylabel('best frequency (Hz)')
set(gca, 'XTick', XTickIdx);
set(gca, 'XTickLabel', XTickIdx.*3);
xlabel('Time (ms)');
caxis([0 8e4]); %set color 
colorbar;

%calculate and plot summarized autocorrelation, code from Ray
method.dt=1/actualsignal.sfreq;
method.segmentNo=1;
method.nonlinCF=savedBFlist;

minPitch=	80; maxPitch=	4000; numPitches=100;    % specify lags
pitches=10.^ linspace(log10(minPitch), log10(maxPitch),numPitches);
pitches=fliplr(pitches);
filteredSACFParams.lags=1./pitches;     % autocorrelation lags vector
filteredSACFParams.acfTau=	.003;       % time constant of running ACF
filteredSACFParams.lambda=	0.12;       % slower filter to smooth ACF
filteredSACFParams.lambda=	0.01;       % slower filter to smooth ACF

filteredSACFParams.plotACFs=0;          % special plot (see code)
filteredSACFParams.plotFilteredSACF=0;  % 0 plots unfiltered ACFs
filteredSACFParams.plotMoviePauses=.5;%.3          % special plot (see code)

filteredSACFParams.usePressnitzer=0; % attenuates ACF at  long lags
filteredSACFParams.lagsProcedure=  'useAllLags';
filteredSACFParams.criterionForOmittingLags=3;
filteredSACFParams.plotACFsInterval=50;%200;

% compute ACF
%switch saveAN_spikesOrProbability
%    case 'probability'
        inputToACF=ANprobRateOutput.^0.5;
%    otherwise
%        inputToACF=ANoutput;
%end

disp ('computing ACF...')
t=method.dt*(1:length(actualsignal.waveform));
[P, BFlist, sacf, boundaryValue] = ...
    filteredSACF(inputToACF, method, filteredSACFParams);
P = real(P); %dont know why sometimes P can be (very slightly) complex
disp(' ACF done.')

% SACF
set(gcf,'CurrentAxes',handles.axes8);
imagesc(P)
ylabel('periodicities (Hz)')
xlabel('time (s)')
%title(['running smoothed (root) SACF. ' saveAN_spikesOrProbability ' input'])
pt=[1 get(gca,'ytick')]; % force top xtick to show
set(gca,'ytick',pt)
set(gca,'ytickLabel', round(pitches(pt)))
tt=get(gca,'xtick');
tt=tt(tt<length(t));
set(gca,'xtickLabel', round(100*t(tt))/100)
colorbar;

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% SHOW THE  (PREVIOUSLY STORED) AUDITORY PROFILE CORRESPONDING TO THE
%% PARAMETER FILE
addpath(['..' filesep 'multithreshold 1.46' filesep]);
% read the parameter file
tmp = get(handles.edit7,'String');
parameterfile=tmp;
if isempty(parameterfile)
    error('Please specify parameter file');
end

%what kind of profiles are in the profiles directory?
files_in_profiles = dir(['..' filesep 'profiles' filesep '*.m']);%if exist('',file)
done = 0;

%try to map the specified parameter file to these found in the directory
%and
for iCounter = 1:length(files_in_profiles)
    if strfind(files_in_profiles(iCounter).name,parameterfile)
        %ask the user if this is the correct file - because you'll never
        %know ;)
        stringtoask = ['Is ' files_in_profiles(iCounter).name ' the correct profile file?'];
        Answer = questdlg(stringtoask,'Possible match found','Yes','No','Yes');
        switch Answer,
            case 'Yes'
                plotprofile(files_in_profiles(iCounter).name(1:end-2),'profile_CMA_L');
                done = 1;
            case 'No'
        end %switch
    end %if
end
if ~done
    msgbox('No match for this parameter file found.');
end

function multipleplot(actualsignal,handles)
%% FUNCTION TO DO 3 PLOTS
time_axis = [0:1/actualsignal.sfreq:(length(actualsignal.waveform)-1)/actualsignal.sfreq];

%waveform plot (temporal)
plot(handles.axes1,time_axis,actualsignal.waveform);
set(get(handles.axes1,'XLabel'),'String','Time (s)');
set(get(handles.axes1,'YLabel'),'String','Amplitude');
set(handles.axes1,'XLim',[0 length(actualsignal.waveform)/actualsignal.sfreq]);
highestabsampl = max(abs(actualsignal.waveform))*1.05;
set(handles.axes1,'YLim',[-highestabsampl highestabsampl]);

%average spectrum plot
complspectrum = fft(actualsignal.waveform)/(length(actualsignal.waveform));
frequency = [0:1/time_axis(end):1/(2*(time_axis(2)-time_axis(1)))];
plot(handles.axes2,frequency,20*log10(sqrt(2)*abs(complspectrum(1:round(length(complspectrum)/2)))));
set(get(handles.axes2,'XLabel'),'String','frequency (Hz)');
set(get(handles.axes2,'YLabel'),'String','fourier amplitude (dB)');
set(handles.axes2,'XLim',[100 10000]);
set(handles.axes2,'XScale','log');

%spectrogram plot (10 ms temporal resolution)
[s,f,t] = spectrogram(actualsignal.waveform,hann(round(0.01*actualsignal.sfreq)),[],[],actualsignal.sfreq); %10ms short term windows
set(gcf,'CurrentAxes',handles.axes3);
if max(f)>10000
    [tmp,idx] = min(abs(f-10000));
    imagesc(t,f(1:idx),20*log10(abs(s(1:idx,:))));
else
    imagesc(t,f,20*log10(abs(s)));
end
axis xy;
set(get(handles.axes3,'YLabel'),'String','frequency (Hz)');
set(get(handles.axes3,'XLabel'),'String','Time (s)');


function plotrateOutput(ratepattern,sfreq,axeshandle,colorrange)
%% calculate rate representation and plot, code from Nick
global savedBFlist
%calculate rate representation
ANsmooth = [];%Cannot pre-allocate a size as it is unknown until the enframing
hopSize = 10; %ms
winSize = 25; %ms
winSizeSamples = round(winSize*sfreq/1000);
hann = hanning(winSizeSamples);
hopSizeSamples = round(hopSize*sfreq/1000);
for chan = 1:size(ratepattern,1)
    f = enframe(ratepattern(chan,:), hann, hopSizeSamples);
    ANsmooth(chan,:) = mean(f,2)';
end

%plot rate representation
time_axis_rate=[0:hopSize/1000:size(ratepattern,2)/sfreq];
set(gcf,'CurrentAxes',axeshandle);
YTickIdx = 1:floor(numel(savedBFlist)/6):numel(savedBFlist);
XTickIdx = 1:floor(numel(time_axis_rate)/6):numel(time_axis_rate);
imagesc(ANsmooth);
axis xy;
set(gca, 'YTick', YTickIdx);
set(gca, 'YTickLabel', num2str(  savedBFlist(YTickIdx)', '%0.0f' ));
ylabel('best frequency (Hz)')
set(gca, 'XTick', XTickIdx);
set(gca, 'XTickLabel', XTickIdx.*10);
xlabel('Time (ms)');
caxis(colorrange); %set color from average spontaneous rate to a maximum of 600
colorbar;


% --- Executes if the radiobuttons hsr and lsr are changed
function radiobuttonselected(source, eventdata)

global rateaxeshandle
global ANprobRateOutput
global actualsignal

selected = get(get(source,'SelectedObject'),'String');

if isempty(ANprobRateOutput)
    msgbox('Please calculate AN pattern first!');
elseif strcmp(selected,'HSR')
    AN_HSRoutput = ANprobRateOutput(size(ANprobRateOutput)/2+1:end,:);
    colorrange = [33.7 300];
    plotrateOutput(AN_HSRoutput,actualsignal.sfreq,rateaxeshandle,colorrange);
else
    AN_LSRoutput = ANprobRateOutput(1:size(ANprobRateOutput)/2,:);
    colorrange = [10.7 100];
    plotrateOutput(AN_LSRoutput,actualsignal.sfreq,rateaxeshandle,colorrange);
end


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%plot the auditory nerve firing probability as line plot after grabbing the
%time value from edit11
global AN_HSRoutput
global actualsignal
global savedBFlist

middletime = str2num(get(handles.edit11,'String'))/1000;

start_time = min([size(AN_HSRoutput,2)/actualsignal.sfreq-0.01 max([0 middletime-0.025])]); %start time to plot in s
end_time = max([start_time+0.01 min([start_time+0.05 size(AN_HSRoutput,2)/actualsignal.sfreq-0.01])]); %plot max 50ms
plotIFRAN(AN_HSRoutput,start_time,end_time,actualsignal.sfreq,savedBFlist,handles.axes4);
