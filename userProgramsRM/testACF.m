% function [LP_SACF dt lags SACF]= testACF
% testACF is a *script* to demonstrate the smoothed ACF of
% Balaguer-Ballestera, E. Denham, S.L. and Meddis, R. (2008).
%
% Convert this to a *function* by uncommenting the first line
%  The function returns the LP_SACF matrix plotted in Figure 96.
%  If a function is used, the following outputs are returned:
%   LP_SACF:  smoothed SACF (lags x time matrix)
%   dt: time interval between successive columns of LP_SACF
%   lags: lags used in computing LP_SACF
%   SACF: unsmoothed SACFs
%
% A range of options are supplied in the early part of the program
%
% #1
% Identify the model parameter file (in 'MAPparamsName')
%
% #2
% Identify the kind of model required (in 'AN_spikesOrProbability')
%  'probability' is recommended for ACF work
%
% #3
% Choose between a harmonic complex or file input
%  by commenting out unwanted code
%
% #4
% Set the signal rms level (in leveldBSPL)
%
% #5
% Identify the model channel BFs in the vector 'BFlist'.
%
% #6
% Last minute changes to the model parameters can be made using
%  the cell array of strings 'paramChanges'.
%  This is used here to control the details of the ACF computations
%  Read the notes in this section for more information
%
% displays:
% Figure 97 shows the AN response to the stimulus. this is a channel x time
%  display. The z-axis (and colour) is the AN fiber firing rate
%
% Figure 96 shows the LP_SACF-matrix, the smoothed SACF.
%
% Figure 89 shows a summary of the evolution of the unsmoothed SACF
%  over time. If you wish to take a snapshot of the LP_SACF-matrix at a
%  particular time, this figure can help identify when to take it.
%  The index on the y-axis, identifies the required row numbers
%    of the LP_SACF or SACF matrix, e.g. LP_SACF(:,2000)
%
%  On request, (filteredSACFParams.plotACFs=1) Figure 89 shows the channel
%   by channel ACFs at intervals during the computation as a movie.
%  The number of ACF displays is controlled by 'plotACFsInterval'
%   and the movie can be slowed or speeded up using 'plotMoviePauses'
%   (see paramChanges section below).

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% This global will find results from MAP1_14
global savedInputSignal ANprobRateOutput ANoutput dt dtSpikes savedBFlist
% This global,from model parameter file
global filteredSACFParams

% User sets up requirements
%%  #1 parameter file name
MAPparamsName='Normal';                 % recommended


%% #2 probability (fast) or spikes (slow) representation: select one
% AN_spikesOrProbability='spikes';
%   or
AN_spikesOrProbability='probability';   % recommended

%% #3 A. harmonic sequence or B. speech file input
% Comment out unwanted code
% A. harmonic tone (Hz) - useful to demonstrate a broadband sound
sampleRate= 44100;              % recommended 44100
signalType= 'tones';
duration=0.100;                 % seconds
beginSilence=0.020;
endSilence=0.020;
rampDuration=.005;              % raised cosine ramp (seconds)

% toneFrequency is a vector of component frequencies
F0=120;
toneFrequency= [3*F0 4*F0 5*F0];

%   or
% B. file input
% signalType= 'file';
% fileName='Oh No';
% fileName='twister_44kHz';

%% #4 rms level
leveldBSPL= 100;                  % dB SPL (80 for Lieberman)

%% #5 number of channels in the model
%   21-channel model (log spacing of BFs)
numChannels=21;
lowestBF=250; 	highestBF= 5000;
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

%% #6 change model parameters
% Parameter changes can be used to change one or more model parameters
%  *after* the MAPparams file has been read (see manual)

% Take control of ACF parameters
%  The filteredACF parameters are set in the MAPparamsNormal file
%  However, it is convenient to change them here leving the file intacta
minPitch=	400; maxPitch=	3000; numPitches=200;
maxLag=1/minPitch; minLag=1/maxPitch;
lags= linspace(minLag, maxLag, numPitches);

paramChanges={...
    'filteredSACFParams.lags=lags;     % autocorrelation lags vector;',...
    'filteredSACFParams.acfTau=	2;     % (Wiegrebe) time constant ACF;',...
    'filteredSACFParams.lambda=	0.12;  % slower filter to smooth ACF;',...
    'filteredSACFParams.plotACFs=1;    % plot ACFs while computing;',...
    'filteredSACFParams.plotACFsInterval=0.01;',...
    'filteredSACFParams.plotMoviePauses=.1;  ',...
    'filteredSACFParams.usePressnitzer=0; % attenuates ACF at  long lags;',...
    'filteredSACFParams.lagsProcedure=  ''useAllLags'';',...
    };

% Notes:
% acfTau:   time constant of unsmoothed ACF
% lambda:   time constant of smoothed ACFS
% plotACFs: plot ACFs during computation (0 to switch off, for speed)
% plotACFsInterval: sampling interval for plots
% plotMoviePauses:  pause duration between frames to allow viewing
% usePressnitzer:   gives low weights to long lags
% lagsProcedure:    used to fiddle with output (ignore)

%% delare 'showMap' options to control graphical output
% see UTIL_showMAP for more options
showMapOptions=[];
% showMapOptions.showModelOutput=0;     % plot of all stages
showMapOptions.surfAN=1;                % surface plot of HSR response
showMapOptions.PSTHbinwidth=0.001;      % smoothing for PSTH

if exist('fileName','var')
    % needed for labeling plot
    showMapOptions.fileName=fileName;
end

%% Generate stimuli
switch signalType
    case 'tones'
        % Create tone stimulus
        dt=1/sampleRate; % seconds
        time=dt: dt: duration;
        inputSignal=sum(sin(2*pi*toneFrequency'*time), 1);
        amp=10^(leveldBSPL/20)*28e-6;   % converts to Pascals (peak)
        inputSignal=amp*inputSignal;
        % apply ramps
        % catch rampTime error
        if rampDuration>0.5*duration, rampDuration=duration/2; end
        rampTime=dt:dt:rampDuration;
        ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
            ones(1,length(time)-length(rampTime))];
        inputSignal=inputSignal.*ramp;
        ramp=fliplr(ramp);
        inputSignal=inputSignal.*ramp;
        % add silence
        intialSilence= zeros(1,round(beginSilence/dt));
        finalSilence= zeros(1,round(endSilence/dt));
        inputSignal= [intialSilence inputSignal finalSilence];

    case 'file'
        %% file input simple or mixed
        [inputSignal sampleRate]=wavread(fileName);
        dt=1/sampleRate;
        inputSignal=inputSignal(:,1);
        targetRMS=20e-6*10^(leveldBSPL/20);
        rms=(mean(inputSignal.^2))^0.5;
        amp=targetRMS/rms;
        inputSignal=inputSignal*amp;
end

wavplay(inputSignal, sampleRate)

%% run the model
dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])

fprintf('\n')
disp(['Signal duration= ' num2str(length(inputSignal)/sampleRate)])
disp([num2str(numChannels) ' channel model: ' AN_spikesOrProbability])
disp('Computing MAP ...')

MAP1_14(inputSignal, sampleRate, BFlist, ...
    MAPparamsName, AN_spikesOrProbability, paramChanges);


%% The model run is now complete. Now display the results
% display the AN response
UTIL_showMAP(showMapOptions)

% compute ACF
switch AN_spikesOrProbability
    case 'probability'
        % use only HSR fibers
        inputToACF=ANprobRateOutput(end-length(savedBFlist)+1:end,:);
    otherwise
        inputToACF=ANoutput;
        dt=dtSpikes;
end

disp ('computing ACF...')

% read paramChanges to get new filteredSACFParams
for i=1:length(paramChanges)
    eval(paramChanges{i});
end

[LP_SACF BFlist SACF]= filteredSACF(inputToACF, dt, savedBFlist, ...
    filteredSACFParams);
disp(' ACF done.')

%% plot original waveform on summary/smoothed ACF plot
figure(96), clf
subplot(3,1,3)
t=dt*(1:length(savedInputSignal));
plot(t,savedInputSignal, 'k')
xlim([0 t(end)])
title(['stimulus: ' num2str(leveldBSPL, '%4.0f') ' dB SPL']);

% plot SACF
figure(96)
subplot(2,1,1)
imagesc(LP_SACF)
colormap bone
ylabel('periodicities (Hz)'), xlabel('time (s)')
title(['smoothed SACF. (periodicity x time)'])
% y-axis specifies pitches (1/lags)
% Force MATLAB to show the lowest pitch
postedYvalues=[1 get(gca,'ytick')]; set(gca,'ytick',postedYvalues)
pitches=1./filteredSACFParams.lags;
set(gca,'ytickLabel', round(pitches(postedYvalues)))
% x-axis is time at which LP_SACF is samples
[nCH nTimes]=size(LP_SACF);
t=dt:dt:dt*nTimes;
tt=get(gca,'xtick');
set(gca,'xtickLabel', round(100*t(tt))/100)

%% On a new figure show a cascade of SACFs
figure(89), clf
% select 100 samples;
[r c]=size(SACF);
step=round(c/100);
idx=step:step:c;

UTIL_cascadePlot(SACF(:,idx)', 1./pitches)

xlabel('lag (s)'), ylabel('time pointer -->')
title(' SACF summary over time')
yValues=get(gca,'yTick');
set(gca,'yTickLabel', num2str(yValues'*100))

path(restorePath)

