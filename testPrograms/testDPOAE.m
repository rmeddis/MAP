function [frequencies fft_ampdB]= ...
    testDPOAE (leveldBSPL, toneFrequencies)
% MAPdemo runs the MATLAB auditory periphery model
%
% The OAE is simulated by combining the output from all DRNL channels
%
% arguments leveldBSPL and toneFrequencies are optional
%  defaults are 70 and [5000 6000]
%
% e.g. 
%  MAPdemoMultiChOAE (60, [3000 4000])


global dt DRNLoutput
dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'], ['..' filesep 'utilities'], ...
    ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'],...
    ['..' filesep 'testPrograms'])

% set parameter file here
paramsName='Normal';
% choose probability because spikes not used to evaluate BM
AN_spikesOrProbability='probability';
% add parameter changes here. paramchanges is a cell array of command
% strings
paramChanges={};

% DRNL channels
lowestBF=1000; 	highestBF= 8000; 	numChannels=41;
% includes BFs at 250 500 1000 2000 4000 8000 (for 11, 21, 31 BFs)
%  the output from all these filters will be combined to form the OAE
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

if nargin<2
    toneFrequencies= 2000;            % single pure tone test
    toneFrequencies=[ 2000 3000]; 	%  F1 F2 for DPOAEs
    toneFrequencies=[ 5000 6000]; 	%  F1 F2
end
duration=0.05;                  % seconds
duration=0.05;                  % seconds
rampDuration=.005;

if nargin<1
    leveldBSPL=70;                  % dB SPL
end
amp=10^(leveldBSPL/20)*28e-6;   % converts to Pascals (peak level)

% Create pure stimulus
sampleRate= 100000;
dt=1/sampleRate;
time=dt: dt: duration;
inputSignal=sum(sin(2*pi*toneFrequencies'*time), 1);
inputSignal=amp*inputSignal;

% apply ramps
if rampDuration>0.5*duration, rampDuration=duration/2; end
rampTime=dt:dt:rampDuration;
ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
    ones(1,length(time)-length(rampTime))];
inputSignal=inputSignal.*ramp;  % at the beginning
ramp=fliplr(ramp);
inputSignal=inputSignal.*ramp;  % and at the end

% add 10 ms silence
silenceDuration=0.01;
silence= zeros(1,round(silenceDuration/dt));
inputSignal= [silence inputSignal silence];
time=dt: dt: dt*length(inputSignal);

%% delare 'showMap' options to control graphical output
showMapOptions.printModelParameters=0;   % prints all parameters
showMapOptions.showModelOutput=1;       % plot of all stages
showMapOptions.printFiringRates=1;      % prints stage activity levels
showMapOptions.showACF=0;               % shows SACF (probability only)
showMapOptions.showEfferent=0;          % tracks of AR and MOC
showMapOptions.surfAN=0;       % 2D plot of HSR response
showMapOptions.surfSpikes=0;            % 2D plot of spikes histogram
showMapOptions.ICrates=0;               % IC rates by CNtauGk


MAP1_14(inputSignal, 1/dt, BFlist, ...
    paramsName, AN_spikesOrProbability, paramChanges);

UTIL_showMAP(showMapOptions)
pause(0.1)

% use this to produce a comnplete record of model parameters
% UTIL_showAllMAPStructures

OAE=sum(DRNLoutput);
figure(5),subplot(2,1,1)
plot(time,OAE)
title(['F=' num2str(toneFrequencies)])
[fft_powerdB, fft_phase, frequencies, fft_ampdB]= UTIL_FFT(OAE, dt, 1e-15);
idx=find(frequencies<1e4);

figure(5),subplot(2,1,2)
plot(frequencies(idx),fft_ampdB(idx))
title ('FFT of OAE')
ylabel('dB')
ylim([0 100])
grid on

path(restorePath);
