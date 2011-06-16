function demoTwisterProbability

% MAPdemo runs the MATLAB auditory periphery model (MAP1_14) as far as
%  the AN (probabilities) or IC (spikes) with graphical output

restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])

%%  #1 parameter file name
MAPparamsName='Normal';


%% #2 probability (fast) 
AN_spikesOrProbability='probability';


%% #3  speech file input
signalType= 'file';
fileName='twister_44kHz';


%% #4 rms level
leveldBSPL=60;        % dB SPL


%% #5 number of channels in the model
%   21-channel model (log spacing)
numChannels=21;
lowestBF=250; 	highestBF= 8000;
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

%% #6 no change to model parameters
paramChanges=[];


%% delare showMap options
showMapOptions.printModelParameters=1;
showMapOptions.showModelOutput=1;
showMapOptions.printFiringRates=1;
showMapOptions.showACF=0;
showMapOptions.showEfferent=0;
showMapOptions.surfProbability=1;       % 2D plot of HSR response 

%% Generate stimuli

dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'])
switch signalType
    case 'tones'
        inputSignal=createMultiTone(sampleRate, toneFrequency, ...
            leveldBSPL, duration, rampDuration);

    case 'file'
        [inputSignal sampleRate]=wavread(fileName);
        inputSignal(:,1);
        targetRMS=20e-6*10^(leveldBSPL/20);
        rms=(mean(inputSignal.^2))^0.5;
        amp=targetRMS/rms;
        inputSignal=inputSignal*amp;
end


%% run the model
tic

MAP1_14(inputSignal, sampleRate, BFlist, ...
    MAPparamsName, AN_spikesOrProbability, paramChanges);

toc

% the model run is now complete. Now display the results
UTIL_showMAP(showMapOptions)

toc
path(restorePath)


function inputSignal=createMultiTone(sampleRate, toneFrequency, ...
    leveldBSPL, duration, rampDuration)
% Create pure tone stimulus
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

% add 10 ms silence
silence= zeros(1,round(0.03/dt));
% inputSignal= [silence inputSignal silence];

