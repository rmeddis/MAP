function demoTwisterProbability

% MAPdemo runs the MATLAB auditory periphery model (MAP1_14) as far as
%  the AN (probabilities) or IC (spikes) with graphical output

% Things you might want to change; #1 - #5

%%  #1 parameter file name
MAPparamsName='Normal';


%% #2 probability (fast) or spikes (slow) representation
% AN_spikesOrProbability='spikes';
% or
AN_spikesOrProbability='probability';


%% #3 pure tone, harmonic sequence or speech file input
signalType= 'tones';
duration=0.100;                 % seconds
duration=0.020;                 % seconds
sampleRate= 64000;
% toneFrequency= 250:250:8000;    % harmonic sequence (Hz)
toneFrequency= 2000;            % or a pure tone (Hz8

rampDuration=.005;              % seconds

% or
signalType= 'file';
fileName='twister_44kHz';
% fileName='new-da-44khz';


%% #4 rms level
% signal details
leveldBSPL=70;                  % dB SPL


%% #5 number of channels in the model
%   21-channel model (log spacing)
numChannels=21;
lowestBF=250; 	highestBF= 8000; 
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

%   or specify your own channel BFs
% BFlist=toneFrequency;


%% #6 change model parameters
paramChanges=[];

% or
% Parameter changes can be used to change one or more model parameters
%  *after* the MAPparams file has been read
% This example declares only one fiber type with a calcium clearance time
% constant of 80e-6 s (HSR fiber) when the probability option is selected.
% switch AN_spikesOrProbability
%     case 'probability'
%         paramChanges={'IHCpreSynapseParams.tauCa=80e-6;'};
%     otherwise
%         paramChanges=[];
% end

%% delare showMap options
showMapOptions=[];  % use defaults

% or (example: show everything including an smoothed SACF output
    showMapOptions.showModelParameters=1;
    showMapOptions.showModelOutput=1;
    showMapOptions.printFiringRates=1;
    showMapOptions.showACF=0;
    showMapOptions.showEfferent=1;

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
showMAP(showMapOptions)

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

