function test_MAP1_14Hopkins
% test_MAP1_14 is a general purpose test routine that can be adjusted to
% test a number of different applications of MAP1_14
%
% A range of options are supplied in the early part of the program
%
% One use of the function is to create demonstrations; filenames <demoxx>
%  to illustrate particular features
%
% #1
% Identify the file (in 'MAPparamsName') containing the model parameters
%
% #2
% Identify the kind of model required (in 'AN_spikesOrProbability').
%  A full brainstem model (spikes) can be computed or a shorter model
%  (probability) that computes only so far as the auditory nerve
%
% #3
% Choose between a tone signal or file input (in 'signalType')
%
% #4
% Set the signal rms level (in leveldBSPL)
%
% #5
% Identify the channels in terms of their best frequencies in the vector
%  BFlist.
%
% Last minute changes to the parameters fetched earlier can be made using
%  the cell array of strings 'paramChanges'.
%  Each string must have the same format as the corresponding line in the
%  file identified in 'MAPparamsName'
%
% When the demonstration is satisfactory, freeze it by renaming it <demoxx>

dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])

%%  #1 parameter file name
MAPparamsName='Normal';


%% #2 probability (fast) or spikes (slow) representation
AN_spikesOrProbability='spikes';
%   or
AN_spikesOrProbability='probability';


%% #3 pure tone, harmonic sequence or speech file input
signalType= 'tones';
sampleRate= 50000;
duration=0.100;                 % seconds
rampDuration=.005;              % raised cosine ramp (seconds)
beginSilence=0.050;               
endSilence=0.050;                  
toneFrequency= 1000;            % or a pure tone (Hz)

%   or
% harmonic sequence (Hz)
F0=1000/11;
toneFrequency= 10*F0:F0:12*F0;    

%   or
% signalType= 'file';
% fileName='twister_44kHz';



%% #4 rms level
% signal details
leveldBSPL= [44 50 44];                  % dB SPL (80 for Lieberman)
leveldBSPL= [50 50 50];                  % dB SPL (80 for Lieberman)
noiseLevel=-35;

%% #5 number of channels in the model
%   21-channel model (log spacing)
numChannels=21;
lowestBF=500; 	highestBF= 2000;
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

%   or specify your own channel BFs
numChannels=1;
BFlist=1000;


%% #6 change model parameters

paramChanges={};

% Parameter changes can be used to change one or more model parameters
%  *after* the MAPparams file has been read
% This example declares only one fiber type with a calcium clearance time
% constant of 80e-6 s (HSR fiber) when the probability option is selected.
% paramChanges={'AN_IHCsynapseParams.ANspeedUpFactor=5;', ...
%     'IHCpreSynapseParams.tauCa=86e-6; '};



%% delare 'showMap' options to control graphical output
showMapOptions.printModelParameters=1;   % prints all parameters
showMapOptions.showModelOutput=1;       % plot of all stages
showMapOptions.printFiringRates=1;      % prints stage activity levels
showMapOptions.showACF=0;               % shows SACF (probability only)
showMapOptions.showEfferent=0;          % tracks of AR and MOC
showMapOptions.surfProbability=1;       % 2D plot of HSR response 
showMapOptions.surfSpikes=0;            % 2D plot of spikes histogram
showMapOptions.ICrates=0;               % IC rates by CNtauGk
showMapOptions.PSTHbinwidth=0.001;
showMapOptions.colorbar=0;
showMapOptions.view=[0 90];

% disable certain silly options
if strcmp(AN_spikesOrProbability, 'spikes')
    % avoid nonsensical options
    showMapOptions.surfProbability=0;
    showMapOptions.showACF=0;
end

if strcmp(signalType, 'file')
    % needed for labeling plot
    showMapOptions.fileName=fileName;
else
    showMapOptions.fileName=[];
end

%% Generate stimuli

switch signalType
    case 'tones'
        % Create pure tone stimulus
        dt=1/sampleRate; % seconds
        time=dt: dt: duration;
        amp=10.^(leveldBSPL/20)*28e-6;   % converts to Pascals (peak)
        inputSignal=sin(2*pi*toneFrequency'*time);
        amps=repmat(amp',1,length(time));
        inputSignal=amps.*inputSignal;
        inputSignal=sum(inputSignal, 1);
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

        %% TEN noise input simple or mixed
        [noise sampleRate]=wavread('TEN.wav');
        dt=1/sampleRate;
        noise=noise(:,1);
        targetRMS=20e-6*10^(noiseLevel/20);
        rms=(mean(noise.^2))^0.5;
        amp=targetRMS/rms;
        noise=noise*amp;
        inputSignal=inputSignal+noise(1:length(inputSignal))';
end


%% run the model
tic
fprintf('\n')
disp(['Signal duration= ' num2str(length(inputSignal)/sampleRate)])
disp([num2str(numChannels) ' channel model: ' AN_spikesOrProbability])
disp('Computing ...')

MAP1_14(inputSignal, sampleRate, BFlist, ...
    MAPparamsName, AN_spikesOrProbability, paramChanges);


%% the model run is now complete. Now display the results
UTIL_showMAP(showMapOptions, paramChanges)
figure(97)
title (['tones / noise levels: ' num2str([leveldBSPL noiseLevel])])
if strcmp(signalType,'tones')
    disp(['duration=' num2str(duration)])
    disp(['level=' num2str(leveldBSPL)])
    disp(['toneFrequency=' num2str(toneFrequency)])
    global DRNLParams
    disp(['attenuation factor =' ...
        num2str(DRNLParams.rateToAttenuationFactor, '%5.3f') ])
    disp(['attenuation factor (probability)=' ...
        num2str(DRNLParams.rateToAttenuationFactorProb, '%5.3f') ])
    disp(AN_spikesOrProbability)
end
disp(paramChanges)
toc
path(restorePath)

