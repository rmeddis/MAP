function runMAP1_14
% runMAP1_14 is a general purpose test routine that can be adjusted to
% explore different applications of MAP1_14
%
% It is also designed as 'starter' code for building new applications.
%
% A range of options are supplied in the early part of the program
%
% #1
% Identify the file (in 'MAPparamsName') containing the model parameters
%
% #2
% Identify the kind of model required (in 'AN_spikesOrProbability').
%  A full brainstem model ('spikes') can be computed or a shorter model
%  ('probability') that computes only so far as the auditory nerve
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
% #6
% Last minute changes to the model parameters can be made using
%  a cell array of strings, 'paramChanges'.

dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])

%%  #1 parameter file name
MAPparamsName='Normal';


%% #2 probability (fast) or spikes (slow) representation: select one
% AN_spikesOrProbability='spikes';
%   or
AN_spikesOrProbability='probability';


%% #3 A. pure tone, B. harmonic sequence or C. speech file input
% comment out unwanted code

% A. tone
sampleRate= 95000;
signalType= 'tones';
toneFrequency= 1000;            % or a pure tone (Hz)
duration=0.500;                 % seconds
beginSilence=0.010;
endSilence=0.020;
rampDuration=.005;              % raised cosine ramp (seconds)

%   or
% B. harmonic tone (Hz) - useful to demonstrate a broadband sound
% sampleRate= 44100;
% signalType= 'tones';
% toneFrequency= F0:F0:8000;
% duration=0.500;                 % seconds
% beginSilence=0.250;
% endSilence=0.250;
% F0=210;
% rampDuration=.005;              % raised cosine ramp (seconds)

%   or
% C. signalType= 'file';
% fileName='twister_44kHz';

%% #4 rms level
% signal details
leveldBSPL= 80;                  % dB SPL (80 for Lieberman)

%% #5 number of channels in the model
%   21-channel model (log spacing)
numChannels=21;
lowestBF=100; 	highestBF= 6000;
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

%   or specify your own channel BFs
% numChannels=1;
% BFlist=toneFrequency;


%% #6 change model parameters

paramChanges={};    % no changes

% Parameter changes can be used to change one or more model parameters
%  *after* the MAPparams file has been read
%  Each string must have the same format as the corresponding line in the
%   file identified in 'MAPparamsName'
% This example declares only one fiber type with a calcium clearance time
% constant of 80e-6 s (HSR fiber) when the probability option is selected.
%   paramChanges={'AN_IHCsynapseParams.ANspeedUpFactor=5;', ...
%     'IHCpreSynapseParams.tauCa=86e-6; '};
paramChanges={ 'IHCpreSynapseParams.tauCa=86e-6; ',...
    'DRNLParams.rateToAttenuationFactorProb = 0;'};



%% delare 'showMap' options to control graphical output
% see UTIL_showMAP for more options
showMapOptions.printModelParameters=1;   % prints all parameters
showMapOptions.showModelOutput=0;       % plot of all stages
showMapOptions.printFiringRates=1;      % prints stage activity levels
showMapOptions.showEfferent=0;          % tracks of AR and MOC
showMapOptions.surfAN=1;       % 2D plot of HSR response

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
        intialSilence= zeros(1,round(0.1/dt));
        finalSilence= zeros(1,round(0.2/dt));
        inputSignal= [intialSilence inputSignal' finalSilence];
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
UTIL_showMAP(showMapOptions)

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
disp('paramChanges')
for i=1:length(paramChanges)
    disp(paramChanges{i})
end

toc
path(restorePath)

