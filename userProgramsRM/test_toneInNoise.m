function test_toneInNoise
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

global dt dtSpikes  savedBFlist saveAN_spikesOrProbability saveMAPparamsName...
    savedInputSignal OMEextEarPressure TMoutput OMEoutput ARattenuation ...
    DRNLoutput IHC_cilia_output IHCrestingCiliaCond IHCrestingV...
    IHCoutput ANprobRateOutput ANoutput savePavailable ANtauCas  ...
    CNtauGk CNoutput  ICoutput ICmembraneOutput ICfiberTypeRates ...
    MOCattenuation
global OMEParams DRNLParams IHC_cilia_RPParams IHCpreSynapseParams
global AN_IHCsynapseParams MacGregorParams MacGregorMultiParams
global ICrate


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
sampleRate= 44100;          % must agree with noise
duration=0.010;                 % seconds
rampDuration=.001;              % raised cosine ramp (seconds)
beginSilence=0.010;               
endSilence=0.010;                  
toneFrequency= 2000;            % or a pure tone (Hz)

%   or
% harmonic sequence (Hz)
% F0=210;
% toneFrequency= F0:F0:8000;    

%   or
% signalType= 'file';
% fileName='twister_44kHz';



%% #4 rms level
% signal details
leveldBSPL= 70;                  % dB SPL (80 for Lieberman)
leveldBSPLNoise=30;

%% #5 number of channels in the model
%   21-channel model (log spacing)
numChannels=21;
lowestBF=250; 	highestBF= 8000;
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

% %   or specify your own channel BFs
% numChannels=1;
% BFlist=toneFrequency;


%% #6 change model parameters

paramChanges={};

% Parameter changes can be used to change one or more model parameters
%  *after* the MAPparams file has been read
% This example declares only one fiber type with a calcium clearance time
% constant of 80e-6 s (HSR fiber) when the probability option is selected.
% paramChanges={'AN_IHCsynapseParams.ANspeedUpFactor=5;', ...
%     'IHCpreSynapseParams.tauCa=86e-6; '};
paramChanges={'DRNLParams.MOCtauProb =.25;', ...
    'DRNLParams.rateToAttenuationFactorProb = 0.02; '};
% paramChanges={'DRNLParams.MOCtauProb =.15;', ...
%     'DRNLParams.rateToAttenuationFactorProb = 0.00; '};


%% delare 'showMap' options to control graphical output
showMapOptions.printModelParameters=1;   % prints all parameters
showMapOptions.showModelOutput=1;       % plot of all stages
showMapOptions.printFiringRates=1;      % prints stage activity levels
showMapOptions.showACF=0;               % shows SACF (probability only)
showMapOptions.showEfferent=1;          % tracks of AR and MOC
showMapOptions.surfProbability=1;       % 2D plot of HSR response 
showMapOptions.surfSpikes=1;            % 2D plot of spikes histogram
showMapOptions.ICrates=0;               % IC rates by CNtauGk

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
        
%         [inputNoise sampleRateN]=wavread('babble');
        [inputNoise sampleRateN]=wavread('white noise');
        inputNoise=inputNoise(1:length(inputSignal));
        inputNoise=inputNoise(:,1);
        targetRMS=20e-6*10^(leveldBSPLNoise/20);
        rms=(mean(inputNoise.^2))^0.5;
        amp=targetRMS/rms;
        inputNoise=inputNoise*amp;
        inputSignal=inputSignal+inputNoise';
        figure(2), subplot(2,1,1)
        time=dt:dt:dt*length(inputSignal);
        plot(time,inputSignal,'k')


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
UTIL_showMAP(showMapOptions, paramChanges)

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

figure(2), subplot(2,1,2), plot(ANprobRateOutput(13+21,:))
disp([ 'peak channel 13: ' num2str(max(ANprobRateOutput(13+21,:)))])

for i=1:length(paramChanges)
    disp(paramChanges{i})
end

path(restorePath)

