function test_MAP1_14
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
% Indentify the channels in terms of their best frequencies in the vector
%  BFlist.
%
% Last minute changes to the parameters fetched earlier can be made using
%  the cell array of strings 'paramChanges'.
%  Each string must have the same format as the corresponding line in the
%  file identified in 'MAPparamsName'
%
% When the demonstration is satisfactory, freeze it by renaming it <demoxx>


%%  #1 parameter file name
MAPparamsName='Normal';


%% #2 probability (fast) or spikes (slow) representation
AN_spikesOrProbability='spikes';
% or
AN_spikesOrProbability='probability';


%% #3 pure tone, harmonic sequence or speech file input
signalType= 'tones';
duration=0.100;                 % seconds
% duration=0.020;                 % seconds
sampleRate= 64000;
% toneFrequency= 250:250:8000;    % harmonic sequence (Hz)
toneFrequency= 2000;            % or a pure tone (Hz8

rampDuration=.005;              % seconds

% or
signalType= 'file';
fileName='twister_44kHz';
% fileName='new-da-44khz';

% mix with an optional second file?
mixerFile=[];
%or
mixerFile='babble';
leveldBSPL2=-60;

%% #4 rms level
% signal details
leveldBSPL= 60;                  % dB SPL


%% #5 number of channels in the model
%   21-channel model (log spacing)
numChannels=21;
lowestBF=250; 	highestBF= 8000; 
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

%   or specify your own channel BFs
BFlist=toneFrequency;


%% #6 change model parameters
paramChanges=[];

% or
% Parameter changes can be used to change one or more model parameters
%  *after* the MAPparams file has been read
% This example declares only one fiber type with a calcium clearance time
% constant of 80e-6 s (HSR fiber) when the probability option is selected.
% paramChanges={'AN_IHCsynapseParams.ANspeedUpFactor=5;', ...
%     'IHCpreSynapseParams.tauCa=86e-6;'};
% paramChanges={'AN_IHCsynapseParams.ANspeedUpFactor=5;', ...
%     'DRNLParams.rateToAttenuationFactorProb = 0;'};

%% delare showMap options
showMapOptions=[];  % use defaults

% or (example: show everything including an smoothed SACF output
showMapOptions.showModelParameters=1;
showMapOptions.showModelOutput=1;
showMapOptions.printFiringRates=1;
showMapOptions.showACF=0;
showMapOptions.showEfferent=1;
if strcmp(AN_spikesOrProbability, 'probability')
    showMapOptions.surfProbability=1;
end
if strcmp(signalType, 'file')
    showMapOptions.fileName=fileName;
end

%% Generate stimuli

dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'])
switch signalType
    case 'tones'
        inputSignal=createMultiTone(sampleRate, toneFrequency, ...
            leveldBSPL, duration, rampDuration);
        
    case 'file'
        %% file input simple or mixed
        [inputSignal sampleRate]=wavread(fileName);       
        inputSignal=inputSignal(:,1);
        targetRMS=20e-6*10^(leveldBSPL/20);
        rms=(mean(inputSignal.^2))^0.5;
        amp=targetRMS/rms;
        inputSignal=inputSignal*amp;
        if ~isempty(mixerFile)
            [inputSignal2 sampleRate]=wavread(mixerFile);
            inputSignal2=inputSignal2(:,1);
            [r c]=size(inputSignal);
            inputSignal2=inputSignal2(1:r);
            targetRMS=20e-6*10^(leveldBSPL2/20);
            rms=(mean(inputSignal2.^2))^0.5;
            amp=targetRMS/rms;
            inputSignal2=inputSignal2*amp;
        inputSignal=inputSignal+inputSignal2;
        end
end


%% run the model
tic

fprintf('\n')
disp(['Signal duration= ' num2str(length(inputSignal)/sampleRate)])
disp([num2str(numChannels) ' channel model'])
disp('Computing ...')

restorePath=path;
addpath (['..' filesep 'MAP'])

MAP1_14(inputSignal, sampleRate, BFlist, ...
    MAPparamsName, AN_spikesOrProbability, paramChanges);
path(restorePath)
toc

% the model run is now complete. Now display the results
showMAP(showMapOptions)
for i=1:length(paramChanges)
disp(paramChanges{i})
end

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
inputSignal= [silence inputSignal silence];

