function pitchModel_RM
% Modification of testMAP_14 to replicate the pitch model published
%  in JASA 2006.
%
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

% Pitch model modification here
global ICrate           % used to collect rate profile from showMAP temporary
rates=[]; F0count=0;

% F0s=[150 200 250];          % fundamental frequency
% harmonics= 3:5;

% F0s=[3000];          % fundamental frequency
F0s=50:5:1000;
harmonics= 1;
% F0s=150;
for F0=F0s
    F0count=F0count+1;
    

%%  #1 parameter file name
MAPparamsName='Normal';


%% #2 probability (fast) or spikes (slow) representation
AN_spikesOrProbability='spikes';

% or
% NB probabilities are not corrected for refractory effects
% AN_spikesOrProbability='probability';


%% #3 pure tone, harmonic sequence or speech file input
signalType= 'tones';
sampleRate= 50000;
duration=0.50;                 % seconds
% toneFrequency= 1000;            % or a pure tone (Hz8

% F0=210;
toneFrequency= F0*harmonics;  % harmonic sequence (Hz)

rampDuration=.005;              % raised cosine ramp (seconds)

% or

% signalType= 'file';
% fileName='twister_44kHz';


%% #4 rms level
% signal details
leveldBSPL= 50;                  % dB SPL


%% #5 number of channels in the model
%   21-channel model (log spacing)
numChannels=21;
lowestBF=250; 	highestBF= 8000;
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

%   or specify your own channel BFs
% numChannels=1;
BFlist=toneFrequency;
% BFlist=500;


%% #6 change model parameters
paramChanges=[];

% or
% Parameter changes can be used to change one or more model parameters
%  *after* the MAPparams file has been read
% This example declares only one fiber type with a calcium clearance time
% constant of 80e-6 s (HSR fiber) when the probability option is selected.

% paramChanges={'AN_IHCsynapseParams.ANspeedUpFactor=5;', ...
%     'IHCpreSynapseParams.tauCa=86e-6;'};

% paramChanges={'DRNLParams.rateToAttenuationFactorProb = 0;'};

% paramChanges={'IHCpreSynapseParams.tauCa=86e-6;',
%     'AN_IHCsynapseParams.numFibers=	1000;'};

% fixed MOC attenuation(using negative factor)
% paramChanges={'DRNLParams.rateToAttenuationFactorProb=-0.005;'};

% slow the CN chopping rate
% paramChanges={'IHCpreSynapseParams.tauCa= 70e-6;'...'
%     'MacGregorMultiParams.tauGk=	[0.75e-3:.0001 : 3e-3];'...
%     ' MacGregorParams.dendriteLPfreq=4000;'...
%     'MacGregorParams.tauGk=	1e-4;'...
%     'MacGregorParams.currentPerSpike=220e-8;'...
% };
paramChanges={...
    'MacGregorMultiParams.currentPerSpike=25e-9;'...
    'MacGregorMultiParams.tauGk=	[0.1e-3:.00005 : 1e-3];'...
'MacGregorParams.currentPerSpike=40e-9;'...
};

%% delare 'showMap' options to control graphical output

showMapOptions.printModelParameters=0;   % prints all parameters
showMapOptions.showModelOutput=1;       % plot of all stages
showMapOptions.printFiringRates=1;      % prints stage activity levels
showMapOptions.showACF=0;               % shows SACF (probability only)
showMapOptions.showEfferent=0;          % tracks of AR and MOC
showMapOptions.surfProbability=0;       % 2D plot of HSR response 
showMapOptions.surfSpikes=0;            % 2D plot of spikes histogram
showMapOptions.ICrates=1;               % IC rates by CNtauGk

% disable certain silly options
if strcmp(AN_spikesOrProbability, 'spikes')
    % avoid nonsensical options
    showMapOptions.surfProbability=0;
    showMapOptions.showACF=0;
else
    showMapOptions.surfSpikes=0;
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
        inputSignal=createMultiTone(sampleRate, toneFrequency, ...
            leveldBSPL, duration, rampDuration);

    case 'file'
        %% file input simple or mixed
        [inputSignal sampleRate]=wavread(fileName);
        dt=1/sampleRate;
        inputSignal=inputSignal(:,1);
        targetRMS=20e-6*10^(leveldBSPL/20);
        rms=(mean(inputSignal.^2))^0.5;
        amp=targetRMS/rms;
        inputSignal=inputSignal*amp;
        silence= zeros(1,round(0.1/dt));
        inputSignal= [silence inputSignal' silence];
end


%% run the model
tic

fprintf('\n')
disp(['Signal duration= ' num2str(length(inputSignal)/sampleRate)])
disp([num2str(numChannels) ' channel model'])
disp([num2str(F0) ' F0'])
disp('Computing ...')

MAP1_14(inputSignal, sampleRate, BFlist, ...
    MAPparamsName, AN_spikesOrProbability, paramChanges);


% the model run is now complete. Now display the results
UTIL_showMAP(showMapOptions, paramChanges)

%% pitch model Collect and analyse data
% ICrate is global and computed in showMAP
%  a vector of 'stage4' rates; one value for each tauCNGk
rates=[rates; ICrate];
figure(92), imagesc(rates)
ylabel ('F0 no'), xlabel('tauGk')
% figure(92), plot(rates), ylim([0 inf])

h=figure(99); CNmovie(F0count)=getframe(h);
figure(91), plot(rates'),ylim([0 inf])
pause (0.1)
path(restorePath)

end
%% show results
toc
figure(91), plot(F0s,rates'), xlabel('F0'), ylabel('rate'),ylim([0 inf])
% figure(99),clf,movie(CNmovie,1,4)


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
silence= zeros(1,round(0.005/dt));
inputSignal= [silence inputSignal silence];

