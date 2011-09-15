function demoTwisterSpikes

% MAPdemo runs the MATLAB auditory periphery model (MAP1_14) as far as
%  IC (spikes) with graphical output

restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])


%%  #1 parameter file name
MAPparamsName='Normal';


%% #2 spikes (slow) representation
AN_spikesOrProbability='spikes';


%% #3 speech file input
fileName='twister_44kHz';


%% #4 rms level
% signal details
leveldBSPL=70;                  % dB SPL


%% #5 number of channels in the model
%   21-channel model (log spacing)
numChannels=21;
lowestBF=250; 	highestBF= 8000;
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));


%% #6 change model parameters
paramChanges=[];

%% delare showMap options
showMapOptions=[];  % use defaults

% or (example: show everything including an smoothed SACF output
showMapOptions.printModelParameters=1;
showMapOptions.showModelOutput=1;
showMapOptions.printFiringRates=1;
showMapOptions.showACF=0;
showMapOptions.showEfferent=0;
showMapOptions.surfSpikes=0;
showMapOptions.surfProbability=0;       % 2D plot of HSR response

%% Generate stimuli
[inputSignal sampleRate]=wavread(fileName);
inputSignal(:,1);
targetRMS=20e-6*10^(leveldBSPL/20);
rms=(mean(inputSignal.^2))^0.5;
amp=targetRMS/rms;
inputSignal=inputSignal*amp;


%% run the model
tic

fprintf('\n')
disp(['Signal duration= ' num2str(length(inputSignal)/sampleRate)])
disp([num2str(numChannels) ' channel model'])
disp('Computing ...')

MAP1_14(inputSignal, sampleRate, BFlist, ...
    MAPparamsName, AN_spikesOrProbability, paramChanges);


% the model run is now complete. Now display the results
UTIL_showMAP(showMapOptions, paramChanges)

toc
path(restorePath)

