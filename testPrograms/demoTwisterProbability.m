function demoTwisterProbability

% MAPdemo runs the MATLAB auditory periphery model (MAP1_14) as far as
%  the AN (probabilities) or IC (spikes) with graphical output

dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])

%%  #1 parameter file name
MAPparamsName='Normal';


%% #2 probability (fast) 
AN_spikesOrProbability='probability';


%% #3  speech file input
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


%% the model run is finished. Now display the results

%% delare showMap options
showMapOptions.printModelParameters=1;
showMapOptions.showModelOutput=1;
showMapOptions.printFiringRates=1;
showMapOptions.showACF=0;
showMapOptions.showEfferent=0;
showMapOptions.surfProbability=1;       % 2D plot of HSR response 

UTIL_showMAP(showMapOptions, paramChanges)

toc
path(restorePath)

