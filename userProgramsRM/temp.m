function test_speechInNoise

leveldBSPL= 60;                  % dB SPL
leveldBSPLNoise=-55;

paramChanges={};
% no attenuation
paramChanges={'DRNLParams.rateToAttenuationFactorProb = 0.00; '};
% fixed attenuation
% paramChanges={'DRNLParams.rateToAttenuationFactorProb = -0.04;'};
% % dynamic attenuation
paramChanges={'DRNLParams.MOCtauProb =.15;', ...
    'DRNLParams.rateToAttenuationFactorProb = 0.01; '};
% 
fileName='twister_44kHz';
fileName='1o7a_44kHz';

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


%% #3  speech file input

beginSilence=.25;
endSilence=0.25;
noiseRampDuration=0.01;

%% #5 number of channels in the model
%   21-channel model (log spacing)
numChannels=21;
lowestBF=300; 	highestBF= 6000;
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

%   or specify your own channel BFs
% numChannels=1;
% BFlist=1000;



%% delare 'showMap' options to control graphical output
showMapOptions.printModelParameters=1;   % prints all parameters
showMapOptions.showModelOutput=1;       % plot of all stages
showMapOptions.printFiringRates=1;      % prints stage activity levels
showMapOptions.showACF=0;               % shows SACF (probability only)
showMapOptions.showEfferent=1;          % tracks of AR and MOC
showMapOptions.surfProbability=1;       % 2D plot of HSR response
showMapOptions.surfSpikes=0;            % 2D plot of spikes histogram
showMapOptions.ICrates=0;               % IC rates by CNtauGk
showMapOptions.PSTHbinwidth=0.002;

% disable certain silly options
if strcmp(AN_spikesOrProbability, 'spikes')
    % avoid nonsensical options
    showMapOptions.surfProbability=0;
    showMapOptions.showACF=0;
else
    showMapOptions.surfSpikes=0;
end
    % needed for labeling plot
    showMapOptions.fileName=fileName;

%% Generate stimuli

        %% file input simple or mixed
        [inputSignal sampleRate]=wavread(fileName);
        dt=1/sampleRate;
        inputSignal=inputSignal(:,1);
        targetRMS=20e-6*10^(leveldBSPL/20);
        rms=(mean(inputSignal.^2))^0.5;
        amp=targetRMS/rms;
        inputSignal=inputSignal*amp;
        
        % add silences
        intialSilence= zeros(1,round(beginSilence*sampleRate));
        finalSilence= zeros(1,round(endSilence*sampleRate));
        inputSignal= [intialSilence inputSignal' finalSilence];
        
        [inputNoise sampleRateN]=wavread('babble');
        inputNoise=inputNoise(1:length(inputSignal));
       inputNoise=inputNoise(:,1);
        targetRMS=20e-6*10^(leveldBSPLNoise/20);
        rms=(mean(inputNoise.^2))^0.5;
        amp=targetRMS/rms;
        inputNoise=inputNoise*amp;
        time=dt: dt: dt*length(inputNoise);
        rampTime=dt:dt:noiseRampDuration;
        ramp=[0.5*(1+cos(2*pi*rampTime/(2*noiseRampDuration)+pi)) ...
            ones(1,length(time)-length(rampTime))];
        inputNoise=inputNoise'.*ramp;
        inputSignal=inputSignal+inputNoise;
        


%% run the model
tic

fprintf('\n')
disp(['Signal duration= ' num2str(length(inputSignal)/sampleRate)])
disp([num2str(numChannels) ' channel model'])
disp('Computing ...')

MAP1_14(inputSignal, sampleRate, BFlist, ...
    MAPparamsName, AN_spikesOrProbability, paramChanges);


%% the model run is now complete. Now display the results
UTIL_showMAP(showMapOptions, paramChanges)
figure(97), view([-3 82])
title(['speech/ noise: ' num2str([leveldBSPL leveldBSPLNoise])])

disp(['level=' num2str(leveldBSPL)])
disp(['noise level=' num2str(leveldBSPLNoise)])

global DRNLParams
disp(['attenuation factor =' ...
    num2str(DRNLParams.rateToAttenuationFactor, '%5.3f') ])
disp(['attenuation factor (probability)=' ...
    num2str(DRNLParams.rateToAttenuationFactorProb, '%5.3f') ])
disp(AN_spikesOrProbability)
disp(paramChanges)
toc
path(restorePath)

