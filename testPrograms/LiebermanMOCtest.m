function LiebermanMOCtest

% In test_MAP1_14.m set these parameters

% AN_spikesOrProbability='spikes';

% global dt ANdt  savedBFlist saveAN_spikesOrProbability saveMAPparamsName...
%     savedInputSignal OMEextEarPressure TMoutput OMEoutput ARattenuation ...
%     DRNLoutput IHC_cilia_output IHCrestingCiliaCond IHCrestingV...
%     IHCoutput ANprobRateOutput ANoutput savePavailable ANtauCas  ...
%     CNtauGk CNoutput  ICoutput ICmembraneOutput ICfiberTypeRates ...
%     MOCattenuation
global dt ANdt   saveAN_spikesOrProbability ANprobRateOutput ICoutput    

global DRNLParams

LibermanData=[
    2	0.2;
2.1	0.19;2.2	0.18;2.3	0.18;2.4	0.16;2.5	0.15;2.6	0.15;2.7 0.15;
2.8	0.12;2.9	0.12;3	0.1;3.1	0.1;3.2	0.05;3.3	0.05;3.4	0;3.5	-0.1;
3.6	-0.4;3.7	-1.2;3.8	-1.6;3.9	-1.8;4	-1.85;4.1	-2;4.2	-2.05;
4.3	-2.05;4.4	-2.15;4.5	-2.2;4.6	-2.15;4.7	-2.1;4.8	-2.15;4.9 -2.2;
5	-2.1;5.1	-2.1;5.2	-2.25;5.3	-2.1;5.4	-2.15;5.5	-2.1;5.6 -2.15;
5.7	-2.1;5.8	-2.2;5.9	-2.05;6	-2.15;6.1	-2.05;6.2 -2;6.3 -2.2;6.4 -2.1;
6.5	-2.05;6.6	-2.05;6.7	-2.05;6.8 -2.2;6.9 -2.1;7	-2.05;7.1 -2.05;7.2	-0.7;
7.3	-0.1;7.4	0;7.5	0.1;7.6	0.2;7.7	0.35;7.8	0.2;7.9	0.15;8	0.2;8.1	0.15;8.2	0.15;
8.3	0.15;8.4	0.12;8.5	0.1;8.6	0.09;8.7	0.08;8.8	0.07;8.9	0.06;9	0.05;
];

restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])

%%  #1 parameter file name
MAPparamsName='Normal';


%% #2 probability (fast) or spikes (slow) representation
AN_spikesOrProbability='spikes';
%   or
% AN_spikesOrProbability='probability';


%% #3 pure tone, harmonic sequence or speech file input
signalType= 'tones';
sampleRate= 50000;
rampDuration=.005;              % raised cosine ramp (seconds)
toneFrequency= 1000;            % or a pure tone (Hz)
duration=3.6;                   % Lieberman test
beginSilence=1;                 % 1 for Lieberman
endSilence=1;                   % 1 for Lieberman

%% #4 rms level
% signal details
leveldBSPL= 80;                  % dB SPL (80 for Lieberman)


%% #5 number of channels in the model

numChannels=1;
BFlist=toneFrequency;


%% #6 change model parameters
paramChanges={};

%% delare 'showMap' options to control graphical output
showMapOptions.printModelParameters=1;   % prints all parameters
showMapOptions.showModelOutput=1;       % plot of all stages
showMapOptions.printFiringRates=1;      % prints stage activity levels
showMapOptions.showACF=0;               % shows SACF (probability only)
showMapOptions.showEfferent=1;          % tracks of AR and MOC
showMapOptions.surfProbability=1;       % 2D plot of HSR response 
showMapOptions.surfSpikes=0;            % 2D plot of spikes histogram
showMapOptions.ICrates=0;               % IC rates by CNtauGk

%% Generate stimuli

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
    disp(['attenuation factor =' ...
        num2str(DRNLParams.rateToAttenuationFactor, '%5.3f') ])
    disp(['attenuation factor (probability)=' ...
        num2str(DRNLParams.rateToAttenuationFactorProb, '%5.3f') ])
    disp(AN_spikesOrProbability)
end
disp(paramChanges)



%% superimpose Lieberman data

global DRNLParams
% scale up DPOAE results to a maximum of whatever MAP_14 uses as its
% maximum
scalar=-DRNLParams.minMOCattenuationdB/min(LibermanData(:,2));

figure(98), subplot(2,1,2), hold on
plot(LibermanData(:,1)-2.5,scalar*LibermanData(:,2),'r:')

PSTHbinwidth=0.001;

if strcmp(saveAN_spikesOrProbability,'probability')
    PSTH=UTIL_PSTHmaker...
        (ANprobRateOutput(2,:), ANdt, PSTHbinwidth)*ANdt/PSTHbinwidth;
else
    PSTH=UTIL_PSTHmaker(ICoutput(2,:), dt, PSTHbinwidth)*dt/PSTHbinwidth;
end

time=PSTHbinwidth:PSTHbinwidth:PSTHbinwidth*length(PSTH);
figure(89)
plot(time, PSTH)
set(gcf,'name', 'Lieberman')
title(saveAN_spikesOrProbability)

toc
path(restorePath)

% figure(88), plot(MOCattenuation)
