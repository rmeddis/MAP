function test_MAP1_14RAM

global dt  DRNLoutput 
global  DRNLParams TMoutput
global peakOutputDisp alignedOutput toneFrequencyList levels leveldBSPL logOutput


dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])
figure(2), clf

%% # BFlist is the BF of the filter to be assessed
BFlist=8000;
numChannels=1;

% define probe frequencies
numFs=6;
lowestF=BFlist/5; 	highestF= BFlist*1.4;
lowestF=4000; 	highestF= 12600;
toneFrequencyList=round(logspace(log10(lowestF), log10(highestF), numFs));

%%  # parameter file name
MAPparamsName='Normal';

%% # probability representation
AN_spikesOrProbability='probability';

%% # tone duration
sampleRate= 100000;
duration=0.0200;                 % seconds
rampDuration=.005;              % raised cosine ramp (seconds)
beginSilence=0.050;
endSilence=0.050;

%% # levels
% probe details
levels=0:10:80;
% levels=80;

%% # change model parameters
% Parameter changes can be used to change one or more model parameters
%  *after* the MAPparams file has been read

% adjust linear path gain (previously g=100)
% switch off all efferent effects
paramChanges={...
    'DRNLParams.rateToAttenuationFactorProb = 0.00; ',...
    'OMEParams.rateToAttenuationFactorProb=0.0;', ...
    'DRNLParams.ctBMdB = -10;'...
    'DRNLParams.g=200;'...
    'DRNLParams.linCFs=6500;'...
    'DRNLParams.linBWs=2000;'...
    };
%     'DRNLParams.a=0;'...


%% delare 'showMap' options to control graphical output
showMapOptions.printModelParameters=0;   % prints all parameters
showMapOptions.showModelOutput=1;       % plot of all stages
showMapOptions.printFiringRates=0;      % prints stage activity levels
showMapOptions.showACF=0;               % shows SACF (probability only)
showMapOptions.showEfferent=0;          % tracks of AR and MOC
showMapOptions.surfProbability=0;       % 2D plot of HSR response
showMapOptions.surfSpikes=0;            % 2D plot of spikes histogram
showMapOptions.ICrates=0;               % IC rates by CNtauGk
showMapOptions.fileName=[];


%% now vary level and frequency
peakOutputDisp=zeros(length(levels),length(toneFrequencyList));
peakStapesDisp=zeros(length(levels),length(toneFrequencyList));
levelNo=0;
for leveldBSPL=levels
    levelNo=levelNo+1;
disp(['level: ' num2str(leveldBSPL)])
    freqNo=0;
    for toneFrequency=toneFrequencyList
        freqNo=freqNo+1;

        %% Generate stimuli
        dt=1/sampleRate;
        time=dt: dt: duration;
        inputSignal=sum(sin(2*pi*toneFrequency'*time), 1);
        amp=10^(leveldBSPL/20)*28e-6;   % converts to Pascals (peak)
        inputSignal=amp*inputSignal;
        % apply ramps
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
%         fprintf('\n')
%         disp([num2str(numChannels) ' channel model: ' AN_spikesOrProbability])
%         disp('Computing ...')

        MAP1_14(inputSignal, sampleRate, BFlist, ...
            MAPparamsName, AN_spikesOrProbability, paramChanges);

        peakOutputDisp(levelNo,freqNo)=max(DRNLoutput);
        peakStapesDisp(levelNo,freqNo)=max(TMoutput);

        %% the model run is complete. Now display the results
%         UTIL_showMAP(showMapOptions, paramChanges)
    end % probe frequencies
figure(2), loglog(toneFrequencyList, peakOutputDisp), hold on
xlabel('frequency')
ylabel('peak DRNL displacement (m)')
    
end  % levels
%% alignmed plot


x=peakOutputDisp./peakStapesDisp;
logOutput=20*log10(x.*repmat((2*pi*toneFrequencyList),length(levels),1));
alignedOutput=logOutput-repmat(logOutput(:,1), 1, length(toneFrequencyList));
figure(3), clf, semilogx(toneFrequencyList, alignedOutput), hold on
ylim([-20 80])
xlim([2000 20000])
% legend(num2str(levels'))
xlabel('frequency')
ylabel('peak DRNL velocity (dB)')
title(['level ' num2str(leveldBSPL) ' dB SPL'])

disp(['level=' num2str(leveldBSPL)])
disp(['toneFrequency=' num2str(toneFrequency)])

disp(['attenuation factor =' ...
    num2str(DRNLParams.rateToAttenuationFactor, '%5.3f') ])
disp(['attenuation factor (probability)=' ...
    num2str(DRNLParams.rateToAttenuationFactorProb, '%5.3f') ])
disp(AN_spikesOrProbability)
disp(paramChanges)

%     UTIL_printTabTable([toneFrequencyList' peakOutputDisp'])

path(restorePath)

