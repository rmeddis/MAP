function MAPrunner(MAPparamsName, AN_spikesOrProbability, ...
    signalCharacteristics, paramChanges, showMapOptions)

dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])


%% #3 pure tone, harmonic sequence or speech file input
signalType= signalCharacteristics.type;
sampleRate= signalCharacteristics.sampleRate;
duration=signalCharacteristics.duration;                 % seconds
rampDuration=signalCharacteristics.rampDuration;              % raised cosine ramp (seconds)
beginSilence=signalCharacteristics.beginSilence;               
endSilence=signalCharacteristics.endSilence;                  
toneFrequency= signalCharacteristics.toneFrequency;            % or a pure tone (Hz)
leveldBSPL=signalCharacteristics.leveldBSPL;

BFlist=-1;


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

MAP1_14(inputSignal, sampleRate, BFlist, ...
    MAPparamsName, AN_spikesOrProbability, paramChanges);

%% the model run is now complete. Now display the results
UTIL_showMAP(showMapOptions, paramChanges)

path(restorePath)

