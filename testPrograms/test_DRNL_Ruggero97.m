function test_DRNL_Ruggero97
% test_DRNL_Ruggero97 attempts to match Ruggero's (1992 and 1997)
%  iso-intensity data by fiddling with the parameters

% # BF is the BF of the filter to be assessed
BF=9000;

% # test frequencies. check that BF is one of them
%    copy Ruggero's test tones as fara as possible
numFs=6; lowestF=4000; 	highestF= 11030;
toneFrequencyList=round(logspace(log10(lowestF), log10(highestF), numFs));

%  # parameter file name. this is the base set of parameters
MAPparamsName='Normal';

% # probability representation (not directly relevant here as only
%    the DRNL output is used
AN_spikesOrProbability='probability';

% # tone characteristics
sampleRate= 100000;
duration=0.0200;                % Ruggero uses 5, 10, 25 ms tones
rampDuration=0.0015;              % raised cosine ramp (seconds)
beginSilence=0.050;
endSilence=0.020;

% # levels
levels=[3 10:10:80];

%% # change model parameters
% Parameter changes can be used to change one or more model parameters
%  *after* the MAPparams file has been read

paramChanges={};

paramChanges={...
    'DRNLParams.ctBMdB = -20;'...
    'DRNLParams.g=1000;'...
    };

% delare 'showMap' options to control graphical output
% (not needed but might be useful
% showMapOptions.printModelParameters=0;   % prints all parameters
% showMapOptions.showModelOutput=1;       % plot of all stages
% showMapOptions.printFiringRates=0;      % prints stage activity levels
% showMapOptions.showACF=0;               % shows SACF (probability only)
% showMapOptions.showEfferent=0;          % tracks of AR and MOC
% showMapOptions.surfAN=0;       % 2D plot of HSR response
% showMapOptions.surfSpikes=0;            % 2D plot of spikes histogram
% showMapOptions.ICrates=0;               % IC rates by CNtauGk
% showMapOptions.fileName=[];

% run the program
global dt  DRNLoutput DRNLParams TMoutput
dbstop if error
restorePath=path; 
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])
figure(2), clf,figure(3), clf,figure(4), clf,

peakOutputDisp=zeros(length(levels),length(toneFrequencyList));
peakStapesDisp=zeros(length(levels),length(toneFrequencyList));

%% now vary level and frequency while measuring the response
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
    
        MAP1_14(inputSignal, sampleRate, BF, ...
            MAPparamsName, AN_spikesOrProbability, paramChanges);

        peakOutputDisp(levelNo,freqNo)=max(DRNLoutput);
        peakStapesDisp(levelNo,freqNo)=max(TMoutput);

        % the model run is complete. Now display the results if debugging
        % UTIL_showMAP(showMapOptions)
    end % probe frequencies
    
    % monitor progress
    figure(2), loglog(toneFrequencyList, peakOutputDisp), hold on
    xlabel('frequency')
    ylabel('peak DRNL displacement (m)')
    title ('DRNL uncorrected displacement')
end  % levels
figure(2),legend(num2str(toneFrequencyList'),'location','northwest')

% convert from model BM displacement to Ruggero's velocity
DRNLvelocity= peakOutputDisp ...
    .*repmat(2*pi*toneFrequencyList,length(levels),1);
% convert from model stapes displacement to Ruggero's velocity
stapesVelocity= peakStapesDisp ...
    .*repmat(2*pi*toneFrequencyList,length(levels),1);
% velocity gain is increased velocity attributable to the DRNL
DRNLvelocityGain = DRNLvelocity./stapesVelocity;
DRNLvelocityGain_dB= 20*log10(DRNLvelocityGain );
% iso-intensity equates all functions at the same outpu for the lowest
%  velocity tested
isoIntensityDRNLvel_dB= DRNLvelocityGain_dB- ...
    repmat(DRNLvelocityGain_dB(:,1),1,length(toneFrequencyList));

% displays
% iso-velocity function (nrmalised by stapes velocity)
figure(3), clf, semilogx(toneFrequencyList, isoIntensityDRNLvel_dB)
ylim([-10 50])
xlim([3000 20000])
xlabel('frequency (Hz)')
set(gca,'Xtick', [1000 ],'xticklabel',{'1000'})

ylabel('peak DRNL velocity gain (dB)')
title(['CF= ' num2str(BF) ' Hz'])
legend(num2str(levels'),'location','northwest')

% velocity I/O function
figure(4), clf, semilogy(levels, (DRNLvelocity*1e6)'), hold on
ylim([5e0 1.2e4])
xlim([0 80])
xlabel('level')
ylabel(' velocity (microm/ s)')
title(['CF= ' num2str(BF) ' Hz'])
legend(num2str(toneFrequencyList'),'location','northwest')

% command window reports
disp(''), disp('***')
disp(AN_spikesOrProbability)

% DRNL parameter set
UTIL_showStructureSummary(DRNLParams, 'DRNLParams', 10)

% stimulus characteristics
disp(['CF=' num2str(BF)])
disp(['tone Duration=' num2str(rampDuration)])
disp(['ramp Duration=' num2str(duration)])

% particular parameter changes used on this run
for i=1:length(paramChanges)
    disp(paramChanges{i})
end

% if required dump one or more matrices in tab spaced format
%  (suitable for pasting directly into EXCEL.
%     UTIL_printTabTable([toneFrequencyList' peakOutputDisp'])

% leave the path as you found it
path(restorePath)


