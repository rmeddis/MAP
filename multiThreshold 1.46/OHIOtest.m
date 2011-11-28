dbstop if error
path=pathdef;
restorePath=path;
addpath (['..' filesep 'utilities']) % model physiology tests

globalStimParams.FS=100000;
doPlot=0;

stim.type='OHIO';
stim.OHIOtype='OHIOrand';
stim.phases='sin';
stim.beginSilence=0;
stim.endSilence=-1;
stim.rampOnDur=.002;
stim.rampOffDur=-1;

% 1. ‘OHIOabs’ paradigm is a baseline procedure for measuring absolute
% thresholds (in dB SPL) of the single tone with 12 frequencies:
%    1   2    3     4    5     6     7     8     9     10     11    12
% 494, 663, 870, 1125, 1442, 1838, 2338, 2957, 3725, 4689, 5866, 7334

% 2. ‘OHIOtemp’ is for measuring thresholds for temporally integrated
% combinations of 2, 4, 8, and 12 tones presented simultaneously.
% In our experiment, we used 4680Hz frequency.

% 3. ‘OHIOspec’ is for measuring thresholds for spectrally integrated
% combinations of 2(7335 and 5866Hz), 4(7334, 5866, 4680, and 3725Hz),
% 8(7334, 5866, 4680, 3725, 2957, 2338, 1838, and
% 1442Hz), and
% 12(all 12 frequencies) tones presented simultaneously.

% 4. ‘OHIOspectemp’ is for measuring thresholds for patterned signals
% differing in both the spectral and temporal domains.
% The frequency conditions are the same as that of ‘OHIOspec’.

% 5. ‘OHIOrand’ is for measuring thresholds for spectrotemporally varying
% signals with random frequency presentation.

nTonesList=[2 4 8 12];
allFreqs=[494, 663, 870, 1125, 1442, 1838, 2338, 2957, 3725, 4689, 5866, 7334];
absThresholds= 50*ones(1,12);   % protem

for nTones=nTonesList
    switch stim.OHIOtype
        case ' OHIOabs'
            % one tone frequency at a time
            stim.frequencies=allFreqs(1);
            stim.amplitudesdB=0;    

        case 'OHIOrand'
            % chose nTones frequencies at random
            x=rand(1,12);
            [sorted idx]=sort(x);
            stim.frequencies=allFreqs(idx(1:nTones));
            stim.amplitudesdB=absThresholds(idx);

        case 'OHIOtemp'
            % 4680 Hz repeated nTones times
            stim.frequencies=4680*ones(1,nTones);
            stim.amplitudesdB=repmat(absThresholds(10),1,nTones);   

        case {'OHIOspect',  'OHIOspectemp'}
            % nTones frequencies either simulataneously or sequentially
            switch nTones
                case 2
                    stim.frequencies=[7335 5866];
                    idx=[12 11];
                    stim.amplitudesdB=absThresholds(idx);    
                case 4
                    stim.frequencies=[7334, 5866, 4680, 3725];
                    idx=[12:-1:9 ];
                    stim.amplitudesdB=absThresholds(idx);    
                case 8
                    stim.frequencies=...
                        [7334, 5866, 4680, 3725, 2957, 2338, 1838, 1442];
                    idx=[12:-1:5 ];
                    stim.amplitudesdB=absThresholds(idx);    
                case 12
                    stim.frequencies=allFreqs;
                    idx=[12:-1:1 ];
                    stim.amplitudesdB=absThresholds(idx);    
            end
    end

    switch stim.OHIOtype
        case {'OHIOabs', 'OHIOspect'}
            stim.toneDuration=.02;
            globalStimParams.overallDuration=stim.toneDuration;
        otherwise
            stim.toneDuration=nTones*0.02;
            globalStimParams.overallDuration=stim.toneDuration;
    end

    disp(num2str(stim.frequencies))

    [audio, msg]=stimulusCreate(globalStimParams, stim, doPlot);
    wavplay(audio,globalStimParams.FS)
end
path=restorePath;
