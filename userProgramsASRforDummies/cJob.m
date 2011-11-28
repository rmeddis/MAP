%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef cJob
    %CJOB Responsible for conversion of wav files into recogniser features
    %   Please see the documentation located in a separate file for further
    %   information.
    
    %%  *********************************************************
    %  properties                      _   _
    %                                 | | (_)
    %  _ __  _ __ ___  _ __   ___ _ __| |_ _  ___  ___
    % | '_ \| '__/ _ \| '_ \ / _ \ '__| __| |/ _ \/ __|
    % | |_) | | | (_) | |_) |  __/ |  | |_| |  __/\__ \
    % | .__/|_|  \___/| .__/ \___|_|   \__|_|\___||___/
    % | |             | |
    % |_|             |_|
    %************************************************************
    
    %% **********************************************************
    % Public properties - can be set by user
    %************************************************************
    properties(Access = public)
        wavFolder
        opFolder
        noiseFolder
        
        wavList
        todoStatus
        
        
        participant        = 'Normal';%'DEMO2_multiSpont';
        noiseName          = '8TalkerBabble';
        numWavs            =  5;
        noiseLevToUse      = -200;
        speechLevToUse     =  50;
        
        speechLevStd       = 0;
        noiseLevStd        = 0;
        freezeNoise        = 0;
        meanSNR            = 20;
        speechDist         = 'None';
        noiseDist          = 'None';
        
        noisePreDur        = 0.0;
        noisePostDur       = 0.0;
        truncateDur        = 0.0; %Dr. RF used 0.550
        
        currentSpeechLevel
        currentNoiseLevel
        
        
        useSpectrogram = false;
        numCoeff = 9;
        
        %************************************************************
        % plotting handles
        %************************************************************
        probHaxes   = [];
        probHaxesSM = [];
        sacfHaxes   = [];
        sacfHaxesSM = [];
        featHaxes   = [];
        reconHaxes  = [];
        
        %************************************************************
        % SACF params
        %************************************************************
        useSACF             = false;
        SACFacfTau          = 2; % > 1 = Wiegrebe mode
        SACFnBins           = 128;
        SACFminLag          = 1 / 4000;
        SACFmaxLag          = 1 / 50;
        SACFlambda          = 10e-3;
        
        %************************************************************
        % MAP params
        %************************************************************
        MAProot                 = fullfile('..');
        MAPplotGraphs           = 0;
        MAPDRNLSave             = 0;
        
        MAPopLSR                = 0;
        MAPopMSR                = 0;
        MAPopHSR                = 1;
        
        MAPparamChanges = {};
        
        %************************************************************
        % HTK stuff - writing to HTK recogniser format
        %************************************************************
        frameshift              = 10;       % shift between frames (ms)
        sampPeriodFromMsFactor  = 10000;    % appropriate for 10ms frame rate
        paramKind               = 9;        % HTK USER format for user-defined features (p73 of HTK book)
        
        removeEnergyStatic      = false;
        doCMN                   = false; % Cepstral mean normalisation
        
        %************************************************************
        % Portable EssexAid params
        %************************************************************
        bwOct = 1/1;
        filterOrder  = 2;
        
        mainGain = [ 1;    1;    1;    1;    1];     % gain in linear units
        TCdBO    = [40;   40;   40;   40;   40];      %Compression thresholds (in dB OUTPUT from 2nd filt)
        TMdBO    = [10;   10;   10;   10;   10];      %MOC thresholds (in dB OUTPUT from 2nd filt)
        DRNLc    = [ 0.2;  0.2;  0.2;  0.2;  0.2;]
        
        ARtau  = 0.03;            % decay time constant
        ARthresholddB = 85;       % dB SPL (input signal level) =>200 to disable
        MOCtau = 0.3;
        MOCfactor = 0.5;   %dB per dB OUTPUT
        
        numSamples = 1024; %MAX=6912
        useAid = 0;
    end
    
    %%  *********************************************************
    % Protected properties - inter/ra class communication
    %************************************************************
    properties(Access = protected)
        jobLockFid % File identifier for mutex
        
        %Nick C comment on this:
        %OK. The big-endian thing used to work because in the config file
        %'config_tr_zcpa12' there was a flag called NATURALREADORDER that was set to
        %FALSE and thus appeared to override x86 standard:
        %little-endian. Endianess has **!nothing!** to do with win vs *nix
        byteOrder = 'le';  % byte order is big endian
    end
    properties(Dependent = true)
        jobLockTxtFile
    end
    
    %%  *********************************************************
    % methods        _   _               _
    %               | | | |             | |
    % _ __ ___   ___| |_| |__   ___   __| |___
    %| '_ ` _ \ / _ \ __| '_ \ / _ \ / _` / __|
    %| | | | | |  __/ |_| | | | (_) | (_| \__ \
    %|_| |_| |_|\___|\__|_| |_|\___/ \__,_|___/
    %************************************************************
    
    methods
        %% **********************************************************
        % Constructor
        %************************************************************
        function obj = cJob(LearnOrRecWavs, jobFolder)
            if nargin < 1
                warning('job:ambiguouscorpus',...
                    'We need to know whether the (L)earning or (R)ecognition corpus is being used, assuming ''L''');
                LearnOrRecWavs = 'L';
            end
            if nargin > 1
                obj.opFolder = jobFolder;
            else
                if isunix
                    if ismac
                        obj.opFolder = '~/ASR/exps/_foo';
                    else
                        obj.opFolder = '/scratch/nrclark/exps/_foo';
                    end
                else
                    obj.opFolder = 'D:\exps\_foo';
                end
            end
            
            obj = obj.assignWavPaths(LearnOrRecWavs);
            obj = obj.initMAP;
            
            
        end % ------ OF CONSTRUCTOR
        
        %% **********************************************************
        % Set Wav Paths
        %************************************************************
        function obj = assignWavPaths(obj, LearnOrRecWavs)
            if isunix
                if ismac
                    lWAVpath = fullfile('demo_wavs', 'TrainingData-Clean');
                    rWAVpath = fullfile('demo_wavs', 'TripletTestData');
                    obj.noiseFolder = fullfile('demo_wavs', 'noises');
                else
                    lWAVpath = fullfile('demo_wavs', 'TrainingData-Clean');
                    rWAVpath = fullfile('demo_wavs', 'TripletTestData');
                    obj.noiseFolder = fullfile('demo_wavs', 'noises');
                end
            else
                lWAVpath = 'D:\AURORA digits (wav)\TrainingData-Clean';
                rWAVpath = 'D:\AURORA digits (wav)\TripletTestData';
                obj.noiseFolder = 'D:\AURORA digits (wav)\noises';
            end
            
            if strcmpi(LearnOrRecWavs, 'l')
                obj.wavFolder = lWAVpath;
            elseif strcmpi(LearnOrRecWavs, 'r')
                obj.wavFolder = rWAVpath;
            else
                error('First argument to constructor must be ''L'' or ''R''');
            end
        end
        
        %% **********************************************************
        % mutex related    _
        %                 | |
        %  _ __ ___  _   _| |_ _____  __
        % | '_ ` _ \| | | | __/ _ \ \/ /
        % | | | | | | |_| | ||  __/>  <
        % |_| |_| |_|\__,_|\__\___/_/\_\
        %************************************************************
        
        %% **********************************************************
        % lockJobList - File mutex protecting from race conditions
        %************************************************************
        function obj = lockJobList(obj)
            lockON = false;
            while (~lockON)
                if numel(dir(obj.jobLockTxtFile)) %Check to see if lock already in place
                    wTime = randi(750)+250; %3,000-10,000 ms
                    disp(['File mutex in place. Retrying in ' num2str(wTime) ' ms'])
                    pause(wTime/1000);
                else
                    obj.jobLockFid = fopen(obj.jobLockTxtFile,'w');
                    disp('Job locked');
                    pause(50/1000);
                    lockON = true;
                end
            end
            fclose(obj.jobLockFid);
        end% ------ OF LOCKJOBLIST
        
        %% **********************************************************
        % unlockJobList - unlocks for other workers
        %************************************************************
        function obj = unlockJobList(obj)
            lockON = logical(numel(dir(obj.jobLockTxtFile)));
            while(lockON)
                try
                    delete(obj.jobLockTxtFile);
                    disp('Job unlocked');
                    pause(50/1000);
                    lockON = false;
                catch %#ok<CTCH>
                    disp('Unjamming lock - retrying immediately')
                    pause(200/1000)
                end
            end
        end% ------ OF UNLOCKJOBLIST
        
        %% **********************************************************
        % storeSelf - save a copy of this object in opFolder
        %************************************************************
        function storeSelf(obj)
            doneSaving = false;
            while(~doneSaving)
                try
                    save(fullfile(obj.opFolder, 'jobObject'), 'obj');
                    doneSaving = true;
                catch  %#ok<CTCH>
                    wTime = randi(3800)+200; %200-4000 ms
                    disp(['Save collision (THIS IS VERY VERY BAD - have you not used the mutex?). Retrying in ' num2str(wTime) ' ms']);
                    pause(wTime/1000);
                end
            end
        end % ------ OF STORESELF
        
        %% **********************************************************
        % loadSelf - load a copy of this object from opFolder
        %************************************************************
        function obj = loadSelf(obj)
            doneLoading = false;
            while(~doneLoading)
                try
                    load(fullfile(obj.opFolder,'jobObject'));
                    doneLoading = true;
                catch  %#ok<CTCH>
                    wTime = randi(3800)+200; %200-4000 ms
                    disp(['Load collision (THIS IS VERY VERY BAD - have you not used the mutex?). Retrying in ' num2str(wTime) ' ms'])
                    pause(wTime/1000);
                end
            end
        end% ------ OF LOADSELF
        
        
        %%  *********************************************************
        %  _                          _                   _
        % | |                        | |                 (_)
        % | |__   ___  _   _ ___  ___| | _____  ___ _ __  _ _ __   __ _
        % | '_ \ / _ \| | | / __|/ _ \ |/ / _ \/ _ \ '_ \| | '_ \ / _` |
        % | | | | (_) | |_| \__ \  __/   <  __/  __/ |_) | | | | | (_| |
        % |_| |_|\___/ \__,_|___/\___|_|\_\___|\___| .__/|_|_| |_|\__, |
        %                                          | |             __/ |
        %                                          |_|            |___/
        %************************************************************
        
        %% **********************************************************
        % get method for dependtent var jobLockTxtFile
        %************************************************************
        function value = get.jobLockTxtFile(obj)
            %Name the MUTEX file here
            value = [fullfile(obj.opFolder, 'jobLock') '.txt'];
        end
        
        %% **********************************************************
        % checkStatus - see how much of the current job is complete
        %************************************************************
        function checkStatus(obj)
            NOWopen = sum(obj.todoStatus==0); %Starting from R2010b, Matlab supports enumerations. For now, we resort to integers for compatibility.
            NOWpend = sum(obj.todoStatus==1);
            NOWdone = sum(obj.todoStatus==2);
            
            zz = clock;
            disp([num2str(zz(3)) '-' num2str(zz(2)) '-' num2str(zz(1)) '   ' num2str(zz(4)) ':' num2str(zz(5))])
            disp(['CURRENT JOB:' obj.opFolder]);
            disp(' ')
            disp(['open - ' num2str(NOWopen) ' || pending - ' num2str(NOWpend) ' || complete - ' num2str(NOWdone)])
            
            % ---- PROGRESS BAR ----
            pcDone = 100*NOWdone/obj.numWavs;
            progBarLength = 40;
            charBars = repmat('=',1,floor(pcDone/100 * progBarLength));
            charWhiteSpace = repmat(' ',1, progBarLength - numel(charBars));
            disp(' ')
            disp([' -[' charBars charWhiteSpace ']-  ' num2str(pcDone, '%0.1f') '%'])
            disp(' ')
            disp(' ')
            % -- END PROGRESS BAR ---
        end% ------ OF CHECKSTATUS
        
        %% **********************************************************
        % initMAP - add MAP stuff to path
        %************************************************************
        function obj = initMAP(obj)
            addpath(...fullfile(obj.MAProot, 'modules'),...
                fullfile(obj.MAProot, 'utilities'),...
                fullfile(obj.MAProot, 'MAP'),...
                fullfile(obj.MAProot, 'parameterStore'));
        end % ------ OF INIT MAP
        
        %% **********************************************************
        % assign files to testing and training sets
        %************************************************************
        function obj = assignFiles(obj)
            speechWavs  = dir(fullfile(obj.wavFolder, '*.wav'));
            assert(obj.numWavs <= size(speechWavs, 1) ,...
                'not enough files available in the corpus');  % make sure we have enough wavs
            
            randomWavs = rand(1, size(speechWavs, 1));
            [~,b] = sort(randomWavs);
            trainFileIdx = b(1:obj.numWavs);
            
            obj.wavList  = speechWavs(trainFileIdx); %This is a record of all of the wavs that should be done
            
            %Starting from R2010b, Matlab should support enumerated types. For now we
            %use integers for compatibility.
            %0=open, 1=processing, 2=done
            obj.todoStatus = zeros(numel(obj.wavList), 1);
            
        end % ------ OF ASSIGN FILES
        
        %% **********************************************************
        % generate  feature
        %************************************************************
        function obj = genFeat(obj, currentWav)
            fprintf(1,'Processing: %s \n', currentWav);
            if strcmpi(obj.speechDist,'Gaussian')
                tempSpeechLev = obj.speechLevToUse + obj.speechLevStd*randn;
            elseif strcmpi(obj.speechDist,'Uniform')
                % for a uniform distribution, the standard deviation is
                % range/sqrt(12)
                % http://mathforum.org/library/drmath/view/52066.html
                tempSpeechLev = obj.speechLevToUse - obj.speechLevStd*sqrt(12)/2 + obj.speechLevStd*sqrt(12)*rand;
            elseif strcmpi(obj.speechDist,'None')
                tempSpeechLev = obj.speechLevToUse;
            end
            
            if strcmpi(obj.noiseDist,'Gaussian')
                tempNoiseLev  = speechLev - obj.meanSNR  + obj.noiseLevStd*randn;
            elseif strcmpi(obj.noiseDist,'Uniform')
                tempNoiseLev  = tempSpeechLev - obj.meanSNR - obj.noiseLevStd*sqrt(12)/2 + obj.noiseLevStd*sqrt(12)*rand;
            elseif strcmpi(obj.noiseDist,'None')
                tempNoiseLev  = obj.noiseLevToUse;
            end
            
            disp(['Current speech level = ' num2str(tempSpeechLev)]);
            disp(['Current noise level  = ' num2str(tempNoiseLev)]);
            
            obj.currentSpeechLevel = tempSpeechLev;
            obj.currentNoiseLevel = tempNoiseLev;
            [finalFeatures, ~] = processWavs(obj, currentWav); %discard the output from ANprobabilityResponse and method using ~
            opForHTK(obj, currentWav, finalFeatures);
        end % ------ OF GENFEAT
        
        %% **********************************************************
        % write o/p in HTK friendly format
        %************************************************************
        function obj = opForHTK(obj, currentWav, featureData)
            
            featureName = strrep(currentWav, '.wav','.map');
            targetFilename = fullfile(obj.opFolder, featureName);
            
            % write in a format HTK compliant for the recogniser to use
            obj.writeHTK(...
                targetFilename,...
                featureData,...
                size(featureData,2),...
                obj.frameshift*obj.sampPeriodFromMsFactor,...
                size(featureData,1)*4,...
                obj.paramKind,...
                obj.byteOrder);            
        end % ------ OF opForHTK
        
        
        %% **********************************************************
        %      _                   _                                     _
        %     (_)                 | |                                   (_)
        %  ___ _  __ _ _ __   __ _| |  _ __  _ __ ___   ___ ___  ___ ___ _ _ __   __ _
        % / __| |/ _` | '_ \ / _` | | | '_ \| '__/ _ \ / __/ _ \/ __/ __| | '_ \ / _` |
        % \__ \ | (_| | | | | (_| | | | |_) | | | (_) | (_|  __/\__ \__ \ | | | | (_| |
        % |___/_|\__, |_| |_|\__,_|_| | .__/|_|  \___/ \___\___||___/___/_|_| |_|\__, |
        %         __/ |               | |                                         __/ |
        %        |___/                |_|                                        |___/
        %************************************************************
        
        %% **********************************************************
        % getStimulus - what it says on the tin
        %************************************************************
        function [stimulus, sampleRate] = getStimulus(obj, currentWav)
            
            % getStimulus.m - NC Aug 2010
            % Modified version of Rob's script to include:
            %
            % 1)Signal and noise samples with different sample rate. The component with
            % lower sample rate is upsampled to the rate of that with the
            % higher rate.
            % 2) Clearer level setting
            % 3) Parameter to change noise intro duration
            % 4) Noise padding at end of stimulus
            %
            % ORIGINAL HEADER:
            % getStimulus.m
            %
            % Robert T. Ferry
            % 13th May 2009
            
            % Set levels
            [speech speechSampleRate] = wavread(fullfile(obj.wavFolder, currentWav ));
            speech = speech./sqrt(mean(speech.^2)); %Normalize RMS to 1
            speech = speech * 20e-6 * 10^(obj.currentSpeechLevel/20); %Convert RMS to pascals at desired level
            %disp(20*log10(sqrt(mean(speech.^2))/20e-6))
            
            [noise noiseSampleRate] = wavread(fullfile(obj.noiseFolder, obj.noiseName ));
            noise = noise./sqrt(mean(noise.^2)); %Normalize RMS to 1
            noise = noise * 20e-6*10^(obj.currentNoiseLevel/20); %Convert RMS to pascals at desired level
            %disp(20*log10(sqrt(mean(noise.^2))/20e-6))
                        
            % Do sample rate conversion if needed
            % Will always convert stimulus component with lower rate up to that with
            % higher rate.
            if speechSampleRate > noiseSampleRate
                %     disp('S>N')
                [p,q] = rat(speechSampleRate/noiseSampleRate,0.0001);
                noise = resample(noise, p, q);
                noiseSampleRate = speechSampleRate;
            elseif noiseSampleRate > speechSampleRate
                %     disp('N>S')
                [p,q] = rat(noiseSampleRate/speechSampleRate,0.0001);
                speech = resample(speech, p, q);
                speechSampleRate = noiseSampleRate; %#ok<NASGU>
            else
                %DO NOTHING BUT ASSERT
                assert(speechSampleRate == noiseSampleRate);
            end
            sampleRate = noiseSampleRate;
            dt = 1/sampleRate;
            
            % mix stimuli
            % Everything from here down (mostly) is RTF's original
            silenceStart = floor(obj.noisePreDur*sampleRate);
            silenceEnd = floor(obj.noisePostDur*sampleRate);
            
            silencePointsStart = zeros(silenceStart,1);
            silencePointsEnd = zeros(silenceEnd,1);
            
            speech = [silencePointsStart; speech; silencePointsEnd];
            
            stimLength = length(speech);
            noiseLength = length(noise);
            if obj.freezeNoise
                idx = 1;
            else
                idx = ceil(rand*(noiseLength-stimLength));
            end
            noise = noise(idx:idx+stimLength-1);
            
            stimulus = speech+noise;
            
            % add ramps to noise
            stimInNoiseTime = dt:dt:dt*length(stimulus);
            rampDuration = 0.100;
            rampTime = dt : dt : rampDuration;
            ramp = [0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(stimInNoiseTime)-length(rampTime))];
            stimulus = stimulus'.*ramp;
            ramp = fliplr(ramp);
            stimulus = stimulus.*ramp;
            stimulus = stimulus';
            %disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))
            
            % add additional silent period to start of stimulus for model to 'settle down'
            additionalSilenceLength = round(0.050*sampleRate);
            additionalSilencePoints = zeros(additionalSilenceLength,1);
            stimulus = [additionalSilencePoints; stimulus]'; %now rotated.
            %disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))
        end% ------ OF GETSTIMULUS
                
        
        %% **********************************************************
        % processWavs - do all of the MAP signal processing
        %************************************************************
        function [finalFeatures, ANprobabilityResponse] = processWavs(obj, currentWav)
            
            %**********************************************************
            % FIRST GET THE STIMULUS
            %**********************************************************
            [stimulus, sampleRate] = obj.getStimulus(currentWav);
            
            %**********************************************************
            % NOW TO LOAD IN THE HEARING AID
            %**********************************************************
            if obj.useAid
                stimulus = [stimulus; stimulus]'; %EsxAid requires stereo stim
                stimulus = EssexAid_guiEmulatorWrapper(stimulus, sampleRate, obj);
                stimulus = stimulus(1,:); %convert back to mono
            end
            
            AN_spikesOrProbability = 'probability';
            
            if obj.useSpectrogram
                lowestBF=100; 	highestBF= 4500; 	numChannels=30;
                F=round(logspace(log10(lowestBF),log10(highestBF),numChannels));
                
                nfft = 1024;
                hopSamples = 64;
                noverlap = nfft - hopSamples;
                dt = hopSamples/sampleRate;
                [~,~,~,P] = spectrogram(stimulus,nfft,noverlap,F,sampleRate);
                
                ANprobabilityResponse = 10*log10(  abs(P) /  ((20e-6)^2)  ); %now correct [(a^2)/(b^2) = (a/b)^2]
                
            else
                [ANprobabilityResponse, dt, myBFlist] = MAPwrap(stimulus, sampleRate, -1, obj.participant, AN_spikesOrProbability, obj.MAPparamChanges);
            end
            nChannels = numel(myBFlist);
            
            
            time_ANresponse = dt:dt:dt*size(ANprobabilityResponse,2);
            idx = time_ANresponse > obj.truncateDur; %RTF had this @ 0.550
            ANprobabilityResponse = ANprobabilityResponse(:, idx);
            
            
            if ~obj.useSpectrogram
                if obj.MAPopLSR && ~obj.MAPopHSR
                    ANprobabilityResponse = ANprobabilityResponse(1:nChannels, :); %use LSR
                end
                if ~obj.MAPopLSR && obj.MAPopHSR
                    ANprobabilityResponse = ANprobabilityResponse(end-nChannels+1:end, :); %or use HSR
                end
                if obj.MAPopMSR
                    assert(0,'Not working with MSR at the mo')
                end
            end
            
            % OPTIONAL PLOTTING
            YTickIdx = 1:floor(numel(myBFlist)/6):numel(myBFlist);
            YTickIdxRev = numel(myBFlist)+1-YTickIdx;
            if ~isempty(obj.probHaxes)
                axes(obj.probHaxes);  %#ok<MAXES>
                imagesc(ANprobabilityResponse); colorbar('peer', obj.probHaxes)
                set(obj.probHaxes, 'YTick', YTickIdx);
                set(obj.probHaxes, 'YTickLabel', num2str(    myBFlist(YTickIdxRev)'     ));
                ylabel('cf in Hz')
            end
            
            % OPTIONAL PLOTTING SMOOTHED
            if ~isempty(obj.probHaxesSM)
                axes(obj.probHaxesSM); %#ok<MAXES>
                anSM=flipud(obj.makeANsmooth(ANprobabilityResponse, 1/dt));
                imagesc((1:size(anSM,2))./100,1:size(ANprobabilityResponse,1),anSM);
                set(obj.probHaxesSM, 'YTick', YTickIdx);
                set(obj.probHaxesSM, 'YTickLabel', num2str(    myBFlist(YTickIdxRev)'     ));
                shading(obj.probHaxesSM, 'flat'); colorbar('peer', obj.probHaxesSM)
                ylabel('cf (Hz)')
                xlabel('Time (s)')
            end
            
            
            %**********************************************************
            % optional SACF stage
            %**********************************************************
            if obj.useSACF
                
                % A slightly ugly copying is needed
                SACFparams.lambda = obj.SACFlambda;
                SACFparams.acfTau = obj.SACFacfTau;
                SACFparams.lags = logspace(log10(obj.SACFminLag),log10(obj.SACFmaxLag),obj.SACFnBins);
                SACFparams.lags = linspace(obj.SACFminLag, obj.SACFmaxLag,obj.SACFnBins );
                
                SACFmethod.dt = dt;
                SACFmethod.nonlinCF = myBFlist;
                
                %This is slightly misleading as the ANprob is now a SACF
                [ANprobabilityResponse, ~, ~, ~] = filteredSACF(ANprobabilityResponse, SACFmethod, SACFparams);
                
                % OPTIONAL PLOTTING
                YTickIdx = 1:floor(obj.SACFnBins/6):obj.SACFnBins;
                %YTickIdxRev = obj.SACFnBins+1-YTickIdx;
                if ~isempty(obj.sacfHaxes)
                    axes(obj.sacfHaxes);  %#ok<MAXES>
                    imagesc(flipud(ANprobabilityResponse)); shading(obj.sacfHaxes, 'flat'); colorbar('peer', obj.sacfHaxes)
                    set(obj.sacfHaxes, 'YTick', YTickIdx);
                    set(obj.sacfHaxes, 'YTickLabel', num2str(    1./SACFparams.lags(YTickIdx)'    ,'%0.1f' ));
                    ylabel('Pitch in Hz')
                end
                
                % OPTIONAL PLOTTING SMOOTHED
                if ~isempty(obj.sacfHaxesSM)
                    axes(obj.sacfHaxesSM);  %#ok<MAXES>
                    imagesc(flipud(obj.makeANsmooth(ANprobabilityResponse, 1/dt))); shading(obj.sacfHaxesSM, 'flat'); colorbar('peer', obj.sacfHaxesSM)
                    set(obj.sacfHaxesSM, 'YTick', YTickIdx);
                    set(obj.sacfHaxesSM, 'YTickLabel', num2str(    1./SACFparams.lags(YTickIdx)'    ,'%0.1f' ));
                    ylabel('Pitch in Hz')
                end
                
            end
            
            
            finalFeatures = obj.makeANfeatures(  ...
                obj.makeANsmooth(ANprobabilityResponse, 1/dt), obj.numCoeff  );
            
            if obj.removeEnergyStatic
                finalFeatures = finalFeatures(2:end,:);
                % disp(size(finalFeatures))
            end
            
            if obj.doCMN
                m = mean(finalFeatures,2);
                finalFeatures = finalFeatures - repmat(m,1,size(finalFeatures,2));
            end
            
            % OPTIONAL PLOTTING
            if ~isempty(obj.featHaxes)
                pcolor(obj.featHaxes, finalFeatures); shading(obj.featHaxes, 'flat'); colorbar('peer', obj.featHaxes)
            end
            if ~isempty(obj.reconHaxes)
                reconsData = idct(finalFeatures,obj.SACFnBins);
                axes(obj.reconHaxes);  %#ok<MAXES>
                imagesc(flipud( reconsData ));
            end
            
            opForHTK(obj, currentWav, finalFeatures);
        end % ------ OF PROCESSWAVS
        
    end % ------ OF METHODS
    
    %%  *********************************************************
    %      _        _   _                       _   _               _
    %     | |      | | (_)                     | | | |             | |
    %  ___| |_ __ _| |_ _  ___   _ __ ___   ___| |_| |__   ___   __| |___
    % / __| __/ _` | __| |/ __| | '_ ` _ \ / _ \ __| '_ \ / _ \ / _` / __|
    % \__ \ || (_| | |_| | (__  | | | | | |  __/ |_| | | | (_) | (_| \__ \
    % |___/\__\__,_|\__|_|\___| |_| |_| |_|\___|\__|_| |_|\___/ \__,_|___/
    %************************************************************
    
    methods(Static)
        %% ********************************************************
        % makeANsmooth - smooth the AN response into hanning windowed chunks
        %**********************************************************
        function ANsmooth = makeANsmooth(ANresponse, sampleRate, winSize, hopSize)
            if nargin < 3
                winSize = 25; %default 25 ms window
            end
            if nargin < 4
                hopSize = 10; %default 10 ms jump between windows
            end
            
            winSizeSamples = round(winSize*sampleRate/1000);
            hopSizeSamples = round(hopSize*sampleRate/1000);
            
            % smooth
            hann = cJob.NRC_hanning(winSizeSamples);
            
            ANsmooth = [];%Cannot pre-allocate a size as it is unknown until the enframing
            for chan = 1:size(ANresponse,1)
                f = cJob.enframe(ANresponse(chan,:), hann, hopSizeSamples);
                ANsmooth(chan,:) = mean(f,2)'; %#ok<AGROW> see above comment
            end
        end% ------ OF makeANsmooth
        
        %% ********************************************************
        % makeANfeatures - dct wizardry
        %**********************************************************
        function ANfeatures = makeANfeatures(ANrate, numCoeff)
            % make feature vectors
            features = cJob.GJB_dct(ANrate);
            ANfeatures = features(1:numCoeff,:);
        end % ------ OF makeANfeatures
        
        %% ************************************************************************
        % enframe - AVOID SIGNAL PROCESSING TOOLBOX buffer function
        %**************************************************************************
        function f=enframe(x,win,inc)
            %ENFRAME split signal up into (overlapping) frames: one per row. F=(X,WIN,INC)
            %
            %	F = ENFRAME(X,LEN) splits the vector X(:) up into
            %	frames. Each frame is of length LEN and occupies
            %	one row of the output matrix. The last few frames of X
            %	will be ignored if its length is not divisible by LEN.
            %	It is an error if X is shorter than LEN.
            %
            %	F = ENFRAME(X,LEN,INC) has frames beginning at increments of INC
            %	The centre of frame I is X((I-1)*INC+(LEN+1)/2) for I=1,2,...
            %	The number of frames is fix((length(X)-LEN+INC)/INC)
            %
            %	F = ENFRAME(X,WINDOW) or ENFRAME(X,WINDOW,INC) multiplies
            %	each frame by WINDOW(:)
            
            %	   Copyright (C) Mike Brookes 1997
            %      Version: $Id: enframe.m,v 1.4 2006/06/22 19:07:50 dmb Exp $
            %
            %   VOICEBOX is a MATLAB toolbox for speech processing.
            %   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   This program is free software; you can redistribute it and/or modify
            %   it under the terms of the GNU General Public License as published by
            %   the Free Software Foundation; either version 2 of the License, or
            %   (at your option) any later version.
            %
            %   This program is distributed in the hope that it will be useful,
            %   but WITHOUT ANY WARRANTY; without even the implied warranty of
            %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
            %   GNU General Public License for more details.
            %
            %   You can obtain a copy of the GNU General Public License from
            %   http://www.gnu.org/copyleft/gpl.html or by writing to
            %   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nx=length(x(:));
            nwin=length(win);
            if (nwin == 1)
                len = win;
            else
                len = nwin;
            end
            if (nargin < 3)
                inc = len;
            end
            nf = fix((nx-len+inc)/inc);
            f=zeros(nf,len);
            indf= inc*(0:(nf-1)).';
            inds = (1:len);
            f(:) = x(indf(:,ones(1,len))+inds(ones(nf,1),:));
            if (nwin > 1)
                w = win(:)';
                f = f .* w(ones(nf,1),:);
            end
        end % ------ OF ENFRAME
        
        %% ************************************************************************
        % GJB_dct - AVOID SIGNAL PROCESSING TOOLBOX
        %**************************************************************************
        function b=GJB_dct(a,n)
            
            if nargin == 0,
                error('Not enough input arguments.');
            end
            
            if isempty(a)
                b = [];
                return
            end
            
            % If input is a vector, make it a column:
            do_trans = (size(a,1) == 1);
            if do_trans, a = a(:); end
            
            if nargin==1,
                n = size(a,1);
            end
            m = size(a,2);
            
            % Pad or truncate input if necessary
            if size(a,1)<n,
                aa = zeros(n,m);
                aa(1:size(a,1),:) = a;
            else
                aa = a(1:n,:);
            end
            
            % Compute weights to multiply DFT coefficients
            ww = (exp(-1i*(0:n-1)*pi/(2*n))/sqrt(2*n)).';
            ww(1) = ww(1) / sqrt(2);
            
            if rem(n,2)==1 || ~isreal(a), % odd case
                % Form intermediate even-symmetric matrix
                y = zeros(2*n,m);
                y(1:n,:) = aa;
                y(n+1:2*n,:) = flipud(aa);
                
                % Compute the FFT and keep the appropriate portion:
                yy = fft(y);
                yy = yy(1:n,:);
                
            else % even case
                % Re-order the elements of the columns of x
                y = [ aa(1:2:n,:); aa(n:-2:2,:) ];
                yy = fft(y);
                ww = 2*ww;  % Double the weights for even-length case
            end
            
            % Multiply FFT by weights:
            b = ww(:,ones(1,m)) .* yy;
            
            if isreal(a), b = real(b); end
            if do_trans, b = b.'; end
        end % ----- of GJB_DCT
        
        %% ************************************************************************
        % NRC_hanning - AVOID SIGNAL PROCESSING TOOLBOX
        %**************************************************************************
        function w=NRC_hanning(n)
            calc_hanning = @(m,n)0.5*(1 - cos(2*pi*(1:m)'/(n+1))); %cheeky anonymous function - I <3 these
            if ~rem(n,2)
                % Even length window
                half = n/2;
                w = calc_hanning(half,n);
                w = [w; w(end:-1:1)];
            else
                % Odd length window
                half = (n+1)/2;
                w = calc_hanning(half,n);
                w = [w; w(end-1:-1:1)];
            end
        end % ------ of NRC_HANNING
        
        %% ************************************************************************
        % writeHTK - convert data to htk format -> by anonymous Sue (2001)
        %**************************************************************************
        function retcode = writeHTK(filename, htkdata, nFrames, sampPeriod, SampSize, ParamKind, byte_order)
            % Write an HTK format file.
            %
            % Input parameters:
            %    filename		HTK data file
            %    htkdata      HTK data read: an m x n matrix with
            %                    m = no. of channels
            %                    n = no. of frames
            %  The following are from the HTK header (see HTK manual):
            %    nFrames      no. of frames (samples)
            %    sampPeriod   sample period (in 100 ns units?)
            %    SampSize     sample size
            %    ParamKind    parameter kind code
            %
            %    byteorder    'be' for big-endian (typical for Unix) (default)
            %                 'le' for little-endian (typical for MSWindows)
            %
            % Output parameters:
            %    retcode      0 if successful
            
            % Written by Sue 17/12/01
            
            retcode=-1;	% initialise in case of error
            if nargin < 6
                fprintf('Usage: %s(filename, htkdata, nFrames, sampPeriod, SampSize, ParamKind [, byte_order])', mfilename);
            end;
            
            % Default to big-endian (HTK format)
            if nargin < 7
                byte_order = 'be';
            end;
            
            fid = fopen (filename, 'w', sprintf('ieee-%s', byte_order));
            if fid < 1
                fprintf('%s: can''t open output file %s\n', mfilename, filename);
                return
            end
            
            % Write header
            fwrite (fid, nFrames, 'int32'); %nSamples in HTK
            fwrite (fid, sampPeriod, 'int32');
            fwrite (fid, SampSize, 'int16');
            fwrite (fid, ParamKind, 'int16');
            
            % Write actual data
            fwrite(fid, htkdata, 'float32');
            
            fclose(fid);
            
            retcode=0;
        end% ------ OF WRITEHTK
        
        %% ************************************************************************
        % readHTK - just incase you ever want to go backwards
        %**************************************************************************
        function [htkdata,nframes,sampPeriod,sampSize,paramKind] = readHTK(filename,byte_order)
            
            if nargin<2
                byte_order = 'be';
            end
            
            fid = fopen(filename,'r',sprintf('ieee-%s',byte_order));
            
            nframes = fread(fid,1,'int32');
            sampPeriod = fread(fid,1,'int32');
            sampSize = fread(fid,1,'int16');
            paramKind = fread(fid,1,'int16');
            
            % read the data
            
            htkdata = fread(fid,nframes*(sampSize/4),'float32');
            htkdata = reshape(htkdata,sampSize/4,nframes);
            fclose(fid);
        end % ------ OF READHTK
        
    end % ------ OF STATIC METHODS
    
end % ------ OF CLASS
