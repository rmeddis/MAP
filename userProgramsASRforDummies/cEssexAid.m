classdef cEssexAid
    %ESSEXAID_WRAPCLASS Wrapper for the EssexAid - Nick Clark July 2011
    %   This class wraps up the EssexAid algorithm function that processes
    %   each block of samples. This wrapper closely emulates the GUI used
    %   in the lab and runs stimuli through the exact same algorithm used
    %   in the lab. It even includes a helper function to generate C code
    %   from the algorithm for use in a real-time framework.

    
    %% *********************************************************
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
        sr         = 48e3;
        numSamples = 1024; %MAX=6912, LAB_USE=48
        stimulusUSER                   
        
        %------------------------------------------------------------------
        % Params for audiometric freqs 250, 500, 1000, 2000, 4000, 8000 Hz
        %------------------------------------------------------------------
        audiometry_dB= [ 0;    0;    0;    0;    0;   0];   %Pure tone threshold in dB SPL
        mainGain_dB  = [ 0;    0;    0;    0;    0;   0];   %Gain applied at audiometric frequencies
        TC_dBHL      = [40;   40;   40;   40;   40;  40];   %Compression thresholds (in dB HL from 2nd filt)
        TM_dBHL      = [10;   10;   10;   10;   10;  10];   %MOC thresholds (in dB OUTPUT from 2nd filt)
        DRNLc        = [ 0.2;  0.2;  0.2;  0.2;  0.2; 0.2]; %Compression exponent at audiometric frequencies
        
        %------------------------------------------------------------------
        % Dynamic compression properties
        %------------------------------------------------------------------
        ARtau = 60e-3;       %decay time constant
        ARthreshold_dB = 85; %dB SPL (input signal level) =>200 to disable
        MOCtau = 450e-3;     %Time constant in Seconds
        MOCfactor = 0.5;     %dB attenuation applied to the input per dB exceeding output threshold
        
        %------------------------------------------------------------------
        % Band filtering properties
        %------------------------------------------------------------------
        bwOct = 1/2; %1/1, 1/2, 1/3, 1/4, 1/5
        filterOrder = 2 %BUTTER=2, GTF=3
        useGTF = false; %If false, revert to butterworth
    end
    
    %% **********************************************************
    % Read only properties that are not dependent
    %************************************************************
    properties(SetAccess = private)          
        MOCrecord
    end
    
    %% **********************************************************
    % Constant properties 
    %************************************************************
    properties(Constant = true, Hidden = true)  
        numAudiometricFreqs = 6;
    end
    
    %% **********************************************************
    % Dependent visable properties - calculated at runtime
    %************************************************************
    properties(Dependent = true, Hidden = false)
        channelBFs %= 250 * 2.^((0:fNmax)'*params.bwOct);        
        numChannels %= numel(channelBFs);        
        aidOPnice %aid output reformatted to be exactly the same dimensions as the input stimulus                               
    end
    
    %% **********************************************************
    % Dependent invisable properties - calculated at runtime
    %************************************************************
    properties(Dependent = true, Hidden = true)
        TC_dBO_INTERP % Compression threshold in terms of 2nd filter o/p in dB SPL
        TM_dBO_INTERP % MOC threshold in terms of 2nd filter o/p in dB SPL                    
        bwOct_INTERP
        DRNLb_INTERP %=  ( 2e-5 .* 10.^(TCdBO/20)) .^ (1-DRNLc)  ; 
        DRNLc_INTERP               
        mainGain_INTERP %Interp'd and in linear units 
        
        ARthresholdPa %=  20e-6*10^(ARthreshold_dB/20);% Pa thresh for triggering AR        
        stimulusINTERNAL %input stimulus in correct format for the Aid algo
    end
    
    %% **********************************************************
    % Protected properties - The user never needs to set
    %************************************************************
    properties(Access = protected)  
        aidOP
        emlc_z
        
        %--------------------------------------------------------------
        % ENUMERATIONS USED IN THE FRAME PROCESSOR
        %--------------------------------------------------------------
        enumC_ARb  = 0;
        enumC_ARa  = 2;
        enumC_MOCb = 4;
        enumC_MOCa = 6;
        
        % enumC_BPb1 = 8;
        % enumC_BPa1 = 13;
        % enumC_BPb2 = 18;
        % enumC_BPa2 = 23;
        % enumC_BPb3 = 28;
        % enumC_BPa3 = 33;
        % enumC_BPb4 = 38;
        % enumC_BPa4 = 43;
        
        enumS_AR  = 0;
        
        % enumS_MOC1  = 1;
        % enumS_BPin_1_1 = 2;
        % enumS_BPin_2_1 = 6;
        % enumS_BPout_1_1 = 10;
        % enumS_BPout_2_1 = 14;
        %
        % enumS_MOC2 = 18;
        % enumS_BPin_1_2 = 19;
        % enumS_BPin_2_2 = 23;
        % enumS_BPout_1_2 = 27;
        % enumS_BPout_2_2 = 31;
        % ...
    end
        
    %% **********************************************************
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
        function obj = EssexAid_WrapClass(sr, stimulus)          
            
            if nargin > 0
                obj.sr = sr;
            end
            
            if nargin > 1
                obj.stimulusUSER = stimulus;
            else         
                obj.stimulusUSER = obj.pipSequence(obj.sr);
            end
        end
        
        %% **********************************************************
        % Get method for channelBFs 
        %************************************************************
        function value = get.channelBFs(obj)
            fNmax = 5/obj.bwOct;
            value = 250 * 2.^((0:fNmax)'*obj.bwOct);
        end
        
        %% **********************************************************
        % Get method for numChannels 
        %************************************************************
        function value = get.numChannels(obj)
            value = numel(obj.channelBFs);
        end
        
        %% **********************************************************
        % Get method for ARthresholdPa 
        %************************************************************
        function value = get.ARthresholdPa(obj)
            value = 20e-6*10^(obj.ARthreshold_dB/20);% Pa thresh for triggering AR
        end
        
        %% **********************************************************
        % Get method for TC_dBO_INTERP
        %************************************************************
        function value = get.TC_dBO_INTERP(obj)
            TC_dBO = obj.audiometry_dB - obj.mainGain_dB + obj.TC_dBHL;
            value  = obj.interpPars(TC_dBO, obj.numChannels);
        end
        
        %% **********************************************************
        % Get method for TM_dBO_INTERP
        %************************************************************
        function value = get.TM_dBO_INTERP(obj)
            TM_dBO = obj.audiometry_dB - obj.mainGain_dB + obj.TM_dBHL;
            value  = obj.interpPars(TM_dBO, obj.numChannels);
        end
        
        %% **********************************************************
        % Get method for bwOct_INTERP 
        %************************************************************
        function value = get.bwOct_INTERP(obj)
            value = repmat(obj.bwOct, 1, obj.numChannels);
        end
        
        %% **********************************************************
        % Get method for DRNLb_INTERP 
        %************************************************************
        function value = get.DRNLb_INTERP(obj)
            value = ( 2e-5 .* 10.^(obj.TC_dBO_INTERP/20)) .^ (1-obj.DRNLc_INTERP); 
        end
        
        %% **********************************************************
        % Get method for DRNLc_INTERP 
        %************************************************************
        function value = get.DRNLc_INTERP(obj)
            value  = obj.interpPars(obj.DRNLc, obj.numChannels);
        end
        
        %% **********************************************************
        % Get method for mainGain_INTERP 
        %************************************************************
        function value = get.mainGain_INTERP(obj)
            mainGainLin = 10.^(obj.mainGain_dB/20); %lin units
            value  = obj.interpPars(mainGainLin, obj.numChannels);
        end                
        
        %% ***********************************************************
        % Get method for stimulus
        % -----------------------
        % The hearing aid expects a stereo signal, as the MOC control is
        % linked for left and right channels. It would be more efficient to
        % use a mono version of the aid for simulation in Matlab. However,
        % I always want to use the exact same code for the hardware in the
        % lab and current simulations. This code will make a mono signal
        % stereo if needs be and/or rotate to 2xN array.
        %*************************************************************
        function value = get.stimulusINTERNAL(obj)            
            [nRows, nCols] = size(obj.stimulusUSER);
            
            % Assume that the stimulus duration is greater than 2 samples.
            % Therefore the number of channels is the min dim.            
            [nChans, I] = min([nRows nCols]);                                                
            
            if nChans == 2
                if I == 2
                    value = obj.stimulusUSER;
                else
                    value = obj.stimulusUSER';
                end
            elseif nChans == 1 %Just to be explicit
                if I == 2
                    value = [obj.stimulusUSER obj.stimulusUSER];   
                else
                    value = [obj.stimulusUSER; obj.stimulusUSER]';
                end
            end
        end
        
        %% ***********************************************************
        % Get method for aid output
        % -----------------------
        % This get method is linked to the above internal stimulus method
        % and allows the user to extract the hearing aid output in exactly
        % the same shape and size as the original input stimulus. This is
        % very useful for the speech recognition work and presumably
        % for multithreshold also.
        %*************************************************************
        function value = get.aidOPnice(obj)
            if ~isempty(obj.aidOP)
                [nRows, nCols] = size(obj.stimulusUSER);
                
                % Assume that the stimulus duration is greater than 2 samples.
                % Therefore the number of channels is the min dim.
                [nChans, I] = min([nRows nCols]);
                
                %** The aid output will ALWAYS be a 2xN array **
                %The fist job is to remove trailing zeros that may have been
                %introduced by the framing process
                aidOPtruncated = obj.aidOP(:, 1:max([nRows nCols]));
                
                %The next task is to arrange the op like the ip
                if nChans == 2
                    if I == 1
                        value = aidOPtruncated;
                    else
                        value = aidOPtruncated';
                    end
                elseif nChans == 1 %Just to be explicit
                    if I == 1
                        value = aidOPtruncated(1,:);
                    else
                        value = aidOPtruncated(1,:)';
                    end
                end
            else % ---- of if isempty statement
                value = [];
            end
        end
        
        %% ***********************************************************
        % *** Set methods ***
        % -----------------------
        % This is a bunch of unexciting error hunting functions. They also
        % flush the aid output if any parameters change. Therefore,
        % processStim will have to be called explicity by the user once
        % again. 
        %*************************************************************
        function obj = set.stimulusUSER(obj,value)
            [nRows, nCols] = size(value);
            
            % Assume that the stimulus duration is greater than 2 samples.
            % Therefore the number of channels is the min dim.            
            nChans = min([nRows nCols]);            
            assert(nChans<3 && nChans, 'Number of stimulus channels must be 1 or 2')
            
            obj = obj.flushAidData; %flush any previous hearing aid data if the input stimulus changes
            obj.stimulusUSER = value;
        end                     
        function obj = set.sr(obj,value)
            assert(value>=20e3 && value<=192e3, 'sr must be between 20 and 192 kHz')            
            obj = obj.flushAidData;
            obj.sr = value;
        end
        function obj = set.numSamples(obj,value)
            assert(value>=48 && value<=6912, 'must be between 48 and 6912 samples')            
            obj = obj.flushAidData;
            obj.numSamples = value;
        end
        function obj = set.audiometry_dB(obj,value)
            [nRows,nCols] = size(value);
            assert(nRows==obj.numAudiometricFreqs && nCols==1, 'must be 6x1 column vector') %#ok<MCSUP>
            obj = obj.flushAidData;
            obj.audiometry_dB = value;
        end
        function obj = set.mainGain_dB(obj,value)
            [nRows,nCols] = size(value);
            assert(nRows==obj.numAudiometricFreqs && nCols==1, 'must be 6x1 column vector') %#ok<MCSUP>
            obj = obj.flushAidData;
            obj.mainGain_dB = value;
        end
        function obj = set.TC_dBHL(obj,value)
            [nRows,nCols] = size(value);
            assert(nRows==obj.numAudiometricFreqs && nCols==1, 'must be 6x1 column vector') %#ok<MCSUP>
            obj = obj.flushAidData;
            obj.TC_dBHL = value;
        end
        function obj = set.TM_dBHL(obj,value)
            [nRows,nCols] = size(value);
            assert(nRows==obj.numAudiometricFreqs && nCols==1, 'must be 6x1 column vector') %#ok<MCSUP>
            obj = obj.flushAidData;
            obj.TM_dBHL = value;
        end
        function obj = set.DRNLc(obj,value)
            [nRows,nCols] = size(value);
            assert(nRows==obj.numAudiometricFreqs && nCols==1, 'must be 6x1 column vector') %#ok<MCSUP>
            assert(all(value)>=0 && all(value)<=1, 'all DRNLc values must be between 0 and 1')
            obj = obj.flushAidData;
            obj.DRNLc = value;
        end
        function obj = set.ARtau(obj,value)
            assert(value>=1e-3 && value<=1, 'must be between 1e-3 and 1s')            
            obj = obj.flushAidData;
            obj.ARtau = value;
        end
        function obj = set.ARthreshold_dB(obj,value)
            assert(value>0, 'set AR to a high value to disable it')            
            obj = obj.flushAidData;
            obj.ARthreshold_dB = value;
        end
        function obj = set.MOCtau(obj,value)
            assert(value>=1e-3 && value<=2, 'must be between 1e-3 and 2s')            
            obj = obj.flushAidData;
            obj.MOCtau = value;
        end
        function obj = set.MOCfactor(obj,value)
            assert(value>=0 && value<=1, 'must be between 0 and 1')
            obj = obj.flushAidData;
            obj.MOCfactor = value;
        end
        function obj = set.bwOct(obj,value)
            assert(value==1/1 || value==1/2 || value==1/3 || value==1/4 || value==1/5, 'must be one of 1./(1:5)')
            obj = obj.flushAidData;
            obj.bwOct = value;
        end
        function obj = set.filterOrder(obj,value)
            assert(value>0 && value<5, 'must be one of 1:4')
            obj = obj.flushAidData;
            obj.filterOrder = value;
        end
        function obj = set.useGTF(obj,value)
            obj = obj.flushAidData;
            obj.useGTF = value;
        end                    
        
        %% **********************************************************
        % flushAidData        
        % This second function is a workaround allowing a set method to
        % change another property value.
        %************************************************************
        function obj = flushAidData(obj) 
            obj.aidOP = [];
            obj.MOCrecord = [];
        end 
        
        
        %% **********************************************************
        % OVERLOADED plot method 
        %************************************************************
        function plot(obj)
            clf
            sig2dBSPL = @(sig)20*log10(abs(sig/20e-6)+(1/(2^32)));
            dt = 1/obj.sr;
            tAxis = dt:dt:dt*size(obj.stimulusINTERNAL,1);
            
            subplot(2,1,1)
            plot(tAxis(1:length(obj.stimulusUSER)), sig2dBSPL(obj.stimulusUSER), 'k')
            if ~isempty(obj.aidOPnice)
                hold on
                plot(tAxis(1:length(obj.stimulusUSER)), sig2dBSPL(obj.aidOPnice), 'r')
            end                                   
            ylim([0 100])
            xlim([0 tAxis(length(obj.stimulusUSER))])
            title('Level response')
            xlabel('Time in seconds')
            ylabel('Level in dB SPL')
            
            subplot(2,1,2) 
            if ~isempty(obj.MOCrecord)
                imagesc(tAxis, 1:obj.numChannels, flipud(-20*log10(obj.MOCrecord)))
                colorbar                               
            end         
            title('MOC attenuation')
            xlabel('Time in seconds')
            ylabel('Band frequency in Hz')
            numSpacers = 1 + (obj.numChannels-numel(obj.DRNLc)) / (numel(obj.DRNLc)-1);
            set(gca, 'YTick', 1:numSpacers:obj.numChannels);
            set(gca, 'YTickLabel', num2str(flipud([250; 500; 1000; 2000; 4000; 8000])));
        end% ------ OVERLOADED plot method
        
        %% **********************************************************
        % OVERLOADED soundsc method 
        %************************************************************
        function soundsc(obj)
            soundsc(obj.aidOPnice, obj.sr)
        end                 
        
        %% **********************************************************
        % processStim
        %************************************************************
        function obj = processStim(obj)
            %--------------------------------------------------------------
            % EMULATION OF THE GUI PARAMETER CONVERSIONS
            %--------------------------------------------------------------
            biggestNumSamples = obj.numSamples; 
            
            filterStatesL = (zeros(3000,1));
            filterStatesR = filterStatesL;
            filterCoeffs = (zeros(5000,1));
            
            %filter coefficients
            ARcutOff=1/(2*pi*obj.ARtau);
            [b,a] = butter(1,ARcutOff/(obj.sr/2));
            filterCoeffs(obj.enumC_ARb+1:obj.enumC_ARb+2) = b;
            filterCoeffs(obj.enumC_ARa+1:obj.enumC_ARa+2) = a;
            
            MOCcutOff=1/(2*pi*obj.MOCtau);
            [bMOC,aMOC] = butter(1,MOCcutOff/(obj.sr/2));
            filterCoeffs(obj.enumC_MOCb+1:obj.enumC_MOCb+2) = bMOC;
            filterCoeffs(obj.enumC_MOCa+1:obj.enumC_MOCa+2) = aMOC;
            
            
            for filterCount = 1:obj.numChannels
                %-----------------------------------
                % nonlinear path - filter bws
                %-----------------------------------
                lowerCutOff=obj.channelBFs(filterCount)*2^(-obj.bwOct_INTERP(filterCount)/2);
                upperCutOff=obj.channelBFs(filterCount)*2^( obj.bwOct_INTERP(filterCount)/2);
                
                if obj.useGTF
                    bwHz = upperCutOff - lowerCutOff;
                    [b_DRNL,a_DRNL] = obj.gammatone(bwHz, obj.channelBFs(filterCount), 1/obj.sr);
                    filterCoeffs(10*(filterCount-1)+9 :10*(filterCount-1)+10) = b_DRNL;
                    filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+16) = a_DRNL;
                else                                             
                    [b_DRNL,a_DRNL] = butter(2,[lowerCutOff upperCutOff]/(obj.sr/2));
                    filterCoeffs(10*(filterCount-1)+9 :10*(filterCount-1)+13) = b_DRNL;
                    filterCoeffs(10*(filterCount-1)+14:10*(filterCount-1)+18) = a_DRNL;
                end
            end
            
            %--------------------------------------------------------------
            % EMULATION OF THE IO CALLBACK THREAD
            %--------------------------------------------------------------
            frameBufferL = buffer(obj.stimulusINTERNAL(:,1), obj.numSamples);
            frameBufferR = buffer(obj.stimulusINTERNAL(:,2), obj.numSamples);
            nFrames = size(frameBufferL,2);
            
            pad = zeros(1,biggestNumSamples-obj.numSamples);
            ARampL=ones(1,biggestNumSamples);
            ARampR = ARampL;
            MOCcontrol = ones(obj.numChannels, biggestNumSamples);
            
            peakIPL = zeros(5,1);
            peakOPL = peakIPL;
            rmsIPL  = peakIPL;
            rmsOPL  = peakIPL;
            
            peakIPR = peakIPL;
            peakOPR = peakIPL;
            rmsIPR  = peakIPL;
            rmsOPR  = peakIPL;
            
            MOCend = zeros(obj.numChannels,1);
            
            op = [];
            moc= [];
            for nn = 1:nFrames
                frameBufferPadL = [frameBufferL(:,nn)' pad];
                frameBufferPadR = [frameBufferR(:,nn)' pad];
                
                [ outBufferL, outBufferR, filterStatesL, filterStatesR,  ARampL, ARampR, MOCend, peakIPL, peakOPL, rmsIPL, rmsOPL, peakIPR, peakOPR, rmsIPR, rmsOPR, MOCcontrol ] =...
                    EssexAidProcessVFrameSwitchable( ...
                    frameBufferPadL,...
                    frameBufferPadR,...
                    filterStatesL,...
                    filterStatesR,...
                    filterCoeffs,...
                    obj.numChannels,...
                    obj.numSamples,...
                    ARampL,...
                    ARampR,...
                    obj.ARthresholdPa,...
                    obj.filterOrder,...
                    obj.DRNLb_INTERP,...
                    obj.DRNLc_INTERP,...
                    obj.TM_dBO_INTERP,...
                    obj.MOCfactor,...
                    peakIPL,...
                    peakOPL,...
                    rmsIPL,...
                    rmsOPL,...
                    peakIPR,...
                    peakOPR,...
                    rmsIPR,...
                    rmsOPR,...
                    MOCend,...
                    MOCcontrol,...
                    obj.mainGain_INTERP,...
                    obj.useGTF);
                                
                
                outBuffer = ( [outBufferL(:, 1:obj.numSamples); outBufferR(:, 1:obj.numSamples)] );                
                op = [op outBuffer]; %#ok<AGROW>   
                moc= [moc MOCcontrol]; %#ok<AGROW>
                
            end %End of frame processing emulation loop
            obj.aidOP = op;
            obj.MOCrecord=moc;
                       
            
        end %End of process stim method          
        
    end %End of methods block
    
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
        % pipOut - sequence of tone pips at various levels
        %**********************************************************        
        function pipOut = pipSequence(sampleRate, freq, dBlevs, pulseDur, silDur)
            if nargin < 5
                silDur = 0.3;
            end
            if nargin < 4
                pulseDur = 0.1;
            end
            if nargin < 3
                dBlevs = 20:20:100;
            end
            if nargin < 2
                freq = 500;
            end
            if nargin < 1
                sampleRate = 48e3;
            end
            
            dt = 1/sampleRate;
            tAxis = dt:dt:pulseDur;
            sPulse = sin(2*pi*freq*tAxis);
            sPulse = sPulse./sqrt(mean(sPulse.^2));
            rms2dBspl = @(dBspl)20e-6*10^(dBspl/20); %sneaky short-hand function by (ab)using function handles
            zPad = zeros(1,ceil(sampleRate*silDur));
            
            pipOut = [];
            for nn = 1:numel(dBlevs)                
                pipOut = [ pipOut sPulse*rms2dBspl(dBlevs(nn))  zPad]; %#ok<AGROW>
            end

        end% ------ OF pipSequence
        
        %% ********************************************************
        % interpPars - Linear interpolation of given parameter to mimic GUI
        % fitting functionality.
        %**********************************************************        
        function fullArray = interpPars(shortArray, numBands)
            nGUIbands = numel(shortArray);
            if numBands == nGUIbands
                fullArray = shortArray;
            else
                numSpacers = (numBands-nGUIbands) / (nGUIbands-1);
                fullArray = shortArray(1);
                for nn = 2:nGUIbands
                    fullArray = [fullArray,...
                        repmat(mean([shortArray(nn) shortArray(nn-1)]),1,numSpacers),...
                        shortArray(nn)]; %#ok<AGROW>
                end                    
            end                        
        end% ----- OF interpPars
        
        %% ********************************************************
        % gammatone - get filter coefficients
        %********************************************************** 
        function [b,a] = gammatone(bw, cf, dt)
            phi = 2 * pi * bw * dt;
            theta = 2 * pi * cf * dt;
            cos_theta = cos(theta);
            sin_theta = sin(theta);
            alpha = -exp(-phi) * cos_theta;
            b0 = 1.0;
            b1 = 2 * alpha;
            b2 = exp(-2 * phi);
            z1 = (1 + alpha * cos_theta) - (alpha * sin_theta) * 1i;
            z2 = (1 + b1 * cos_theta) - (b1 * sin_theta) * 1i;
            z3 = (b2 * cos(2 * theta)) - (b2 * sin(2 * theta)) * 1i;
            tf = (z2 + z3) / z1;
            a0 = abs(tf);
            a1 = alpha * a0;
            
            a = [b0, b1, b2];
            b = [a0, a1];
        end% ------ OF gammatone
    end% ------ OF static methods
    
end %End of classdef

