function [audio, msg, time]=stimulusCreate(globalStimParams,...
    stimComponents, doPlot)
% updated June 2007
% the only authorised version of stimulus create is the version to be found
% in MAP1_6.  Other versions are copies!!
%
% for a simple tone you need
%
% % Mandatory structure fields
%  globalStimParams.FS=100000;
%  globalStimParams.overallDuration=.1;  % s
% doPlot=1;
%
% stim.type='tone';
% stim.phases='sin';
% stim.toneDuration=.05;;
% stim.frequencies=500;
% stim.amplitudesdB=50;
% stim.beginSilence=.01;
% stim.endSilence=-1;
% stim.rampOnDur=.002;
% stim.rampOffDur=-1;
% [audio, msg]=stimulusCreate(globalStimParams, stim, doPlot);
%
% % or for multi-component stimuli
% % Mandatory structure fields
%  globalStimParams.FS=100000;
%  globalStimParams.overallDuration=.1;  % s
% ear=1;
% componentNo=1;
%
% stimComponents(ear, componentNo).type='tone';
% stimComponents(ear, componentNo).phases='sin';
% stimComponents(ear, componentNo).toneDuration=.05;;
% stimComponents(ear, componentNo).frequencies=500;
% stimComponents(ear, componentNo).amplitudesdB=50;
% stimComponents(ear, componentNo).beginSilence=.01;
% stimComponents(ear, componentNo).endSilence=-1;
% stimComponents(ear, componentNo).rampOnDur=.002;
% stimComponents(ear, componentNo).rampOffDur=-1;
%
% % All components are forced to have the same overall duration and sample rate
%
%
% [audio, msg]=stimulusCreate(globalStimParams, stimComponents);
%
%  Optional fields
%  .ears overides ear setting by component
%  globalStimParams.ears='monoticL'; % 'monoticL', 'monoticR', 'diotic', 'dichotic'
%
%  globalStimParams.filter = [leftfl leftfu rightfl right fu]
%    filter is applied separately to left and right combined sounds
%
%  correction factor is applied to all signals to compensate for differences in output devices.
% audioOutCorrection is a scalar
%  globalStimParams.audioOutCorrection=2;
%
%  globalStimParams.FFT= 1; % {0, 1}
%  globalStimParams.doPlot=1; % {0, 1}
%  globalStimParams.doPlay=1; % {0, 1}
%
%  stimComponents(ear, componentNo).filter=[100 10000 2] % [lower, upper, order] applies to an
% individual component
%
% Output arguments:
%  audio is a stereo signal, a 2-column vector
%
%
% stim.type='noise';    % {'IRN', 'irn', 'noise', 'pinkNoise'}
% stim.toneDuration=.05;;
% stim.amplitudesdB=50;
% stim.beginSilence=.01;
% stim.endSilence=-1;
% stim.rampOnDur=.002;
% stim.rampOffDur=-1;
% [audio, msg]=stimulusCreate(globalStimParams, stim);
%
% % for IRN only
% stim.type='IRN';
% stim.toneDuration=.05;;
% stim.amplitudesdB=50;
% stim.beginSilence=.01;
% stim.endSilence=-1;
% stim.rampOnDur=.002;
% stim.rampOffDur=-1;
% stim.niterations = 8;   %0 for white noise
% stim.delay = 1/150;
% stim.irnGain = 1;
% [audio, msg]=stimulusCreate(globalStimParams, stim);stimComponents.clickWidth;
%
% stim.type='clicks';    % uses frequencies for duty cycle
% stim.toneDuration=.05;;
% stim.frequencies=500;
% stim.amplitudesdB=50;
% stim.beginSilence=.01;
% stim.endSilence=-1;
% stim.rampOnDur=.002;
% stim.rampOffDur=-1;
% [audio, msg]=stimulusCreate(globalStimParams, stim);
% stimComponents.clickWidth=.0001;    % default= dt
% stimComponents.clickHeight= 20e-6;  % default= 28e-6 * 10^(stimComponents.amplitudesdB/20);
% stim.clickTimes=[.4 .6];            % times in cylce when clicks occur, default = 1
%
%

msg=''; %error messages; no message is a good message

% plotting can be set as a separate argument or as a globalstimParams
% variable. this is for backwards compatibility only
if nargin>2,
    globalStimParams.doPlot=doPlot;
end

% stimComponents(1,1).endSilence=-1;  % end silence is always computed

% perform checks and set defaults
[globalStimParams, stimComponents]=checkDescriptors(globalStimParams, stimComponents);


% create empty stereo audio of appropriate length
audio=zeros(globalStimParams.nSignalPoints, 2);
% signal=zeros(globalStimParams.nSignalPoints, 1);
dt=globalStimParams.dt;

[Nears nComponentSounds]=size(stimComponents);
for ear=1:Nears % left/ right
    % combinedSignal is the sum of all sounds in one ear
    % it is a column vector
    combinedSignal=zeros(globalStimParams.nSignalPoints,1);

    % find valid components
    % if a component is empty, it is not a validComponent and is ignored
    validStimComponents=[];
    for i=1:nComponentSounds
        if ~isempty(stimComponents(ear,i).type)
            validStimComponents=[validStimComponents i];
        end
    end

    for componentNo=validStimComponents
        % compute individual components before adding
        stim=stimComponents(ear,componentNo);
        switch stim.type
            case 'tone'
                stimulus=maketone(globalStimParams, stim);

            case 'fmTone'
                stimulus=makeFMtone(globalStimParams, stim);

            case 'OHIO'
                stim.beginSilence=0;
                stimulus=makeOHIOtones(globalStimParams, stim);

            case 'transposedStimulus'
                stim.beginSilence=0;    % necessary because of recursion
                stimulus=makeTransposedStimulus(globalStimParams, stim);

            case { 'noise', 'pinkNoise'}
                stimulus=makeNoise(globalStimParams, stim);

            case { 'whiteNoise'}
                stimulus=makeWhiteNoise(globalStimParams, stim);

            case {'IRN', 'irn'}
                stimulus=makeIRN(globalStimParams, stim);

            case {'RPN'}
                stimulus=makeRPN(globalStimParams, stim);

            case 'clicks'
                stimulus=makeClicks(globalStimParams, stim);

            case 'PressnitzerClicks'
                % kxx clicks
                % k is 1/clickRepeatFrequency
                stimulus=makePressnitzerClicks(globalStimParams, stimComponents);

            case 'PressnitzerABxClicks'
                % kxx clicks
                % k is 1/clickRepeatFrequency
                stimulus=makePressnitzerABxClicks(globalStimParams, stimComponents);

            case 'ABxClicks'
                % A=rand*k,  B=k-A, x=rand*k.
                stimulus=makeABxClicks(globalStimParams, stimComponents);

            case 'YostClicks'
                % kxx clicks
                % k is 1/clickRepeatFrequency
                stimulus=makeYostClicks(globalStimParams, stimComponents);

            case 'kxxClicks'
                % kxx clicks
                % k is 1/clickRepeatFrequency
                stimulus=makeKxxClicks(globalStimParams, stimComponents);


%             case 'babble'
%                 % NB random start in a long file
%                 [stimulus sampleRate]= wavread('babble');
%                 nPoints=round(sampleRate*...
%                     stimComponents(ear,componentNo).toneDuration);
%                 start=round(rand(1,1)*(length(stimulus)-nPoints));
%                 stimulus=stimulus(start:start+nPoints-1);
%                 rms= 20*log10((mean(stimulus.^2).^0.5)/20e-6);
%                 dBSPLrms=stimComponents(ear,componentNo).amplitudesdB;
%                 gain=10.^((dBSPLrms-rms)/20);
%                 stimulus=stimulus'*gain;

            case 'speech'
                [stimulus sampleRate]= wavread('speech');
                stimulus=stimulus';
                nPoints=sampleRate*stimComponents(ear,componentNo).toneDuration;
                if nPoints > length(stimulus)
                    initialSilence=zeros(1,nPoints-length(stimulus));
                else
                    initialSilence=[];
                    start=round(rand(1,1)*(length(stimulus)-nPoints));
                    stimulus=stimulus(start:start+nPoints-1);
                end
                rms= 20*log10((mean(stimulus.^2).^0.5)/20e-6);
                dBSPLrms=stimComponents(ear,componentNo).amplitudesdB;
                gain=10.^((dBSPLrms-rms)/20);
                stimulus=stimulus*gain;
                stimulus=[stimulus initialSilence ];


            case 'file'
                % component already read from file and stored in stimulus. Insert it here
                % additional code for establishing signal rms level
                % NB signal is always mono at this stage

                stimulus=stim.stimulus;
                dBSPL=stim.amplitudesdB;

                nPoints=round(stim.toneDuration/dt);
                [r c]=size(stimulus);
                if r>c, stimulus=stimulus'; end    % secure horizontal vector
                stimulus=stimulus(1,1:nPoints); % only mono taken from file

                try
                    % dB rms
                    rms= 20*log10((mean(stimulus.^2).^0.5)/20e-6);
                    % special request to set fixed rms for stimulus
                    dBSPLrms=stimComponents(ear,componentNo).amplitudesdBrms;
                    if ~(dBSPLrms==-1)
                        gain=10.^((dBSPLrms-rms)/20);
                        stimulus=stimulus*gain;
                    end
                catch
                    % If no request for rms level is made
                    % set dB as peak amp
                    [stimulus gain]= normalize(stimulus);
                    dBSPL=stimComponents(ear,componentNo).amplitudesdB;
                    if ~(dBSPL==-1)
                        amplitude=28e-6*10.^(dBSPL/20);
                        stimulus=stimulus*amplitude;
                    end
                end

            case 'none'
                % no stimulus
                stimulus=zeros(1,round(stim.toneDuration/dt));

            case 'digitStrings'
                % select a digit string at random anduse as target
                files=dir(['..' filesep '..' filesep 'multithresholdResources\digitStrings']);
                files=files(3:end);
                nFiles=length(files);
                fileNo=ceil(nFiles*rand);
                fileData=files(fileNo);
                fileName=['..\..\multithresholdResources\digitStrings\' fileData.name];
                [stimulus sampleRate]=wavread(fileName);
                stimulus=stimulus(:,1)';  % make into a row vector
                % estimate the extend of endsilence padding
                nPoints=sampleRate*...
                    stimComponents(ear,componentNo).toneDuration;
                if nPoints > length(stimulus)
                    endSilence=zeros(1,nPoints-length(stimulus));
                else
                    % or truncate the file
                    endSilence=[];
                    stimulus=stimulus(1:nPoints);
                end
                % compute rms before adding silence
                rms= 20*log10((mean(stimulus.^2).^0.5)/20e-6);
                dBSPLrms=stimComponents(ear,componentNo).amplitudesdB;
                gain=10.^((dBSPLrms-rms)/20);
                stimulus=stimulus*gain;
                stimulus=[stimulus endSilence];
                global stimulusParameters
                stimulusParameters.digitString=fileName(end-7:end-5);

            otherwise
                switch stim.type(end-5:end)
                    % any file name with 'Babble' at the end is a
                    % multiThreshold file
                    case 'Babble'
                        % one of the many babbles is selected.
                        % NB random start in a long file
                        %  stim.type should contain the name of the babble file
                        fileName=['..' filesep '..' filesep ...
                            'multithresholdResources' filesep ...
                            'backgrounds and maskers'...
                             filesep stim.type];
                        
                        [stimulus sampleRate]= wavread(fileName);
                        if ~isequal(sampleRate, globalStimParams.FS)
                            % NB the file will read but will disagree with
                            % tone stimuli or other files read
                            msg= ['error: file sample rate disagrees ' ...
                                'with sample rate requested in paradigm'...
                                ' file (' ...
                            num2str(globalStimParams.FS) ').'];
                            error(msg);                            
                        end
                        nPoints=round(sampleRate*...
                            stimComponents(ear,componentNo).toneDuration);
                        start=round(rand(1,1)*(length(stimulus)-nPoints));
                        stimulus=stimulus(start:start+nPoints-1);
                        rms= 20*log10((mean(stimulus.^2).^0.5)/20e-6);
                        dBSPLrms=stimComponents(ear,componentNo).amplitudesdB;
                        gain=10.^((dBSPLrms-rms)/20);
                        stimulus=stimulus'*gain;
                        
                    otherwise
                        % stim.type may be the name of a file to be read
                        % play from beginning for stimulus duration
                        try
                            [stimulus sampleRate]= wavread([stim.type '.wav']);
                        catch
                            error(['stimulusCreate: unrecognised stimulus type -> '...
                                stim.type])
                        end
                        if ~isequal(sampleRate, globalStimParams.FS)
                            % NB the file will read but will disagree with
                            % tone stimuli or other files read
                            msg= ['error: file sample rate disagrees ' ...
                                'with sample rate requested in paradigm'...
                                ' file (' ...
                            num2str(globalStimParams.FS) ').'];
                            error(msg);
                        end
                        stimulus=stimulus';  % make into a row vector
                        % estimate the extend of endsilence padding
                        nPoints=sampleRate*...
                            stimComponents(ear,componentNo).toneDuration;
                        if nPoints > length(stimulus)
                            endSilence=zeros(1,nPoints-length(stimulus));
                        else
                            % or truncate the file
                            endSilence=[];
                            stimulus=stimulus(1:nPoints);
                        end
                        % compute rms before adding silence
                        rms= 20*log10((mean(stimulus.^2).^0.5)/20e-6);
                        dBSPLrms=stimComponents(ear,componentNo).amplitudesdB;
                        gain=10.^((dBSPLrms-rms)/20);
                        stimulus=stimulus*gain;
                        stimulus=[stimulus endSilence];
                end
        end

        % NB stimulus is a row vector now!
        % audio and combinedSignal were column vectors
        % signal will be a row vector

        % filter stimulus
        try
            % if filter field is present, [lower upper order]
            if stim.filter(1)>0	% 0 means don't filter
                stimulus=Butterworth (stimulus, dt, stim.filter(1), ...
                    stim.filter(2), stim.filter(3));

            end
        catch
        end


        % apply amplitude modulation
        if isfield(stim,'AMfrequency') & isfield(stim,'AMdepth')
            if stim.AMfrequency>0 & stim.AMdepth>0
                time=dt:dt:dt*length(stimulus);
                modulator=sin(2*pi*stim.AMfrequency*time);
                modulator=modulator*stim.AMdepth/100 + 1; % depth is percent
                stimulus=stimulus.*modulator/2;
            end
        end

        % Basic stimulus is now created.
        % Add trappings, ramps, silences to main stimulus
        rampOnTime=0;		 %ramp starts at the beginning of the stimulus
        rampOffTime=stim.toneDuration-stim.rampOffDur;
        if stim.rampOnDur>0.0001
            stimulus=applyRampOn(stimulus, stim.rampOnDur, rampOnTime, 1/dt);
            stimulus=applyRampOff(stimulus, stim.rampOffDur, rampOffTime, 1/dt);
        end
        if stim.rampOnDur<-0.0001   % apply Gaussian ramp
            stimulus=applyGaussianRamps(stimulus, -stim.rampOnDur, 1/dt);
        end

        % begin silence
        % start with a signal of the right length consisting of zeros
        signal=zeros(1, globalStimParams.nSignalPoints);
        % identify start of stimulus
        insertStimulusAt=round(stim.beginSilence/dt)+1;
        % add stimulus
        endOfStimulus=insertStimulusAt+length(stimulus)-1;
        if endOfStimulus<=globalStimParams.nSignalPoints
            signal(1, insertStimulusAt: endOfStimulus)=stimulus;
        else
            error('stimulusCreate: tone too long to fit into the overall duration')
        end

        % time signature
        time=dt:dt:dt*length(signal);
        % figure(22), plot(signal), title([num2str(ear) ' - ' num2str(componentNo)]),pause (1)

        try
            % create a column vector and trap if no signal has been created
            signal=signal';
            % also traps if signals are not the same length
            combinedSignal=combinedSignal+signal;
            %             figure(21), plot(combinedSignal), title([num2str(ear) ' - ' num2str(componentNo)]),pause (1)
        catch
            % dump everything because there is a problem
            globalStimParams
            [ear  componentNo]
            stim
            size(combinedSignal)
            size(signal)
            [   size(initialSilence)            size(signal)            size(endSilence)]
            [ ear  componentNo...
                round(globalStimParams.overallDuration*globalStimParams.FS)...
                round(stim.beginSilence*globalStimParams.FS)...
                round(stim.toneDuration*globalStimParams.FS) ...
                round(stim.endSilence*globalStimParams.FS)...
                (                       round(globalStimParams.overallDuration*globalStimParams.FS)...
                -round(stim.beginSilence*globalStimParams.FS)...
                -round(stim.toneDuration*globalStimParams.FS) ...
                -round(stim.endSilence*globalStimParams.FS))...
                ]
            error(' trouble in stimulusCreate: signals are the wrong length ')
        end

    end % component no

    audio(:,ear)= combinedSignal;

    % FFT
    try
        if globalStimParams.FFT
            showFFT (audio, dt)
        end
    catch
    end


end % ear

switch globalStimParams.ears
    % normally the signals are created in appropriate ears but .ears can
    % overide this to produce a mono signal.
    case 'monoticL';
        % combine left and right ears to make a mono signal in the left ear
        audio(:,1)=audio(:,1)+audio(:,2);
        audio(:,2)=zeros(globalStimParams.nSignalPoints, 1);

    case 'monoticR';
        % combine left and right ears to make a mono signal in the right ear
        audio(:,2)=audio(:,1)+audio(:,2);
        audio(:,1)=zeros(globalStimParams.nSignalPoints, 1);

    case 'diotic';
        % combine left and right ears to make a mono signal in both ears
        bothEarsCombined=audio(:,1)+audio(:,2);
        audio(:,2)=bothEarsCombined;
        audio(:,1)=bothEarsCombined;

    otherwise
        % Any other denomination produces no effect here
end

% Plotting as required
if globalStimParams.doPlot
    figure(9), clf
    plot(time,audio(:,1)'), hold on,
    if Nears>1
        offSet=(max(audio(:,1))+max(audio(:,2)))/10;
        offSet=2*offSet+max(audio(:,1))+max(audio(:,2));
        plot(time,audio(:,2)'+offSet,'r'), hold off
    end
    ylabel('left                             right')
end

% Calibration
% all signals are divided by this correction factor
% peakAmp=globalStimParams.audioOutCorrection; % microPascals = 100 dB SPL

audio=audio/globalStimParams.audioOutCorrection;

if globalStimParams.doPlay
        if ispc
            wavplay(audio,globalStimParams.FS)
        else
            sound(audio,globalStimParams.FS)
        end

    wavplay(audio,globalStimParams.FS)
end
% all Done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%---------------------------------------------------- maketone
function tone=maketone(globalStimParams, stimComponents)
% maketone generates a stimComponents tone
% tone=maketone(dt, frequencies, toneDuration, dBSPL, phases)
% tone is returned in Pascals
% frequencies is a list of frequencies
% phase is list of phases or 'sin', 'cos', 'alt'
%
% defaults:
%  phase = sin
%  dBSPL=56 dB SPL

dt=globalStimParams.dt;
frequencies=stimComponents.frequencies;
toneDuration=stimComponents.toneDuration;
dBSPL=stimComponents.amplitudesdB;
phases=stimComponents.phases;

if ischar(phases)
    switch phases
        case 'sin'
            phases= zeros(1,length(frequencies));
        case 'cos'
            phases= pi/2*ones(1,length(frequencies));
        case 'alt'
            phases= repmat([0 pi/2], 1, floor(length(frequencies)/2));
            if length(phases)<length(frequencies)
                phases=[phases 0];
            end
        case {'ran', 'rand'}
            phases= 2*pi*rand(1,length(frequencies));
    end
end

if length(phases)==1, phases=repmat(phases(1), 1, length(frequencies)); end
if length(phases)<length(frequencies)
    error('makeTone:phase specification must have the same length as frequencies')
end

if length(dBSPL)==1, dBSPL=repmat(dBSPL(1), 1, length(frequencies)); end
if length(dBSPL)<length(dBSPL)
    error('makeTone:dBSPL specification must have the same length as frequencies')
end

time=dt:dt:toneDuration;
amplitudes=28e-6* 10.^(dBSPL/20);

tone=zeros(size(time));
for i=1:length(frequencies)
    frequency=frequencies(i);
    phase=phases(i);
    tone=tone+amplitudes(i)*sin(2*pi*frequency*time+phase);
end

% figure(1), plot(time, signal)


%---------------------------------------------------- makeOHIOtones
function stimulus=makeOHIOtones(globalStimParams, stimComponents)

% Generates a stimulus consisting of one or more 10-ms tones
%  The length of the list of frequencies determines the numberof tones
% Tones are either presented at 10-ms intervals or simultaneously
% all tones are individually ramped
% Each tone has its own amplitude and its own ramp.

frequencies=stimComponents.frequencies;
amplitudesdB=stimComponents.amplitudesdB;
nFrequencies=length(frequencies);

if amplitudesdB==-100
    % catch trial
    amplitudesdB=repmat(-100,1,nFrequencies);
end

dt=globalStimParams.dt;

toneDuration=.010;
time=dt:dt:toneDuration;

% fixed ramp, silenceDuration, toneDuration
rampDuration=0.005;
rampTime=dt:dt:rampDuration;
ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
    ones(1,length(time)-length(rampTime))];
ramp=ramp.*fliplr(ramp);

silenceDuration=0.010;
silenceDurationLength=round(silenceDuration/dt);
initialSilence=zeros(1,silenceDurationLength);

silenceToneDuration=toneDuration + silenceDuration;
silenceToneDurationLength=round(silenceToneDuration/dt);

% OHIO spect produces simultaneous tones
if strcmp(stimComponents.OHIOtype,'OHIOspect')
    totalDuration=silenceToneDuration;
else
    totalDuration=silenceToneDuration*nFrequencies;
end

totalDurationLength=round(totalDuration/dt);
stimulus=zeros(1,totalDurationLength);
toneBeginPTR=1;

for i=1:nFrequencies
    frequency=frequencies(i);
    dBSPL=amplitudesdB(i);
    amplitude=28e-6* 10.^(dBSPL/20);
    tone=amplitude*sin(2*pi*frequency*time);
    tone=tone.*ramp;
    % stimulus is normally zeros except for OHIOspect
    stimulus(toneBeginPTR:toneBeginPTR+silenceToneDurationLength-1)=...
        [initialSilence tone]+...
        stimulus(toneBeginPTR:toneBeginPTR+silenceToneDurationLength-1);    
    if ~strcmp(stimComponents.OHIOtype,'OHIOspect')
        toneBeginPTR=toneBeginPTR+silenceToneDurationLength;
    end
end
% figure(2), plot( stimulus')


%---------------------------------------------------- makeFMtone
function tone=makeFMtone(globalStimParams, stimComponents)
% maketone generates a stimComponents tone
% tone=maketone(dt, frequencies, toneDuration, dBSPL, phases)
% tone is returned in Pascals
% frequencies is a list of frequencies
% phase is list of phases or 'sin', 'cos', 'alt'
%
% defaults:
%  phase = sin
%  dBSPL=56 dB SPL

dt=globalStimParams.dt;
frequencies=stimComponents.frequencies;
toneDuration=stimComponents.toneDuration;
dBSPL=stimComponents.amplitudesdB;
phases=stimComponents.phases;
fmDepth=stimComponents.fmDepth;
fmFrequency=stimComponents.fmFrequency;

if ischar(phases)
    switch phases
        case 'sin'
            phases= zeros(1,length(frequencies));
        case 'cos'
            phases= pi/2*ones(1,length(frequencies));
        case 'alt'
            phases= repmat([0 pi/2], 1, floor(length(frequencies)/2));
            if length(phases)<length(frequencies)
                phases=[phases 0];
            end
    end
end

if length(phases)==1, phases=repmat(phases(1), 1, length(frequencies)); end
if length(phases)<length(frequencies)
    error('makeTone:phase specification must have the same length as frequencies')
end

if length(dBSPL)==1, dBSPL=repmat(dBSPL(1), 1, length(frequencies)); end
if length(dBSPL)<length(dBSPL)
    error('makeTone:dBSPL specification must have the same length as frequencies')
end

time=dt:dt:toneDuration;
amplitudes=28e-6* 10.^(dBSPL/20);

tone=zeros(size(time));
for i=1:length(frequencies)
    frequency=frequencies(i);
    phase=phases(i);
    tone=tone+amplitudes(i)*sin(2*pi*frequency*time+phase + fmDepth*sin(2*pi*fmFrequency*time));
end

% figure(1), plot(time, signal)

%----------------------------------------------------makeTransposedStimulus
function [transposedStimulus]=makeTransposedStimulus(globalStimParams, stim)
% globalStimParams.FS=100000;
% globalStimParams.overallDuration=.1;  % s
% globalStimParams.doPlot=1;
%
% stim.type='transposedStimulus';
% stim.phases='sin';
% stim.toneDuration=.05;;
% stim.frequencies=500;
% stim.amplitudesdB=50;
% stim.beginSilence=.01;
% stim.endSilence=-1;
% stim.rampOnDur=.002;
% stim.rampOffDur=-1;
%
% stim.transCarrierFreq=4000;
% stim.transModFreq=100;

transposedStimulus=[];
% make envelope of transposed tone
for i=1:length(stim.transModFreq)
    stim.type='tone';
    stim.frequencies=stim.transModFreq(i);
    stim.endsilence=-1;     stim.beginSilence=0;
    [envelope, msg]=stimulusCreate(globalStimParams, stim);    % NB recursive
    envelope=envelope(:,1);   % mono
    % HW rectify
    envelope(find(envelope<0))=0;
    % LP filter
    maxEnvelope=max(envelope);
    envelope=UTIL_Butterworth (envelope, globalStimParams.dt, 10, .2*stim.transModFreq(i), 2);
    envelope=envelope*maxEnvelope/max(envelope);

    % make the carrier
    stim.frequencies=stim.transCarrierFreq(i);
    stim.endsilence=-1;    stim.beginSilence=0;
    [audio, msg]=stimulusCreate(globalStimParams, stim);
    carrier=audio(:,1);
    x= (carrier.*envelope)';
    % base amplitude on peak of unmodulated carrier
    x=x/max(carrier);
    transposedStimulus=[transposedStimulus; x];
end
transposedStimulus=sum(transposedStimulus,1);
% clf,plot(transposedStimulus)

%--------------------------------------------------------------------makeClicks
function clickTrain=makeClicks(globalStimParams, stimComponents)
% makeClicks(F0, clickTimes, duration, FS);
% F0 specifies the repetition rate of the click sequence
% If F0=-1, a single click is generated at the start of the duration of the signal
%
% clickTimes a are fractions of the period
%  and specify when the click appears in the period
% A click time of 1 is reset to zero.
% if the clicktime plus the click width is greater than the period, no click is generated
% clicks are treated as 20 microPascal high before amplification
%  unless otherwise specified in stimComponents.clickHeight
% click width is dt unless specified in stimComponents.clickWidth
%
% for regular pulse train set clicktimes=1 or zero;
% FS is the sampling interval;
% CarlyonClickTrain(100, [.4 1], 40, 441000);

FS=globalStimParams.FS; % sample rate
dt=1/FS;

try,clickWidth=stimComponents.clickWidth;catch, clickWidth=dt; end
try,clickHeight=stimComponents.clickHeight; catch, clickHeight=28e-6 * 10^(stimComponents.amplitudesdB/20); end
try, clickTimes=stimComponents.clickTimes; catch, clickTimes=1; end

% clickTimes are the times in a cycle that the click
% occurs
checkClickTimes=clickTimes-1;
if max(checkClickTimes)>1
    msg= 'click times must be <= than 1 (period)';
    return
end

if clickWidth>=1/stimComponents.clickRepeatFrequency
    msg= 'click width is too great for frequency';
    return
end

duration=stimComponents.toneDuration;
totalLength=round(stimComponents.toneDuration/dt);
F0=stimComponents.clickRepeatFrequency;
F0=round(F0/dt)*dt;
if F0==-1 % single click required
    F0=1/duration;
    repetitions=1;
    clickStartTimes=1; %clicktimes are fractions of a period
else
    repetitions=round(F0*duration)-1;
    duration=repetitions/F0;
    clickStartTimes=clickTimes;
end
% if a clickTime=1 (end of duty period) set it to the beginning
clickStartTimes(clickStartTimes==1)=0;

period=1/F0;
time=dt:dt:period;
nPoints=length(time);
signal=zeros(1,nPoints);
dBSPL=stimComponents.amplitudesdB;

% compute click train for a single cycle
clickWidthNpoints=round(clickWidth*FS);
for i=1:length(clickStartTimes)
    %     if clickStartTimes(i)<clickWidth
    %         clickStartTimes(i)=dt;
    %     end
    clickTime=round(period*clickStartTimes(i)/dt -dt);
    % clicks are treated as 20 microPascal high
    if clickTime+clickWidthNpoints<length(signal)
        signal(clickTime+1:clickTime+clickWidthNpoints)=clickHeight;
    end
end

clickTrain=repmat(signal, 1, repetitions);

if length(clickTrain)>totalLength
    clickTrain=clickTrain(1:totalLength);
elseif length(clickTrain)<totalLength
    timeToAdd=zeros(1,round((totalLength-length(clickTrain))));
    clickTrain=[clickTrain timeToAdd];
    % figure (1), plot(clickTrain)
end

%----------------------------------------------------------------makePressnitzerClicks
function signal=makePressnitzerClicks(globalStimParams, stimComponents)
% PressnitzerClicks(k,duration,dt)
% Generates a sequence of clicks with intervals kxx
%  where x= rand*k/2
% This is not the same as Kaernbach and Demany clicks

FS=globalStimParams.FS; % sample rate
dt=1/FS;

% obligatory parameter
try
    k=stimComponents.k;
catch
    error('PressnitzerClicks: field ''k'' is missing from stimcomponent')
end

% optional parameters
if isfield(stimComponents,'clickWidth')
    clickWidth=stimComponents.clickWidth;
else
    clickWidth=dt;
end
clickWidthNpoints=round(clickWidth*FS);

if isfield(stimComponents,'clickHeight')
    clickHeight=stimComponents.clickHeight;
else
    clickHeight=28e-6 * 10^(stimComponents.amplitudesdB/20);
end

duration=stimComponents.toneDuration;

signalLength=round(duration/dt);
signal=zeros(1,signalLength);
kInterval=round(k/dt);
halfk=k/2;
signal(1)=clickHeight;
timeIdx=0;
while timeIdx<signalLength
    % first interval = k
    clickTime=timeIdx+kInterval;
    signal(clickTime:clickTime+clickWidthNpoints)=clickHeight;
    timeIdx=timeIdx+kInterval;

    % second interval = 0 : k/2
    intervalx1=round(rand*halfk/dt);
    clickTime=timeIdx+intervalx1;
    signal(clickTime:clickTime+clickWidthNpoints)=clickHeight;
    timeIdx=timeIdx+intervalx1;

    % third interval = 0 : k/2
    intervalx1=round(rand*halfk/dt);
    clickTime=timeIdx+intervalx1;
    signal(clickTime:clickTime+clickWidthNpoints)=clickHeight;
    timeIdx=timeIdx+intervalx1;
end

signal=signal(1:signalLength);
% figure(1),	plot(dt:dt:duration,signal)

%----------------------------------------------------------------makePressnitzerABXClicks
function signal=makePressnitzerABxClicks(globalStimParams, stimComponents)
% Generates a sequence of clicks with intervals ABx
% AB interval is 2*k
% where A= rand* k
%       B= k-A
%       x= k/2
% These are second order clicks

FS=globalStimParams.FS; % sample rate
dt=1/FS;

% obligatory parameter
try
    k=stimComponents.k;
catch
    error('PressnitzerClicks: field ''k'' is missing from stimcomponent')
end

% optional parameters
if isfield(stimComponents,'clickWidth')
    clickWidth=stimComponents.clickWidth;
else
    clickWidth=dt;
end
clickWidthNpoints=round(clickWidth*FS);

if isfield(stimComponents,'clickHeight')
    clickHeight=stimComponents.clickHeight;
else
    clickHeight=28e-6 * 10^(stimComponents.amplitudesdB/20);
end

duration=stimComponents.toneDuration;

signalLength=round(duration/dt);
signal=zeros(1,2*signalLength); % allow for overrun
ABinterval=k/dt;                % i.e. the number of dt steps
randomInterval=ABinterval/2;
signal(1)=clickHeight;
time=0;
while time<signalLength
    % first interval = A
    intervalA=rand*ABinterval;
    clickTime=round(time+intervalA)+1;   % +1 avoids zero index
    signal(clickTime:clickTime+clickWidthNpoints)=clickHeight;
    time=time+intervalA;

    % second interval = B
    intervalB=ABinterval-intervalA;
    clickTime=round(time+intervalB)+1;
    signal(clickTime:clickTime+clickWidthNpoints-1)=clickHeight;
    time=time+intervalB;

    % third interval = 0 : k/2
    intervalx1=rand*randomInterval;    % mean random interval=k
    clickTime=round(time+intervalx1)+1;
    signal(clickTime:clickTime+clickWidthNpoints-1)=clickHeight;
    time=time+intervalx1;
end

signal=signal(1:signalLength);
% figure(1),	plot(dt:dt:duration,signal)

%-----------------------------------------------------makeABxClicks
function signal=makeABxClicks(globalStimParams, stimComponents)
% Generates a sequence of clicks with intervals ABx
% AB interval is 2*k
% where A= rand* k
%       B= k-A
%       x= rand*2*k
% These are second order clicks

FS=globalStimParams.FS; % sample rate
dt=1/FS;

% obligatory parameter
try
    k=stimComponents.k;
catch
    error('PressnitzerClicks: field ''k'' is missing from stimcomponent')
end

% optional parameters
if isfield(stimComponents,'clickWidth')
    clickWidth=stimComponents.clickWidth;
else
    clickWidth=dt;
end
clickWidthNpoints=round(clickWidth*FS);

if isfield(stimComponents,'clickHeight')
    clickHeight=stimComponents.clickHeight;
else
    clickHeight=28e-6 * 10^(stimComponents.amplitudesdB/20);
end

duration=stimComponents.toneDuration;

signalLength=round(duration/dt);
signal=zeros(1,2*signalLength); % allow for overrun
ABinterval=2*k/dt;
randomInterval=ABinterval;
signal(1)=clickHeight;
timeIdx=0;
while timeIdx<signalLength
    % first interval = A
    intervalA=round(rand*ABinterval);
    clickTime=timeIdx+intervalA+1;
    signal(clickTime:clickTime+clickWidthNpoints-1)=clickHeight;
    timeIdx=timeIdx+intervalA;

    % second interval = B
    intervalB=round(ABinterval-intervalA);
    clickTime=timeIdx+intervalB;
    signal(clickTime:clickTime+clickWidthNpoints-1)=clickHeight;
    timeIdx=timeIdx+intervalB;

    % third interval = 0 : k
    intervalx1=round(rand*randomInterval);    % mean random interval=k
    clickTime=timeIdx+intervalx1;
    signal(clickTime:clickTime+clickWidthNpoints-1)=clickHeight;
    timeIdx=timeIdx+intervalx1;
end

signal=signal(1:signalLength);
% figure(1),	plot(dt:dt:duration,signal)

%----------------------------------------------------------------makeYostClicks
function signal=makeYostClicks(globalStimParams, stimComponents)
% Generates a shuffled sequence of clicks with intervals kxxxx
%  where max(x)= 2*k
%  and there are n occurrences of x
% this section requires:
%  stimComponents.k
%  stimComponents.nXs
%  stimComponents.toneDuration
% optional:
%  stimComponents.clickWidth	%useful because width depends on dt
%  stimComponents.clickHeight	%best left to amplitude rescaling later

FS=globalStimParams.FS; % sample rate
dt=1/FS;

% obligatory parameters
try
    k=stimComponents.k;
catch
    error('makeYostClicks: field ''k'' is missing from stimComponents')
end

try
    nXs=stimComponents.nXs;
catch
    error('makeYostClicks: field ''nXs'' is missing from stimComponents')
end

try
    shuffled=stimComponents.shuffled;
catch
    error('makeYostClicks: field ''shuffled'' is missing from stimComponents')
end

try
    duration=stimComponents.toneDuration;
catch
    error('makeYostClicks: field ''toneDuration'' is missing from stimComponents')
end

% optional parameters
if isfield(stimComponents,'clickWidth')
    clickWidth=stimComponents.clickWidth;
else
    clickWidth=dt;
end
clickWidthNpoints=round(clickWidth*FS);

if isfield(stimComponents,'clickHeight')
    clickHeight=stimComponents.clickHeight;
else
    clickHeight=28e-6 * 10^(stimComponents.amplitudesdB/20);
end

kInterval=round(k/dt);
twoK=k*2;				% max width of x interval

signalLength=round(duration/dt);
signal=zeros(1,signalLength);

timeIdx=0;
intervalCount=0;
while timeIdx<signalLength
    timeIdx=timeIdx+kInterval;
    if timeIdx>signalLength, break,end
    intervalCount=intervalCount+1;
    intervals(intervalCount)=kInterval;

    % repeat x intervals as required
    if nXs>0
        for nX=1:nXs
            xInterval=round(rand*twoK/dt);
            timeIdx=timeIdx+xInterval;
            if timeIdx>signalLength, break,end
            intervalCount=intervalCount+1;
            intervals(intervalCount)=xInterval;
        end
    end
    if timeIdx>signalLength, break,end
end

% shuffle intervals
if shuffled
    randomNumbers=rand(1,length(intervals));
    [randomNumbers idx]=sort(randomNumbers);
    intervals=intervals(idx);
    idx=intervals>0;
    intervals=intervals(idx);
end

intervalCount=length(intervals);
signal(1)=clickHeight;
clickTime=0;
for i=1:intervalCount
    clickTime=clickTime+intervals(i);
    signal(clickTime:clickTime+clickWidthNpoints)=clickHeight;
end
signal=signal(1:signalLength);
%  figure(1),	plot(dt:dt:duration,signal)

%--------------------------------------------------------------------makeKxxClicks
function signal=makeKxxClicks(globalStimParams, stimComponents)
% Click train consists of kkxxx.. sequences
% k is the duration of a fixed interval (seconds)
% random intervals are distributed 0 : 2* k (NB not like Pressnitzer clicks)
% nKs is the number of successive 'k' intervals
% nXs is the number of random intervals between k sequences
% sequences of 3 x intervals > k are replaced with new sequences
% shuffled causes all intervals to be reordered randomly
% NB 1/k is the mean click rate

FS=globalStimParams.FS; % sample rate
dt=1/FS;

try
    k=stimComponents.k;		% duration (s) of fixed intervals
catch
    error('makeYostClicks: field ''k'' is missing from stimComponents')
end

try
    duration=stimComponents.toneDuration;
catch
    error('makeYostClicks: field ''duration'' is missing from stimComponents')
end

if isfield(stimComponents,'clickWidth')
    clickWidth=stimComponents.clickWidth;
else
    clickWidth=dt;
end

if isfield(stimComponents,'clickHeight')
    clickHeight=stimComponents.clickHeight;
else
    clickHeight=28e-6 * 10^(stimComponents.amplitudesdB/20);
end


if isfield(stimComponents,'order')
    order=stimComponents.order;
else
    order=1;
end

if isfield(stimComponents,'nKs')
    nKs=stimComponents.nKs;
else
    nKs=1;
end

if isfield(stimComponents,'nXs')
    nXs=stimComponents.nXs;
else
    nXs=1;
end

if isfield(stimComponents,'shuffled')
    shuffled=stimComponents.shuffled;
else
    shuffled=1;
end

kLength=round(k/dt);    % fixed interval
xLength=2*kLength;      % maximum random interval
requiredSignalLength=round(duration/dt);
intervalsPerCycle=(nKs+nXs);
cycleLength=nKs*kLength+nXs*xLength;
% more cycles to allow for uncertainty
nCycles=5*round(requiredSignalLength/cycleLength);
nIntervals=nCycles*intervalsPerCycle;

% random intervals
if nXs>0
    xIntervals=floor(rand(1,nIntervals)*2*kLength);
    % weed out triple intervals > 2*k
    rogues=1;
    while sum(rogues)
        y=(xIntervals>kLength);
        rogues=(sum([y(1:end-2)' y(2:end-1)' y(3:end)'],2)>2);
        xIntervals(rogues)=floor(rand*2*kLength);
    end
    xIntervals=reshape(xIntervals,nCycles,intervalsPerCycle);
else
    xIntervals=[];
end

% insert constrained (k) intervals
if nKs>0
    switch order
        case 1
            kIntervals=floor(ones(nCycles,nKs)*kLength);
        case 2
            nKs=1; % force this
            kIntervals=floor(rand(nCycles,1)*kLength);
            kIntervals=[kIntervals kLength-kIntervals];
    end
else
    kIntervals=[];
end

% combine fixed and random
intervals=[kIntervals xIntervals(:,nKs+1:end)];
% make a single array;
[r c]=size(intervals);
intervals=reshape(intervals',1,r*c);

% shuffle intervals
if shuffled
    randomNumbers=rand(1,length(intervals));
    [randomNumbers idx]=sort(randomNumbers);
    intervals=intervals(idx);
    idx=intervals>0;
    intervals=intervals(idx);
end

% convert intervals to clicks
clickTimes=cumsum(intervals);
signal(clickTimes)=clickHeight;
signal=signal(1:requiredSignalLength);
% figure(1), clf, plot(signal)



%--------------------------------------------------------------------makeNoise
function noise=makeNoise(globalStimParams, stimComponents)
% FS in Hz, noiseDuration in s, delay in s;
% noise is returned with overall level dB(rms) = amplitudesdB
%
% % You need
%
% stim.type='noise'; % or 'IRN', or 'pinkNoise'
% stim.toneDuration=.05;
% stim.amplitudesdB=50;
% stim.beginSilence=.01;
% stim.endSilence=-1;
% stim.rampOnDur=.002;
% stim.rampOffDur=-1;
%
% % Mandatory structure fields
% globalStimParams.FS=100000;
% globalStimParams.overallDuration=.1;  % s
% globalStimParams.doPlot=1;
% globalStimParams.doPlay=1;
%
% [audio, msg]=stimulusCreate(globalStimParams, stim, );
%
% % local
% stim.type='noise'; % or 'IRN'
%
FS=globalStimParams.FS;
noiseDuration= stimComponents.toneDuration;
npts=round(noiseDuration*FS);
noise=randn(1,npts);    % NB randn (normally distributed)

switch stimComponents.type
    case  'pinkNoise'
        %         noise=UTIL_lowpassFilterFreq(noise, 100, 1/FS);
        noise=UTIL_bandPassFilter(noise, 1, 100, 200, 1/FS,[]);
end

rms=(mean(noise.^2)^.5);  %should be 20 microPascals for 0 dB SPL
adjust=20e-6/rms;
noise=noise*adjust;
rms=(mean(noise.^2)^.5);
amplitude=10.^(stimComponents.amplitudesdB/20);
noise=amplitude*noise;
% rms=(mean(noise.^2)^.5);
% dBnoise=20*log10(rms/20e-6)


%--------------------------------------------------------------------makeWhiteNoise
function noise=makeWhiteNoise(globalStimParams, stimComponents)
% FS in Hz, noiseDuration in s, delay in s;
% noise is bandpass filtered between 100 and 10000 Hz
% spectrum level (dB/Hz) is 40 dB below nominal level.
% noise is returned with dB(rms) = amplitudesdB
%
% % You need
%
% stim.type='noise'; % or 'IRN', or 'pinkNoise'
% stim.toneDuration=.05;
% stim.amplitudesdB=50;
% stim.beginSilence=.01;
% stim.endSilence=-1;
% stim.rampOnDur=.002;
% stim.rampOffDur=-1;
%
% % Mandatory structure fields
% globalStimParams.FS=100000;
% globalStimParams.overallDuration=.1;  % s
% globalStimParams.doPlot=1;
% globalStimParams.doPlay=1;
%
% [audio, msg]=stimulusCreate(globalStimParams, stim, );
%
% % local
% stim.type='noise'; % or 'IRN'
%
FS=globalStimParams.FS;
noiseDuration= stimComponents.toneDuration;
npts=round(noiseDuration*FS);
noise=randn(1,npts);

noise=UTIL_bandPassFilter (noise, 6, 100, 10000, 1/FS, []);

rms=(mean(noise.^2)^.5);  %should be 20 microPascals for 0 dB SPL
adjust=20e-6/rms;
noise=noise*adjust;
rms=(mean(noise.^2)^.5);
amplitude=10.^(stimComponents.amplitudesdB/20);
noise=amplitude*noise;
% rms=(mean(noise.^2)^.5);
% dBnoise=20*log10(rms/20e-6)


%-----------------------------------------------------------------makeIRN
function noise=makeIRN(globalStimParams, stimComponents)
% FS in Hz, noiseDuration in s, delay in s;
% noise is returned with dB(rms) = amplitudesdB
%
% % You need
%
% stim.type='noise'; % or 'IRN', or 'pinkNoise'
% stim.toneDuration=.05;
% stim.amplitudesdB=50;
% stim.beginSilence=.01;
% stim.endSilence=-1;
% stim.rampOnDur=.002;
% stim.rampOffDur=-1;
%
% % Mandatory structure fields
% globalStimParams.FS=100000;
% globalStimParams.overallDuration=.1;  % s
% globalStimParams.doPlot=1;
% globalStimParams.doPlay=1;
%
% [audio, msg]=stimulusCreate(globalStimParams, stim, );
%
% % local
% stim.type='noise'; % or 'IRN'
% % for IRN only
% stim.niterations = 8;   %0 for white noise
% stim.delay = 1/150;
% stim.irnGain = 1;
%
FS=globalStimParams.FS;
noiseDuration= stimComponents.toneDuration;

nIterations=stimComponents.niterations;
if nIterations==0
    % white noise is specified as nIterations=1
    nIterations=1;
    IRNgain=0;
    delay=0.01; % dummy
else
    % IRN
    delay=stimComponents.delay;
    IRNgain=stimComponents.irnGain;
end

npts=round(noiseDuration*FS);
dels=round(delay*FS);
noise=randn(1,npts);

%fringe=nIterations*dels;
%npts=npts+fringe;

for i=1:nIterations,
    dnoise=[noise(dels+1:npts) noise(1:dels)];
    dnoise=dnoise.*IRNgain;
    noise=noise+dnoise;
end;

switch stimComponents.type
    case  'pinkNoise'
        noise=UTIL_lowpassFilterFreq(noise, 10000, 1/FS);
end

rms=(mean(noise.^2)^.5);  %should be 20 microPascals for 0 dB SPL
adjust=20e-6/rms;
noise=noise*adjust;
rms=(mean(noise.^2)^.5);
amplitude=10.^(stimComponents.amplitudesdB/20);
noise=amplitude*noise;
% rms=(mean(noise.^2)^.5);
% dBnoise=20*log10(rms/20e-6)

%------------------------------------------------------------------ makeRPN
function RPN=makeRPN(globalStimParams, stim)
% 'period' is a collection of samples - AAABCD
% you need
%
% stim.type='RPN';
% stim.toneDuration=.2;
% stim.amplitudesdB=50;
% stim.beginSilence=.01;
% stim.rampOnDur=.002;
% stim.rampOffDur=-1;
%
% stim.sampleDuration=.005;  %200 Hz pitch
% stim.nSimilarSamples=5;   % pitch strength
% stim.nIndependentSamples=1% dilutes strength
%
% % Mandatory structure fields
%  globalStimParams.FS=44100;
%  globalStimParams.overallDuration=.21;  % s
%
% globalStimParams.doPlot=1;
% globalStimParams.doPlay=1;
% [audio, msg]=stimulusCreate(globalStimParams, stim);

FS=globalStimParams.FS;
ptsPerSample=floor(stim.sampleDuration*FS);

samplesPerPeriod=stim.nSimilarSamples+stim.nIndependentSamples;
periodDuration=samplesPerPeriod*stim.sampleDuration;

totalNumPeriods=2*floor(stim.toneDuration/periodDuration);  % longer than necessary
if totalNumPeriods<1
    error('stimulusCreate: RPN, stimulus duration needs to be longer')
end

RPN=[];
for j=1:totalNumPeriods
    noise=randn(1,ptsPerSample);
    for i=1:stim.nSimilarSamples
        RPN=[RPN noise];
    end

    for i=1:stim.nIndependentSamples
        noise=randn(1,ptsPerSample);
        RPN=[RPN noise];
    end
end

targetStimulusLength=round(stim.toneDuration/FS);
RPN=RPN(1:floor(stim.toneDuration*FS));     % take enough for stimulus

rms=(mean(RPN.^2)^.5);  %should be 20 microPascals for 0 dB SPL
adjust=20e-6/rms;
RPN=RPN*adjust;
rms=(mean(RPN.^2)^.5);
amplitude=10.^(stim.amplitudesdB/20);
RPN=amplitude*RPN;
% rms=(mean(noise.^2)^.5);
% dBnoise=20*log10(rms/20e-6)

%--------------------------------------------------------------------applyRampOn
function signal=applyRampOn(signal, rampDur, rampOnTime, sampleRate)
%applyRampOn applies raised cosine ramp
%rampOntime is the time at which the ramp begins
%At all other times the mask has a value of 1
%signal=applyRampOn(signal, rampDur, rampOnTime, sampleRate)

rampDurPoints=round(rampDur*sampleRate);
rampOn= (1+cos(pi:pi/(rampDurPoints-1): 2*pi))/2';

sigDurPoints=length(signal);
mask(1:sigDurPoints)=1;
rampOnStartIndex=round(rampOnTime*sampleRate+1);
mask(rampOnStartIndex: rampOnStartIndex+ rampDurPoints-1)=rampOn;
signal=signal.*mask;
%plot(mask)

%--------------------------------------------------------------------applyRampOff
function signal=applyRampOff(signal, rampDur, rampOffTime, sampleRate)
%applyRampOn applies raised cosine squared ramp
%rampOffTime is the time at which the ramp begins
%At all other times the mask has a value of 1
% signal=applyRampOff(signal, rampDur, rampOffTime, sampleRate)

rampDurPoints=round(rampDur*sampleRate);
rampOff= (1+cos(0:pi/(rampDurPoints-1): pi))/2';

sigDurPoints=length(signal);
mask=ones(1,sigDurPoints);
rampOffStartIndex=round(rampOffTime*sampleRate+1);
mask(rampOffStartIndex: round(rampOffStartIndex+ rampDurPoints-1))=rampOff;
if length(mask)>sigDurPoints, mask=mask(1:sigDurPoints); end
signal=signal.*mask;
%plot(mask)

function signal=applyGaussianRamps(signal, sigma, sampleRate)
dt=1/sampleRate;
time=dt:dt:dt*length(signal);
ramp=1-exp(-time.^2/(2*sigma^2));
% apply onset ramp
signal=signal.*ramp;
% apply offset ramp
ramp=fliplr(ramp);
signal=signal.*ramp;



%--------------------------------------------------------------------checkDescriptors
function [globalStimParams, stimComponents]=...
    checkDescriptors(globalStimParams, stimComponents)

try
    % if FS exists, it takes priority
    globalStimParams.dt=1/globalStimParams.FS;
catch
    % otherwise set FS using dt
    globalStimParams.FS=1/globalStimParams.dt;
end

sampleRate=globalStimParams.FS;

globalStimParams.nSignalPoints=...
    round(globalStimParams.overallDuration* sampleRate);

% optional field (ears)
try
    globalStimParams.ears;
catch
    % default: dichotic.
    globalStimParams.ears='dichotic';
end

% audioOutCorrection is optional
% audioOutCorrection is a scalar for reducing the sound
try
    globalStimParams.audioOutCorrection;
catch
    %      set to 1 if omitted
    globalStimParams.audioOutCorrection=1;
end

try
    globalStimParams.doPlay;
catch
    % default plays sound only if explicitly requested
    globalStimParams.doPlay=0;
end

try
    globalStimParams.doPlot;
catch
    % no plotting unless explicitly requested
    globalStimParams.doPlot=0;
end

[ears nComponentSounds]=size(stimComponents);

for ear=1:2 % 1=left/ 2=right

    % create a list of components whose type is specified
    % if no type is specified assume that it is an empty specification
    % this is allowed
    validStimComponents=[];
    for i=1:nComponentSounds
        try
            if ~isempty(stimComponents(ear,i).type)
                validStimComponents=[validStimComponents i];
            end
        catch
        end
    end

    for componentNo=validStimComponents
        % If no AM filed is present, create it for completeness
        if ~isfield(stimComponents(ear,componentNo),'AMfrequency') |...
                ~isfield(stimComponents(ear,componentNo),'AMdepth')
            stimComponents(ear,componentNo).AMfrequency=0;
            stimComponents(ear,componentNo).AMdepth=0;
        end

        % all signals must have durations, amplitudes and ramps
        if ...
                isempty(stimComponents(ear,componentNo).type) |...
                isempty(stimComponents(ear,componentNo).toneDuration) |...
                isempty(stimComponents(ear,componentNo).amplitudesdB) |...
                isempty(stimComponents(ear,componentNo).rampOnDur)
            descriptorError( 'missing stimComponent descriptor', stimComponents, ear, componentNo)
        end

        try, stimComponents(ear,componentNo).endSilence; catch, stimComponents(ear,componentNo).endSilence=-1; end

        % ramp checks do not apply to file input
        if ~strcmp(stimComponents(ear,componentNo).type, 'file')
            % match offset ramp to onset if not explicitly specified
            if stimComponents(ear,componentNo).rampOffDur==-1,
                stimComponents(ear,componentNo).rampOffDur=stimComponents(ear,componentNo).rampOnDur;
            end
            % ramps must be shorter than the stimulus
            if stimComponents(ear,componentNo).rampOffDur> stimComponents(ear,componentNo).toneDuration | ...
                    stimComponents(ear,componentNo).rampOnDur> stimComponents(ear,componentNo).toneDuration
                descriptorError( 'ramp longer than sound component', stimComponents, ear, componentNo)
            end
        end

        % end silence is measured to fit into the global duration
        if stimComponents(ear,componentNo).endSilence==-1,
            stimComponents(ear,componentNo).endSilence=...
                globalStimParams.overallDuration-...
                stimComponents(ear,componentNo).beginSilence -...
                stimComponents(ear,componentNo).toneDuration;
            
            endSilenceNpoints=stimComponents(ear,componentNo).endSilence...
                *sampleRate;
        end
        if stimComponents(ear,componentNo).endSilence<0
            globalStimParams
            descriptorError( 'component durations greater than overallDuration', stimComponents, ear, componentNo)
        end

        % check overall duration of this component against global duration
        totalDuration= ...
            stimComponents(ear,componentNo).beginSilence+...
            stimComponents(ear,componentNo).toneDuration+...
            stimComponents(ear,componentNo).endSilence;

        % avoid annoying error message for single stimulus component
        if ears==1 && nComponentSounds==1
            globalStimParams.overallDuration= totalDuration;
        end

        % check total duration
        totalSignalPoints= round(totalDuration* sampleRate);
        if totalSignalPoints  >globalStimParams.nSignalPoints
            descriptorError( 'Signal component duration does not match overall duration ', stimComponents, ear, componentNo)
        end

        if isfield(stimComponents(ear,componentNo), 'filter')
            if ~isequal(length(stimComponents(ear,componentNo).filter), 3)
                descriptorError( 'Filter parameter must have three elements ', stimComponents, ear, componentNo)
            end
        end
    end % component
    % ??
    if strcmp(globalStimParams.ears,'monoticL') | strcmp(globalStimParams.ears, 'monoticR'), break, end
end		% ear


%-------------------------------------------------------------------- descriptorError
function descriptorError( msg, stimComponents, ear, componentNo)
stimComponents(ear, componentNo)

disp(' *********** **************** ************')
disp([ '...Error in stimComponents description: '])
disp([msg ])
disp(['Ear = ' num2str(ear) ' component No ' num2str(componentNo)])
disp(' *********** **************** ************')
error('myError ')


%-------------------------------------------------------------------- normalize
function [normalizedSignal, gain]= normalize(signal)
% normalize (signal)
maxSignal=max(max(signal));
minSignal=min(min(signal));
if -minSignal>maxSignal, normalizingFactor=-minSignal; else normalizingFactor=maxSignal; end
normalizingFactor=1.01*normalizingFactor;
gain= 20*log10(normalizingFactor);
normalizedSignal=signal/normalizingFactor;


%--------------------------------------------------------------------Butterworth
function y=Butterworth (x, dt, fl, fu, order)
% Butterworth (x, dt, fu, fl, order)
% Taken from Yuel and beauchamp page 261
% NB error in their table for K (see their text)
% x is original signal
% fu, fl upper and lower cutoff
% order is the number of times the filter is applied

q=(pi*dt*(fu-fl));
J=1/(1+ cot(q));
K= (2*cos(pi*dt*(fu+fl)))/(1+tan(q)*cos(q));
L= (tan(q)-1)/(tan(q)+1);
b=[J -J];
a=[1 -K  -L];
for i=1:order
    y=filter(b, a, x);
    x=y;
end


% -------------------------------------------------------- UTIL_amp2dB
function [y] = UTIL_amp2dB (x, ref)
% Calculates a dB (ref. ref) value 'y' from a peak amplitude number 'x'.
% if ref omitted treat as dB
% Check the number of arguments that are passed in.
if (nargin < 2)
    ref=1;
end
if (nargin > 2)
    error ('Too many arguments');
end

% Check arguments.
if x < 0.0
    error ('Can not calculate the log10 of a negative number');
elseif x == 0.0
    warning ('log10 of zero.  The result is set to -1000.0');
    y = -1000.0;
else
    % Do calculations.
    y = 20.0 * log10(x/(sqrt(2)*ref));

end

%-------------------------------------------------------------------- FFT
function showFFT (getFFTsignal, dt)
color='r';
figure(2), clf,
hold off

% trim initial silence
idx=find(getFFTsignal>0);
if ~isempty(idx)
    getFFTsignal=getFFTsignal(idx(1):end);
end
%trim final silence
getFFTsignal=getFFTsignal(end:-1:1);
idx=find(getFFTsignal>0);
if ~isempty(idx)
    getFFTsignal=getFFTsignal(idx(1):end);
    getFFTsignal=getFFTsignal(end:-1:1);
end

% Analyse make stimComponents length a power of 2
x=length(getFFTsignal);
squareLength=2;
while squareLength<=x
    squareLength=squareLength*2;
end
squareLength=round(squareLength/2);
getFFTsignal=getFFTsignal(1:squareLength);
n=length(getFFTsignal);

minf=100; maxf=20000;

fft_result = fft(getFFTsignal, n);				% Compute FFT of the input signal.
fft_power = fft_result .* conj(fft_result);% / n;	% Compute power spectrum.  Dividing by 'n' we get the power spectral density.
fft_phase = angle(fft_result);			% Compute the phase spectrum.

frequencies = (1/dt)*(1:n/2)/n;
fft_power=fft_power(1:length(fft_power)/2); % remove mirror frequencies
fft_phase=fft_phase(1:length(fft_phase)/2); % remove mirror frequencies
fft_powerdB = UTIL_amp2dB (fft_power, max(fft_power)); % convert to dB
%     jags=find(diff(fft_phase)>0); % unwrap phase
%     for i=jags, fft_phase(i+1:end)=fft_phase(i+1:end)-2*pi; end

xlim([minf maxf])
semilogx(frequencies, fft_powerdB-max(fft_powerdB), color)
ylim([ -20 5])



function y=UTIL_lowpassFilterFreq(x, cutOffFrequency, dt)
% UTIL_lowpassFilterFreq multi-channel filter
%
% Usage:
% output=UTIL_lowpassFilterFreq(input, cutOffFrequency, dt)
%
% cutoff should be around 100 Hz for audio work
% dt should be <1/50000 s for audio work
%
% Attenuation 	       is - 6 dB per octave above cutoff.


sampleRate=1/dt;

if 4*cutOffFrequency>sampleRate
    warning(['UTIL_lowpassFilterFreq: sample rate ' num2str(1/dt) ' is too low for this BF.  Sampling rate should be >' num2str(4*cutOffFrequency) 'or cutoff (' num2str(4*cutOffFrequency) ') should be lower' ])
    cutOffFrequency=sampleRate/4;
end

tau=1/(2*pi*cutOffFrequency);

y=zeros(size(x));
[numChannels numPoints]=size(x);
for i=1:numChannels
    y(i,:)=filter([dt/tau], [1 -(1-dt/tau)], x(i,:));
end


