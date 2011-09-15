function testFM(BFlist,paramsName, AN_spikesOrProbability,...
    paramChanges)
% testFM(1000,'Normal','probability', [])
% testFM(1000,'Normal','spikes', [])

% Demonstration is based on Harris and Dallos
%   specify whether you want AN 'probability' or 'spikes'
%       spikes is more realistic but takes longer
%       refractory effect is included only for spikes
%


global inputStimulusParams outerMiddleEarParams DRNLParams
global IHC_VResp_VivoParams IHCpreSynapseParams  AN_IHCsynapseParams
global  ANprobRateOutput  ANoutput ANtauCas  ANdt
dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'], ['..' filesep 'utilities'], ...
    ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'],...
    ['..' filesep 'testPrograms'])

if nargin==0
    BFlist=1000;
    paramsName=('Normal');
    AN_spikesOrProbability='spikes';
    paramChanges=[];
else
    if nargin<3
    paramChanges=[];
    end
end

% masker and probe levels are relative to this threshold
thresholdAtCF=10; % dB SPL
maskerLevels=[-80   10 20 30 40 60 ] + thresholdAtCF;

showPSTHs=1;

sampleRate=50000;
dt=1/sampleRate;

% fetch BF from GUI: use only the first target frequency
maskerFrequency=BFlist;
maskerDuration=.1;

targetFrequency=maskerFrequency;
probeLeveldB=20+thresholdAtCF;	% H&D use 20 dB SL/ TMC uses 10 dB SL
probeDuration=0.015; % HD use 15 ms probe (fig 3).

rampDuration=.001;  % HD use 1 ms linear ramp
initialSilenceDuration=0.02;
finalSilenceDuration=0.05;      % useful for showing recovery

figure(7), clf
set(gcf,'position',[613    36   360   247])
set(gcf,'name', ['forward masking: thresholdAtCF=' num2str(thresholdAtCF)])

if showPSTHs
    figure(8), clf
    set(gcf,'name', 'Harris and Dallos simulation')
    set(gcf,'position',[980    36   380   249])
end

% Plot Harris and Dallos result for comparison
figure(7)
gapDurations=[0.001	0.002	0.005	0.01	0.02	0.05	0.1	0.3];
HDmaskerLevels=[+10	+20	+30	+40	+60];
HDresponse=[
    1 1 1 1 1 1 1 1;
    0.8  	0.82	0.81	0.83	0.87	0.95	1	    1;
    0.48	0.5	    0.54	0.58	0.7	    0.85	0.95	1;
    0.3	    0.31	0.35	0.4	    0.5	    0.68	0.82	0.94;
    0.2 	0.27	0.27	0.29	0.42	0.64	0.75	0.92;
    0.15	0.17	0.18	0.23	0.3	     0.5	0.6	    0.82];
semilogx(gapDurations,HDresponse,'o'), hold on
legend(strvcat(num2str(maskerLevels')),'location','southeast')
title([ 'masker dB: ' num2str(HDmaskerLevels)])

%% Run the trials
maxProbeResponse=0;
levelNo=0;
resultsMatrix=zeros(length(maskerLevels), length(gapDurations));
summaryFiringRates=[];

% initial silence
time=dt: dt: initialSilenceDuration;
initialSilence=zeros(1,length(time));

% probe
time=dt: dt: probeDuration;
amp=28e-6*10^(probeLeveldB/20);
probe=amp*sin(2*pi.*targetFrequency*time);
% ramps
% catch rampTime error
if rampDuration>0.5*probeDuration, rampDuration=probeDuration/2; end
rampTime=dt:dt:rampDuration;
% raised cosine ramp
ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
    ones(1,length(time)-length(rampTime))];
%  onset ramp
probe=probe.*ramp;
%  offset ramp makes complete ramp for probe
ramp=fliplr(ramp);
% apply ramp mask to probe. Probe does not change below
probe=probe.*ramp;

% final silence
time=dt: dt: finalSilenceDuration;
finalSilence=zeros(1,length(time));

PSTHplotCount=0;
longestSignalDuration=initialSilenceDuration + maskerDuration +...
    max(gapDurations) + probeDuration + finalSilenceDuration ;
for maskerLeveldB=maskerLevels
    levelNo=levelNo+1;
    allDurations=[];
    allFiringRates=[];

    %masker
    time=dt: dt: maskerDuration;
    masker=sin(2*pi.*maskerFrequency*time);
    % masker ramps
    if rampDuration>0.5*maskerDuration
        % catch ramp duration error
        rampDuration=maskerDuration/2;
    end
    % NB masker ramp (not probe ramp)
    rampTime=dt:dt:rampDuration;
    % raised cosine ramp
    ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi))...
        ones(1,length(time)-length(rampTime))];
    %  onset ramp
    masker=masker.*ramp;
    %  offset ramp
    ramp=fliplr(ramp);
    % apply ramp
    masker=masker.*ramp;
    amp=28e-6*10^(maskerLeveldB/20);
    maskerPa=amp*masker;

    for gapDuration=gapDurations
        time=dt: dt: gapDuration;
        gap=zeros(1,length(time));

        inputSignal=...
            [initialSilence maskerPa gap probe finalSilence];
%         time=dt:dt:length(inputSignal)*dt;
%         figure(99), plot(time,inputSignal)

        % **********************************  run MAP model
%         showPlotsAndDetails=0;
nChanges=length(paramChanges);
paramChanges{nChanges+1}='AN_IHCsynapseParams.numFibers=	500;'; 
        MAP1_14(inputSignal, 1/dt, targetFrequency, ...
            paramsName, AN_spikesOrProbability, paramChanges);

        if strcmp(AN_spikesOrProbability,'probability')
            [nFibers c]=size(ANprobRateOutput);
            nLSRfibers=nFibers/length(ANtauCas);
            ANresponse=ANprobRateOutput(end-nLSRfibers:end,:);
            ANdt=dt; % no adjustment for spikes speedup
        else
            [nFibers c]=size(ANoutput);
            nLSRfibers=nFibers/length(ANtauCas);
            ANresponse=ANoutput(end-nLSRfibers:end,:);
        end

        ANresponse=sum(ANresponse);
%         ANresponseTimes=ANdt:ANdt:length(ANresponse)*ANdt;
%         figure(99), plot(ANresponseTimes,ANresponse)

        % analyse results
        probeStart=initialSilenceDuration+maskerDuration+gapDuration;
        PSTHbinWidth=ANdt;
        responseDelay=0.005;
        probeTimes=probeStart+responseDelay:...
            PSTHbinWidth:probeStart+probeDuration+responseDelay;
        probeIDX=round(probeTimes/PSTHbinWidth);
        probePSTH=ANresponse(probeIDX);
        firingRate=mean(probePSTH);
        % NB this only works if you start with the lowest level masker
        maxProbeResponse=max([maxProbeResponse firingRate]);
        allDurations=[allDurations gapDuration];
        allFiringRates=[allFiringRates firingRate];

        %% show PSTHs
        if showPSTHs
            nLevels=length(maskerLevels);
            nDurations=length(gapDurations);
            figure(8)
            PSTHbinWidth=1e-3;
            PSTH=UTIL_PSTHmaker(ANresponse, ANdt, PSTHbinWidth);
            PSTH=PSTH*ANdt/PSTHbinWidth;
            PSTHplotCount=PSTHplotCount+1;
            subplot(nLevels,nDurations,PSTHplotCount)
            PSTHtime=PSTHbinWidth:PSTHbinWidth:...
                PSTHbinWidth*length(PSTH);
            if strcmp(AN_spikesOrProbability, 'spikes')
            bar(PSTHtime, PSTH/PSTHbinWidth/nFibers)
%                 ylim([0 500])
            else
            bar(PSTHtime, PSTH)
                ylim([0 500])
            end
%             xlim([0 longestSignalDuration])
            grid on

            if PSTHplotCount< (nLevels-1) * nDurations+1
                set(gca,'xticklabel',[])
            end

            if ~isequal(mod(PSTHplotCount,nDurations),1)
                set(gca,'yticklabel',[])
            else
                ylabel([num2str(maskerLevels...
                    (round(PSTHplotCount/nDurations) +1))])
            end

            if PSTHplotCount<=nDurations
                title([num2str(1000*gapDurations(PSTHplotCount)) 'ms'])
            end

%         figure(99),            bar(PSTHtime, PSTH)

        end % showPSTHs

    end     % gapDurations duration
    summaryFiringRates=[summaryFiringRates allFiringRates'];

    figure(7), hold on
    semilogx(allDurations, allFiringRates/maxProbeResponse)
    ylim([0 1]), hold on
    ylim([0 inf]), xlim([min(gapDurations) max(gapDurations)])
    xlim([0.001 1])
    pause(0.1) % to allow for CTRL/C interrupts
    resultsMatrix(levelNo,:)=allFiringRates;
end          % maskerLevel

disp('delay/ rates')
disp(num2str(round( [1000*allDurations' summaryFiringRates])))

% replot with best adjustment
% figure(34), clf% use for separate plot
figure(7), clf
peakProbe=max(max(resultsMatrix));
resultsMatrix=resultsMatrix/peakProbe;
semilogx(gapDurations,HDresponse,'o'), hold on
title(['FrMsk: probe ' num2str(probeLeveldB)...
    'dB SL: peakProbe=' num2str(peakProbe,'%5.0f') ' sp/s'])
xlabel('gap duration (s)'), ylabel ('probe response')
semilogx(allDurations, resultsMatrix'), ylim([0 1])
ylim([0 inf]),
xlim([0.001 1])
legend(strvcat(num2str(maskerLevels'-thresholdAtCF)), -1)

% ------------------------------------------------- display parameters
disp(['parameter file was: ' paramsName])
fprintf('\n')
% UTIL_showStruct(inputStimulusParams, 'inputStimulusParams')
% UTIL_showStruct(outerMiddleEarParams, 'outerMiddleEarParams')
% UTIL_showStruct(DRNLParams, 'DRNLParams')
% UTIL_showStruct(IHC_VResp_VivoParams, 'IHC_VResp_VivoParams')
UTIL_showStruct(IHCpreSynapseParams, 'IHCpreSynapseParams')
UTIL_showStruct(AN_IHCsynapseParams, 'AN_IHCsynapseParams')


path(restorePath);
