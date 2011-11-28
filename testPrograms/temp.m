function vectorStrength=testAN(probeFrequency,BFlist, levels, ...
    paramsName,paramChanges)
% testAN generates rate/level functions for AN and brainstem units.
%  also other information like PSTHs, MOC efferent activity levels.
% Vector strength calculations require the computation of period
% histograms. for this reason the sample rate must always be an integer
% multiple of the tone frequency. This applies to both the sampleRate and
% the spikes sampling rate.
% A full 'spikes' model is used.
% paramChanges is a cell array of strings containing MATLAB statements that
% change model parameters. See MAPparamsNormal for examples.
% e.g.
% testAN(1000,1000, -10:10:80,'Normal',{});


global IHC_VResp_VivoParams  IHC_cilia_RPParams IHCpreSynapseParams
global AN_IHCsynapseParams
global ANoutput dtSpikes CNoutput ICoutput ICmembraneOutput ANtauCas
global ARattenuation MOCattenuation
tic
dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'], ['..' filesep 'utilities'], ...
    ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'],...
    ['..' filesep 'testPrograms'])

if nargin<5, paramChanges={'IHCpreSynapseParams.tauCa= [30e-6 90e-6];'}; end
if nargin<4, paramsName='PL'; end
if nargin<3, levels=-10:10:80; end
if nargin==0, probeFrequency=1000; BFlist=1000; end
nLevels=length(levels);

toneDuration=.2;   rampDuration=0.005;   silenceDuration=.02;
localPSTHbinwidth=0.001;

%% guarantee that the sample rate is an interger multiple 
%   of the probe frequency
% we want 5 bins per period for spikes
spikesSampleRate=5*probeFrequency;
% model sample rate must be an integer multiple of this and in the region
% of 50000
sampleRate=spikesSampleRate*round(50000/spikesSampleRate);
dt=1/sampleRate;
% avoid very slow spikes sampling rate
spikesSampleRate=spikesSampleRate*ceil(10000/spikesSampleRate);

%% pre-allocate storage
AN_HSRonset=zeros(nLevels,1);
AN_HSRsaturated=zeros(nLevels,1);
AN_LSRonset=zeros(nLevels,1);
AN_LSRsaturated=zeros(nLevels,1);
CNLSRrate=zeros(nLevels,1);
CNHSRsaturated=zeros(nLevels,1);
ICHSRsaturated=zeros(nLevels,1);
ICLSRsaturated=zeros(nLevels,1);
vectorStrength=zeros(nLevels,1);

AR=zeros(nLevels,1);
MOC=zeros(nLevels,1);

figure(15), clf
set(gcf,'position',[980   356   401   321])
figure(5), clf
set(gcf,'position', [980 34 400 295])
drawnow


%% main computational loop (vary level)
levelNo=0;
for leveldB=levels
    levelNo=levelNo+1;
    amp=28e-6*10^(leveldB/20);
    fprintf('%4.0f\t', leveldB)

    %% generate tone and silences
    time=dt:dt:toneDuration;
    rampTime=dt:dt:rampDuration;
    ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
        ones(1,length(time)-length(rampTime))];
    ramp=ramp.*fliplr(ramp);
    
    silence=zeros(1,round(silenceDuration/dt));
    
    inputSignal=amp*sin(2*pi*probeFrequency'*time);
    inputSignal= ramp.*inputSignal;
    inputSignal=[silence inputSignal];
    
    %% run the model
    AN_spikesOrProbability='spikes';  
    nExistingParamChanges=length(paramChanges);
    paramChanges{nExistingParamChanges+1}=...
        ['AN_IHCsynapseParams.spikesTargetSampleRate=' ...
        num2str(spikesSampleRate) ';'];

    MAP1_14(inputSignal, 1/dt, BFlist, ...
        paramsName, AN_spikesOrProbability, paramChanges);
    
    nTaus=length(ANtauCas);
    
    %% Auditory nerve evaluate and display (Fig. 5)
    %LSR (same as HSR if no LSR fibers present)
    [nANFibers nTimePoints]=size(ANoutput);
    numLSRfibers=nANFibers/nTaus;
    numHSRfibers=numLSRfibers;
    
    LSRspikes=ANoutput(1:numLSRfibers,:);
    PSTH=UTIL_PSTHmaker(LSRspikes, dtSpikes, localPSTHbinwidth);
    PSTHLSR=mean(PSTH,1)/localPSTHbinwidth;  % across fibers rates
    PSTHtime=localPSTHbinwidth:localPSTHbinwidth:...
        localPSTHbinwidth*length(PSTH);
    AN_LSRonset(levelNo)= max(PSTHLSR); % peak in 5 ms window
    AN_LSRsaturated(levelNo)= mean(PSTHLSR(round(length(PSTH)/2):end));
    
    % AN HSR
    HSRspikes= ANoutput(end- numHSRfibers+1:end, :);
    PSTH=UTIL_PSTHmaker(HSRspikes, dtSpikes, localPSTHbinwidth);
    PSTH=mean(PSTH,1)/localPSTHbinwidth; % sum across fibers (HSR only)
    AN_HSRonset(levelNo)= max(PSTH);
    AN_HSRsaturated(levelNo)= mean(PSTH(round(length(PSTH)/2): end));
    
    figure(5), subplot(2,2,2)
    hold off, bar(PSTHtime,PSTH, 'b')
    hold on,  bar(PSTHtime,PSTHLSR,'r')
    ylim([0 1000])
    xlim([0 length(PSTH)*localPSTHbinwidth])
    grid on
    set(gcf,'name',['PSTH: ' num2str(BFlist), ' Hz: ' num2str(leveldB) ' dB']);
    
    % AN - CV
    %  CV is computed 5 times. Use the middle one (3) as most typical
    cvANHSR=  UTIL_CV(HSRspikes, dtSpikes);
    
    % AN - vector strength
    PSTH=sum(CNoutput (11:20,:));
    [PH, binTimes]=UTIL_periodHistogram...
        (PSTH, dtSpikes, probeFrequency);
    VS=UTIL_vectorStrength(PH);
    vectorStrength(levelNo)=VS;
    disp(['sat rate= ' num2str(AN_HSRsaturated(levelNo)) ...
        ';   phase-locking VS = ' num2str(VS)])
    title(['AN HSR: CV=' num2str(cvANHSR(3),'%5.2f') ...
        'VS=' num2str(VS,'%5.2f')])
    
    %% CN - first-order neurons
    
    % CN driven by LSR fibers
    [nCNneurons c]=size(CNoutput);
    nLSRneurons=round(nCNneurons/nTaus);
    CNLSRspikes=CNoutput(1:nLSRneurons,:);
    PSTH=UTIL_PSTHmaker(CNLSRspikes, dtSpikes, localPSTHbinwidth);
    PSTH=sum(PSTH)/nLSRneurons;
%     CNLSRrate(levelNo)=mean(PSTH(round(length(PSTH)/2):end))/localPSTHbinwidth;
    CNLSRrate(levelNo)=mean(PSTH)/localPSTHbinwidth;
    
    %CN HSR
    MacGregorMultiHSRspikes=...
        CNoutput(end-nLSRneurons+1:end,:);
    PSTH=UTIL_PSTHmaker(MacGregorMultiHSRspikes, dtSpikes, localPSTHbinwidth);
    PSTH=sum(PSTH)/nLSRneurons;
    PSTH=mean(PSTH,1)/localPSTHbinwidth; % sum across fibers (HSR only)
    
%     CNHSRsaturated(levelNo)=mean(PSTH(length(PSTH)/2:end));
    CNHSRsaturated(levelNo)=mean(PSTH);
    
    figure(5), subplot(2,2,3)
    bar(PSTHtime,PSTH)
    ylim([0 1000])
    xlim([0 length(PSTH)*localPSTHbinwidth])
    cvMMHSR= UTIL_CV(MacGregorMultiHSRspikes, dtSpikes);
    title(['CN    CV= ' num2str(cvMMHSR(3),'%5.2f')])
    
    %% IC LSR
    [nICneurons c]=size(ICoutput);
    nLSRneurons=round(nICneurons/nTaus);
    ICLSRspikes=ICoutput(1:nLSRneurons,:);
    PSTH=UTIL_PSTHmaker(ICLSRspikes, dtSpikes, localPSTHbinwidth);
%     ICLSRsaturated(levelNo)=mean(PSTH(round(length(PSTH)/2):end))/localPSTHbinwidth;
    ICLSRsaturated(levelNo)=mean(PSTH)/localPSTHbinwidth;
    
    %IC HSR
    MacGregorMultiHSRspikes=...
        ICoutput(end-nLSRneurons+1:end,:);
    PSTH=UTIL_PSTHmaker(MacGregorMultiHSRspikes, dtSpikes, localPSTHbinwidth);
    ICHSRsaturated(levelNo)= (sum(PSTH)/nLSRneurons)/toneDuration;

    AR(levelNo)=min(ARattenuation);
    MOC(levelNo)=min(MOCattenuation(length(MOCattenuation)/2:end));
    
    time=dt:dt:dt*size(ICmembraneOutput,2);
    figure(5), subplot(2,2,4)
    % plot HSR (last of two)
    plot(time,ICmembraneOutput(end-nLSRneurons+1, 1:end),'k')
    ylim([-0.07 0])
    xlim([0 max(time)])
    title(['IC  ' num2str(leveldB,'%4.0f') 'dB'])
    drawnow
    
    figure(5), subplot(2,2,1)
    plot(20*log10(MOC), 'k'), hold on
    plot(20*log10(AR), 'r'), hold off
    title(' MOC'), ylabel('dB attenuation')
    ylim([-30 0])
    
%     x=input('prompt')
end % level

%% plot with levels on x-axis
figure(5), subplot(2,2,1)
plot(levels,20*log10(MOC), 'k'),hold on
plot(levels,20*log10(AR), 'r'), hold off
title(' MOC'), ylabel('dB attenuation')
ylim([-30 0])
xlim([0 max(levels)])

fprintf('\n')
toneDuration=2;
rampDuration=0.004;
silenceDuration=.02;
nRepeats=200;   % no. of AN fibers
fprintf('toneDuration  %6.3f\n', toneDuration)
fprintf(' %6.0f  AN fibers (repeats)\n', nRepeats)
fprintf('levels')
fprintf('%6.2f\t', levels)
fprintf('\n')

%% ---------------------- display summary results (Fig 15)
figure(15), clf
nRows=2; nCols=2;

% AN rate - level ONSET functions
subplot(nRows,nCols,1)
plot(levels,AN_LSRonset,'ro'), hold on
plot(levels,AN_HSRonset,'ko'), hold off
ylim([0 1000])
if length(levels)>1
    xlim([min(levels) max(levels)])
end

ttl=['tauCa= ' num2str(IHCpreSynapseParams.tauCa)];
title( ttl)
xlabel('level dB SPL'), ylabel('peak rate (sp/s)'), grid on
text(0, 800, 'AN onset', 'fontsize', 14)

% AN rate - level ADAPTED function
subplot(nRows,nCols,2)
plot(levels,AN_LSRsaturated, 'ro'), hold on
plot(levels,AN_HSRsaturated, 'ko'), hold off
ylim([0 400])
set(gca,'ytick',0:50:300)
if length(levels)>1
xlim([min(levels) max(levels)])
end
set(gca,'xtick',[levels(1):20:levels(end)])
%     grid on
ttl=[   'spont=' num2str(mean(AN_HSRsaturated(1,:)),'%4.0f')...
    '  sat=' num2str(mean(AN_HSRsaturated(end,1)),'%4.0f')];
title( ttl)
xlabel('level dB SPL'), ylabel ('adapted rate (sp/s)')
text(0, 340, 'AN adapted', 'fontsize', 14), grid on

% CN rate - level function
subplot(nRows,nCols,3)
plot(levels,CNLSRrate, 'ro'), hold on
plot(levels,CNHSRsaturated, 'ko'), hold off
ylim([0 400])
set(gca,'ytick',0:50:300)
if length(levels)>1
xlim([min(levels) max(levels)])
end
set(gca,'xtick',[levels(1):20:levels(end)])
%     grid on
ttl=[   'spont=' num2str(mean(CNHSRsaturated(1,:)),'%4.0f') '  sat=' ...
    num2str(mean(CNHSRsaturated(end,1)),'%4.0f')];
title( ttl)
xlabel('level dB SPL'), ylabel ('adapted rate (sp/s)')
text(0, 350, 'CN', 'fontsize', 14), grid on

% IC rate - level function
subplot(nRows,nCols,4)
plot(levels,ICLSRsaturated, 'ro'), hold on
plot(levels,ICHSRsaturated, 'ko'), hold off
ylim([0 400])
set(gca,'ytick',0:50:300)
if length(levels)>1
xlim([min(levels) max(levels)])
end
set(gca,'xtick',[levels(1):20:levels(end)]), grid on
ttl=['spont=' num2str(mean(ICHSRsaturated(1,:)),'%4.0f') ...
    '  sat=' num2str(mean(ICHSRsaturated(end,1)),'%4.0f')];
title( ttl)
xlabel('level dB SPL'), ylabel ('adapted rate (sp/s)')
text(0, 350, 'IC', 'fontsize', 14)
set(gcf,'name',' AN CN IC rate/level')

fprintf('\n')
disp('levels vectorStrength')
fprintf('%3.0f \t %6.4f \n', [levels; vectorStrength'])
fprintf('\n')
fprintf('Phase locking, max vector strength=\t %6.4f\n\n',...
    max(vectorStrength))

allData=[ levels'  AN_HSRonset AN_HSRsaturated...
    AN_LSRonset AN_LSRsaturated ...
    CNHSRsaturated CNLSRrate...
    ICHSRsaturated ICLSRsaturated];
fprintf('\n levels \tANHSR Onset \tANHSR adapted\tANLSR Onset \tANLSR adapted\tCNHSR\tCNLSR\tICHSR  \tICLSR \n');
UTIL_printTabTable(round(allData))
fprintf('VS (phase locking)= \t%6.4f\n\n',...
    max(vectorStrength))

UTIL_showStruct(IHC_cilia_RPParams, 'IHC_cilia_RPParams')
UTIL_showStruct(IHCpreSynapseParams, 'IHCpreSynapseParams')
UTIL_showStruct(AN_IHCsynapseParams, 'AN_IHCsynapseParams')

fprintf('\n')
disp('levels vectorStrength')
fprintf('%3.0f \t %6.4f \n', [levels; vectorStrength'])
fprintf('\n')
fprintf('Phase locking, max vector strength= \t%6.4f\n\n',...
    max(vectorStrength))

allData=[ levels'  AN_HSRonset AN_HSRsaturated...
    AN_LSRonset AN_LSRsaturated ...
    CNHSRsaturated CNLSRrate...
    ICHSRsaturated ICLSRsaturated];
fprintf('\n levels \tANHSR Onset \tANHSR adapted\tANLSR Onset \tANLSR adapted\tCNHSR\tCNLSR\tICHSR  \tICLSR \n');
UTIL_printTabTable(round(allData))
fprintf('VS (phase locking)= \t%6.4f\n\n',...
    max(vectorStrength))

path(restorePath)
toc
