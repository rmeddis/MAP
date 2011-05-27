function testAN
% testIHC used either for IHC I/O function ...
%  or receptive field (doReceptiveFields=1)

global experiment method stimulusParameters
global IHC_VResp_VivoParams IHCpreSynapseParams
global AN_IHCsynapseParams
% global saveMembranePotential MacGregorMultiParams
dbstop if error

addpath (['..' filesep 'MAP'], ['..' filesep 'utilities'], ...
    ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'],...
    ['..' filesep 'testPrograms'])

levels=-10:10:80;
% levels=80:10:90;
nLevels=length(levels);

toneDuration=.2;
rampDuration=0.002;
silenceDuration=.02;
localPSTHbinwidth=0.001;

% Use only the first frequency in the GUI targetFrequency box to defineBF
targetFrequency=stimulusParameters.targetFrequency(1);
BFlist=targetFrequency;

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

% ANoutput=zeros(200,200);

figure(15), clf
set(gcf,'position',[980   356   401   321])
figure(5), clf
 set(gcf,'position', [980 34 400 295])
drawnow

levelNo=0;
for leveldB=levels
    levelNo=levelNo+1;

    % sample rate should be amultiple of the targetFrequency for PSTH below
    sampleRate=50000;
    dt=1/sampleRate;
    period=1/targetFrequency;
    dt=dt*(dt/period)*round(period/dt);

    fprintf('%4.0f\t', leveldB)
    amp=28e-6*10^(leveldB/20);

    time=dt:dt:toneDuration;
    rampTime=dt:dt:rampDuration;
    ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
        ones(1,length(time)-length(rampTime))];
    ramp=ramp.*fliplr(ramp);

    silence=zeros(1,round(silenceDuration/dt));

    % create signal (leveldB/ targetFrequency)
    inputSignal=amp*sin(2*pi*targetFrequency'*time);
    inputSignal= ramp.*inputSignal;
    inputSignal=[silence inputSignal];

    %% run the model
    AN_spikesOrProbability='spikes';
    MAPparamsName=experiment.name;
    showPlotsAndDetails=0;

    global ANoutput CNoutput ICoutput ICmembraneOutput tauCas
    global ARattenuation MOCattenuation
    
    MAP1_14(inputSignal, 1/dt, targetFrequency, ...
        MAPparamsName, AN_spikesOrProbability);

    nTaus=length(tauCas);

    %LSR (same as HSR if no LSR fibers present)
    [nANFibers nTimePoints]=size(ANoutput);
    dt=dt* length(inputSignal)/nTimePoints;

    numLSRfibers=nANFibers/nTaus;
    numHSRfibers=numLSRfibers;

    LSRspikes=ANoutput(1:numLSRfibers,:);
    PSTH=UTIL_makePSTH(LSRspikes, dt, localPSTHbinwidth);
    PSTHLSR=mean(PSTH,1)/localPSTHbinwidth;  % across fibers rates
    PSTHtime=localPSTHbinwidth:localPSTHbinwidth:...
        localPSTHbinwidth*length(PSTH);
    AN_LSRonset(levelNo)= max(PSTHLSR); % peak in 5 ms window
    AN_LSRsaturated(levelNo)= mean(PSTHLSR(round(length(PSTH)/2):end));

    % HSR
    HSRspikes= ANoutput(end- numHSRfibers+1:end, :);
    PSTH=UTIL_makePSTH(HSRspikes, dt, localPSTHbinwidth);
    PSTH=mean(PSTH,1)/localPSTHbinwidth; % sum across fibers (HSR only)
    AN_HSRonset(levelNo)= max(PSTH);
    AN_HSRsaturated(levelNo)= mean(PSTH(round(length(PSTH)/2): end));

    figure(5), subplot(2,2,2)
    hold off, bar(PSTHtime,PSTH, 'b')
    hold on,  bar(PSTHtime,PSTHLSR,'r')
    ylim([0 1000])
xlim([0 length(PSTH)*localPSTHbinwidth])
    
    % AN - CV
    %  CV is computed 5 times. Use the middle one (3) as most typical
    cvANHSR=  UTIL_CV(HSRspikes, dt);

    % AN - vector strength
    PSTH=sum(HSRspikes);
    [PH, binTimes]=UTIL_periodHistogram...
        (PSTH, dt, targetFrequency);
    VS=UTIL_vectorStrength(PH);
    vectorStrength(levelNo)=VS;
    disp(['sat rate= ' num2str(AN_HSRsaturated(levelNo)) ...
        ';   phase-locking VS = ' num2str(VS)])
    title(['AN HSR: CV=' num2str(cvANHSR(3),'%5.2f') ...
        'VS=' num2str(VS,'%5.2f')])

    % CN - first-order neurons

    % CN LSR
    [nCNneurons c]=size(CNoutput);
    nLSRneurons=round(nCNneurons/nTaus);
    CNLSRspikes=CNoutput(1:nLSRneurons,:);
    PSTH=UTIL_makePSTH(CNLSRspikes, dt, localPSTHbinwidth);
    PSTH=sum(PSTH)/nLSRneurons;
    CNLSRrate(levelNo)=mean(PSTH(round(length(PSTH)/2):end))/localPSTHbinwidth;

    %CN HSR
    MacGregorMultiHSRspikes=...
        CNoutput(end-nLSRneurons:end,:);
    PSTH=UTIL_makePSTH(MacGregorMultiHSRspikes, dt, localPSTHbinwidth);
    PSTH=sum(PSTH)/nLSRneurons;
    PSTH=mean(PSTH,1)/localPSTHbinwidth; % sum across fibers (HSR only)

    CNHSRsaturated(levelNo)=mean(PSTH(length(PSTH)/2:end));

    figure(5), subplot(2,2,3)
    bar(PSTHtime,PSTH)
    ylim([0 1000])
    xlim([0 length(PSTH)*localPSTHbinwidth])
    cvMMHSR= UTIL_CV(MacGregorMultiHSRspikes, dt);
    title(['CN    CV= ' num2str(cvMMHSR(3),'%5.2f')])

    % IC LSR
    [nICneurons c]=size(ICoutput);
    nLSRneurons=round(nICneurons/nTaus);
    ICLSRspikes=ICoutput(1:nLSRneurons,:);
    PSTH=UTIL_makePSTH(ICLSRspikes, dt, localPSTHbinwidth);
    ICLSRsaturated(levelNo)=mean(PSTH(round(length(PSTH)/2):end))/localPSTHbinwidth;

    %IC HSR
    MacGregorMultiHSRspikes=...
        ICoutput(end-nLSRneurons:end,:);
    PSTH=UTIL_makePSTH(MacGregorMultiHSRspikes, dt, localPSTHbinwidth);
    PSTH=sum(PSTH)/nLSRneurons;
    PSTH=mean(PSTH,1)/localPSTHbinwidth; % sum across fibers (HSR only)

    ICHSRsaturated(levelNo)=mean(PSTH(length(PSTH)/2:end));

    AR(levelNo)=min(ARattenuation);
    MOC(levelNo)=min(MOCattenuation(length(MOCattenuation)/2:end));

    time=dt:dt:dt*size(ICmembraneOutput,2);
    figure(5), subplot(2,2,4)
    plot(time,ICmembraneOutput(2, 1:end),'k')
    ylim([-0.07 0])
    xlim([0 max(time)])
    title(['IC  ' num2str(leveldB,'%4.0f') 'dB'])
    drawnow
    
    figure(5), subplot(2,2,1)
    plot(20*log10(MOC), 'k'),
    title(' MOC'), ylabel('dB attenuation')
    ylim([-30 0])


end % level
figure(5), subplot(2,2,1)
    plot(levels,20*log10(MOC), 'k'),
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


% ---------------------------------------------------- display parameters

figure(15), clf
set(gcf,'position',[1000   356   381   321])
nRows=2; nCols=2;

% AN rate - level ONSET functions
subplot(nRows,nCols,1)
plot(levels,AN_LSRonset,'ro'), hold on
plot(levels,AN_HSRonset,'ko'), hold off
ylim([0 1000]),  xlim([min(levels) max(levels)])
ttl=['tauCa= ' num2str(IHCpreSynapseParams.tauCa)];
title( ttl)
xlabel('level dB SPL'), ylabel('peak rate (sp/s)'), grid on
text(0, 800, 'AN onset', 'fontsize', 16)

% AN rate - level ADAPTED function
subplot(nRows,nCols,2)
plot(levels,AN_LSRsaturated, 'ro'), hold on
plot(levels,AN_HSRsaturated, 'ko'), hold off
ylim([0 400])
set(gca,'ytick',0:50:300)
xlim([min(levels) max(levels)])
set(gca,'xtick',[levels(1):20:levels(end)])
%     grid on
ttl=[   'spont=' num2str(mean(AN_HSRsaturated(1,:)),'%4.0f')...
    '  sat=' num2str(mean(AN_HSRsaturated(end,1)),'%4.0f')];
title( ttl)
xlabel('level dB SPL'), ylabel ('adapted rate (sp/s)')
text(0, 340, 'AN adapted', 'fontsize', 16), grid on

% CN rate - level ADAPTED function
subplot(nRows,nCols,3)
plot(levels,CNLSRrate, 'ro'), hold on
plot(levels,CNHSRsaturated, 'ko'), hold off
ylim([0 400])
set(gca,'ytick',0:50:300)
xlim([min(levels) max(levels)])
set(gca,'xtick',[levels(1):20:levels(end)])
%     grid on
ttl=[   'spont=' num2str(mean(CNHSRsaturated(1,:)),'%4.0f') '  sat=' ...
    num2str(mean(CNHSRsaturated(end,1)),'%4.0f')];
title( ttl)
xlabel('level dB SPL'), ylabel ('adapted rate (sp/s)')
text(0, 350, 'CN', 'fontsize', 16), grid on

% IC rate - level ADAPTED function
subplot(nRows,nCols,4)
plot(levels,ICLSRsaturated, 'ro'), hold on
plot(levels,ICHSRsaturated, 'ko'), hold off
ylim([0 400])
set(gca,'ytick',0:50:300)
xlim([min(levels) max(levels)])
set(gca,'xtick',[levels(1):20:levels(end)]), grid on

ttl=['spont=' num2str(mean(ICHSRsaturated(1,:)),'%4.0f') ...
    '  sat=' num2str(mean(ICHSRsaturated(end,1)),'%4.0f')];
title( ttl)
xlabel('level dB SPL'), ylabel ('adapted rate (sp/s)')
text(0, 350, 'IC', 'fontsize', 16)
set(gcf,'name',' AN CN IC rate/level')

UTIL_showStruct(IHCpreSynapseParams, 'IHCpreSynapseParams')
UTIL_showStruct(AN_IHCsynapseParams, 'AN_IHCsynapseParams')

fprintf('\n')
disp('levels vectorStrength')
fprintf('%3.0f \t %6.4f \n', [levels; vectorStrength'])
fprintf('\n')
fprintf('Phase locking, max vector strength= %6.4f\n\n',...
    max(vectorStrength))

allData=[ levels'  AN_HSRonset AN_HSRsaturated...
    AN_LSRonset AN_LSRsaturated ...
    CNHSRsaturated CNLSRrate...
    ICHSRsaturated ICLSRsaturated];
fprintf('\n levels \tANHSR Onset \tANHSR adapted\tANLSR Onset \tANLSR adapted\tCNHSR\tCNLSR\tICHSR  \tICLSR \n');
UTIL_printTabTable(round(allData))
fprintf('VS (phase locking)= \t%6.4f\n\n',...
    max(vectorStrength))

