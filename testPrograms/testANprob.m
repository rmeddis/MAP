function testANprob(targetFrequency,BFlist, levels, ...
    paramsName, paramChanges)

% testANprob(1000,1000, -10:10:80, 'Normal')

global IHC_VResp_VivoParams  IHC_cilia_RPParams IHCpreSynapseParams
global AN_IHCsynapseParams
global ANprobRateOutput dt ANtauCas
global ARattenuation MOCattenuation

AN_spikesOrProbability='probability';

dbstop if error
addpath (['..' filesep 'MAP'], ['..' filesep 'utilities'], ...
    ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'],...
    ['..' filesep 'testPrograms'])

if nargin<5, paramChanges=[]; end
if nargin<4, paramsName='Normal'; end
if nargin<3, levels=-10:10:80; end
if nargin==0, targetFrequency=1000; BFlist=1000; end

nLevels=length(levels);

toneDuration=.2;
rampDuration=0.002;
silenceDuration=.02;
localPSTHbinwidth=0.001;

% Use only the first frequency in the GUI targetFrequency box to defineBF
% targetFrequency=stimulusParameters.targetFrequency(1);
% BFlist=targetFrequency;

AN_HSRonset=zeros(nLevels,1);
AN_HSRsaturated=zeros(nLevels,1);
AN_LSRonset=zeros(nLevels,1);
AN_LSRsaturated=zeros(nLevels,1);

AR=zeros(nLevels,1);
MOC=zeros(nLevels,1);

figure(15), clf
set(gcf,'position',[980   356   401   321])
drawnow

%% guarantee that the sample rate is at least 10 times the frequency
sampleRate=50000;
while sampleRate< 10* targetFrequency
    sampleRate=sampleRate+10000;
end

%% adjust sample rate so that the pure tone stimulus has an integer
%% numver of epochs in a period
dt=1/sampleRate;
period=1/targetFrequency;

%% main computational loop (vary level)
levelNo=0;
for leveldB=levels
    levelNo=levelNo+1;

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
    showPlotsAndDetails=0;


    MAP1_14(inputSignal, 1/dt, BFlist, ...
        paramsName, AN_spikesOrProbability, paramChanges);

    nTaus=length(ANtauCas);

    %LSR (same as HSR if no LSR fibers present)
    [nANFibers nTimePoints]=size(ANprobRateOutput);

    numLSRfibers=1;
    numHSRfibers=numLSRfibers;

    LSRspikes=ANprobRateOutput(1:numLSRfibers,:);
    PSTH=UTIL_PSTHmaker(LSRspikes, dt, localPSTHbinwidth);
    PSTHLSR=PSTH/(localPSTHbinwidth/dt);  % across fibers rates
    PSTHtime=localPSTHbinwidth:localPSTHbinwidth:...
        localPSTHbinwidth*length(PSTH);
    AN_LSRonset(levelNo)= max(PSTHLSR); % peak in 5 ms window
    AN_LSRsaturated(levelNo)= mean(PSTHLSR(round(length(PSTH)/2):end));

    % HSR
    HSRspikes= ANprobRateOutput(end- numHSRfibers+1:end, :);
    PSTH=UTIL_PSTHmaker(HSRspikes, dt, localPSTHbinwidth);
    PSTH=PSTH/(localPSTHbinwidth/dt); % sum across fibers (HSR only)
    AN_HSRonset(levelNo)= max(PSTH);
    AN_HSRsaturated(levelNo)= mean(PSTH(round(length(PSTH)/2): end));

    figure(15), subplot(2,2,4)
    hold off, bar(PSTHtime,PSTH, 'k')
    hold on,  bar(PSTHtime,PSTHLSR,'r')
    ylim([0 1000])
    xlim([0 length(PSTH)*localPSTHbinwidth])
    set(gcf,'name',[num2str(BFlist), ' Hz: ' num2str(leveldB) ' dB']);

    AR(levelNo)=min(ARattenuation);
    MOC(levelNo)=min(MOCattenuation(length(MOCattenuation)/2:end));


    figure(15), subplot(2,2,3)
    plot(20*log10(MOC), 'k'), hold on
    plot(20*log10(AR), 'r'),  hold off
    title(' MOC/AR'), ylabel('dB attenuation')
    ylim([-30 0])

end % level

figure(15), subplot(2,2,3)
plot(levels,20*log10(MOC), 'k'), hold on
plot(levels,20*log10(AR), 'r'),  hold off
title(' MOC/AR'), ylabel('dB attenuation')
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


nRows=2; nCols=2;

% AN rate - level ONSET functions
subplot(nRows,nCols,1)
plot(levels,AN_LSRonset,'ro'), hold on
plot(levels,AN_HSRonset,'ko'), hold off
ylim([0 1000]),  xlim([min(levels) max(levels)])
ttl=['tauCa= ' num2str(IHCpreSynapseParams.tauCa)];
title( ttl)
xlabel('level dB SPL'), ylabel('peak rate (sp/s)'), grid on
text(0, 800, 'AN onset', 'fontsize', 14)

% AN rate - level ADAPTED function
subplot(nRows,nCols,2)
plot(levels,AN_LSRsaturated, 'ro'), hold on
plot(levels,AN_HSRsaturated, 'ko'), hold off
maxYlim=340;
ylim([0 maxYlim])
set(gca,'ytick',0:50:300)
xlim([min(levels) max(levels)])
set(gca,'xtick',[levels(1):20:levels(end)])
%     grid on
ttl=[   'spont=' num2str(mean(AN_HSRsaturated(1,:)),'%4.0f')...
    '  sat=' num2str(mean(AN_HSRsaturated(end,1)),'%4.0f')];
title( ttl)
xlabel('level dB SPL'), ylabel ('adapted rate (sp/s)')
text(0, maxYlim-50, 'AN adapted', 'fontsize', 14), grid on

allData=[ levels'  AN_HSRonset AN_HSRsaturated...
    AN_LSRonset AN_LSRsaturated ];
fprintf('\n levels \tANHSR Onset \tANHSR adapted\tANLSR Onset \tANLSR adapted\tCNHSR\tCNLSR\tICHSR  \tICLSR \n');
UTIL_printTabTable(round(allData))


UTIL_showStruct(IHC_cilia_RPParams, 'IHC_cilia_RPParams')
UTIL_showStruct(IHCpreSynapseParams, 'IHCpreSynapseParams')
UTIL_showStruct(AN_IHCsynapseParams, 'AN_IHCsynapseParams')


allData=[ levels'  AN_HSRonset AN_HSRsaturated...
    AN_LSRonset AN_LSRsaturated ];
fprintf('\n levels \tANHSR Onset \tANHSR adapted\tANLSR Onset \tANLSR adapted\tCNHSR\tCNLSR\tICHSR  \tICLSR \n');
UTIL_printTabTable(round(allData))

