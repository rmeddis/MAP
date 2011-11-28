function testRF
% testIHC used either for IHC I/O function or receptive field (doReceptiveFields=1)

global experiment method stimulusParameters expGUIhandles
global inputStimulusParams  IHC_ciliaParams
global IHC_VResp_VivoParams IHCpreSynapseParams  AN_IHCsynapseParams
dbstop if error
% set(expGUIhandles.pushbuttonStop, 'backgroundColor', [.941 .941 .941])

addpath (['..' filesep 'MAP'], ['..' filesep 'utilities'], ...
    ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'],...
    ['..' filesep 'testPrograms'])

targetFrequency=stimulusParameters.targetFrequency(1);

sampleRate=50000;
doReceptiveFields=1;

toneDuration=.05;
rampDuration=0.004;
silenceDuration=.02;

nRepeats=100;   % no. of AN fibers

plotGraphsForIHC=1;
% number of MacGregor units is set in the parameter file.

if doReceptiveFields
    % show all receptive field
    frequencies=targetFrequency*    [  0.5         0.7         0.9     1       1.1         1.3         1.6];
    levels=0:20:80; nLevels=length(levels);
    figure(14), clf
    figure(15), clf
else
    % show only I/O function at BF
    frequencies=targetFrequency;
    levels=-20:10:90;
    %         levels=10:.25:13;
    %         levels=-20:1:-15
    nLevels=length(levels);
%     figure(13), clf,
%     set (gcf, 'name', ['IHC/AN  input/output' num2str(AN_IHCsynapseParams.numFibers) ' repeats'])
%     drawnow
end
nFrequencies=length(frequencies);

IHC_RP_peak=zeros(nLevels,nFrequencies);
IHC_RP_min=zeros(nLevels,nFrequencies);
IHC_RP_dc=zeros(nLevels,nFrequencies);
AN_HSRonset=zeros(nLevels,nFrequencies);
AN_HSRsaturated=zeros(nLevels,nFrequencies);
AN_LSRonset=zeros(nLevels,nFrequencies);
AN_LSRsaturated=zeros(nLevels,nFrequencies);
CNLSRsaturated=zeros(nLevels,nFrequencies);
CNHSRsaturated=zeros(nLevels,nFrequencies);
ICHSRsaturated=zeros(nLevels,nFrequencies);
ICLSRsaturated=zeros(nLevels,nFrequencies);


levelNo=0; PSTHplotCount=0;
for leveldB=levels
    fprintf('%4.0f\t', leveldB)
    levelNo=levelNo+1;
    amp=28e-6*10^(leveldB/20);
    
    freqNo=0;
    for frequency=frequencies

        paramFunctionName=['method=MAPparams' experiment.name ...
            '(' num2str(targetFrequency) ');' ];
        eval(paramFunctionName);  % read parameters afresh each pass

        dt=method.dt;
        time=dt:dt:toneDuration;        
        rampTime=dt:dt:rampDuration;
        ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(time)-length(rampTime))];
        ramp=ramp.*fliplr(ramp);
        
        silence=zeros(1,round(silenceDuration/dt));
        
        toneStartptr=length(silence)+1;
        toneMidptr=toneStartptr+round(toneDuration/(2*dt)) -1;
        toneEndptr=toneStartptr+round(toneDuration/dt) -1;
        
        % create signal (leveldB/ frequency)
        freqNo=freqNo+1;
        inputSignal=amp*sin(2*pi*frequency'*time);
        inputSignal= ramp.*inputSignal;
        inputSignal=[silence inputSignal silence];
        
        if doReceptiveFields  % receptive field
            method.plotGraphs=	0; % plot only PSTHs
        else
            method.plotGraphs=	plotGraphsForIHC; % show progress
        end
        
        targetChannelNo=1;
        
        % force parameters
         % the number of AN fibers at each BF
        AN_IHCsynapseParams.numFibers=	nRepeats;
        AN_IHCsynapseParams. mode= 'spikes';
        AN_IHCsynapseParams.plotSynapseContents=0;
        AN_IHCsynapseParams.PSTHbinWidth=.001;
        
        method.DRNLSave=1;
        method.IHC_cilia_RPSave=1;
        method.PSTHbinWidth=1e-3; % useful 1-ms default for all PSTHs
        method.AN_IHCsynapseSave=1;
        method.MacGregorMultiSave=1;
        method.MacGregorSave=1;
        method.dt=dt;
              
        moduleSequence=[1:8];       

        global dtSpikes ARAttenuation TMoutput OMEoutput ...
    DRNLoutput IHC_cilia_output IHCrestingCiliaCond IHCrestingV...
    IHCoutput ANprobRateOutput ANoutput savePavailable ANtauCas  ...
    CNoutput  ICoutput ICmembraneOutput ICfiberTypeRates MOCattenuation

AN_spikesOrProbability='spikes';
AN_spikesOrProbability='probability';
MAPparamsName='Normal';

MAP1_14(inputSignal, 1/dt, targetFrequency, ...
    MAPparamsName, AN_spikesOrProbability);
        
        % RP
        IHC_RPData=IHC_cilia_output;
        IHC_RPData=IHCoutput(targetChannelNo,:);
        IHC_RP_peak(levelNo,freqNo)=max(IHC_RPData(toneStartptr:toneEndptr));
        IHC_RP_min(levelNo,freqNo)=min(IHC_RPData(toneStartptr:toneEndptr));
        IHC_RP_dc(levelNo,freqNo)=mean(IHC_RPData(toneStartptr:toneEndptr));
        
        % AN next
        AN_IHCsynapseAllData=ANoutput;
        method.PSTHbinWidth=0.001;
        
        nTaus=length(ANtauCas);
        numANfibers=size(ANoutput,1);
        numLSRfibers=numANfibers/nTaus;
        
        %LSR (same as HSR if no LSR fibers present)
        channelPtr1=(targetChannelNo-1)*numANfibers+1;
        channelPtr2=channelPtr1+numANfibers-1;
        LSRspikes=AN_IHCsynapseAllData(channelPtr1:channelPtr2,:);
        method.dt=method.AN_IHCsynapsedt;
        PSTH=UTIL_PSTHmaker(LSRspikes, method);
        PSTH=sum(PSTH,1); % sum across fibers (HSR only)
        PSTHStartptr=round(silenceDuration/method.PSTHbinWidth)+1;
        PSTHMidptr=PSTHStartptr+round(toneDuration/(2*method.PSTHbinWidth)) -1;
        PSTHEndptr=PSTHStartptr+round(toneDuration/method.PSTHbinWidth) -1;
        AN_LSRonset(levelNo,freqNo)=max(max(PSTH))/(method.PSTHbinWidth*method.numANfibers);
        AN_LSRsaturated(levelNo,freqNo)=sum(PSTH(PSTHMidptr:PSTHEndptr))/(method.numANfibers*toneDuration/2);
        
        % HSR
        channelPtr1=numLSRfibers+(targetChannelNo-1)*method.numANfibers+1;
        channelPtr2=channelPtr1+method.numANfibers-1;
        HSRspikes=AN_IHCsynapseAllData(channelPtr1:channelPtr2,:);
        method.dt=method.AN_IHCsynapsedt;
        PSTH=UTIL_PSTHmaker(HSRspikes, method);
        PSTH=sum(PSTH,1); % sum across fibers (HSR only)
        PSTHStartptr=round(silenceDuration/method.PSTHbinWidth)+1;
        PSTHMidptr=PSTHStartptr+round(toneDuration/(2*method.PSTHbinWidth)) -1;
        PSTHEndptr=PSTHStartptr+round(toneDuration/method.PSTHbinWidth) -1;
        AN_HSRonset(levelNo,freqNo)=max(max(PSTH))/(method.PSTHbinWidth*method.numANfibers);
        AN_HSRsaturated(levelNo,freqNo)=sum(PSTH(PSTHMidptr:PSTHEndptr))/(method.numANfibers*toneDuration/2);
        [cvANHSR, cvTimes, allTimeStamps, allISIs]=  UTIL_CV(HSRspikes, method.AN_IHCsynapsedt);
        
        PSTHplotCount=PSTHplotCount+1;
        if doReceptiveFields  % receptive field for HSR only
            figure(14), set(gcf,'name','AN')
            plotReceptiveFields(method, PSTH, PSTHplotCount, levels, frequencies)
            ylim([0 method.numANfibers])
            xlabel(['CV= ' num2str(max(cvANHSR),'%4.2f')],'fontsize',8)
        end % doReceptiveFields
        
        % CN
        MacGregorMultiAllData=method.MacGregorMultiData;
        numLSRfibers=method.McGMultinNeuronsPerBF*length(method.nonlinCF)* (nTaus-1);
        
        %LSR (same as HSR if no LSR fibers present)
        channelPtr1=(targetChannelNo-1)*method.McGMultinNeuronsPerBF+1;
        channelPtr2=channelPtr1+method.McGMultinNeuronsPerBF-1;
        MacGregorMultiLSRspikes=MacGregorMultiAllData(channelPtr1:channelPtr2,:);
        method.dt=method.MacGregorMultidt;
        PSTH=UTIL_PSTHmaker(MacGregorMultiLSRspikes, method);
        PSTH=sum(PSTH,1); % sum across fibers (HSR only)
        PSTHStartptr=round(silenceDuration/method.PSTHbinWidth)+1;
        PSTHMidptr=PSTHStartptr+round(toneDuration/(2*method.PSTHbinWidth)) -1;
        PSTHEndptr=PSTHStartptr+round(toneDuration/method.PSTHbinWidth) -1;
        CNLSRsaturated(levelNo,freqNo)=sum(PSTH(PSTHMidptr:PSTHEndptr));
        CNLSRsaturated(levelNo,freqNo)=CNLSRsaturated(levelNo,freqNo)...
            /((toneDuration/2)*method.McGMultinNeuronsPerBF);
        
        %HSR
        channelPtr1=numLSRfibers+(targetChannelNo-1)*method.McGMultinNeuronsPerBF+1;
        channelPtr2=channelPtr1+method.McGMultinNeuronsPerBF-1;
        MacGregorMultiHSRspikes=MacGregorMultiAllData(channelPtr1:channelPtr2,:);
        method.dt=method.MacGregorMultidt;
        PSTH=UTIL_PSTHmaker(MacGregorMultiHSRspikes, method);
        PSTH=sum(PSTH,1); % sum across fibers (HSR only)
        PSTHStartptr=round(silenceDuration/method.PSTHbinWidth)+1;
        PSTHMidptr=PSTHStartptr+round(toneDuration/(2*method.PSTHbinWidth)) -1;
        PSTHEndptr=PSTHStartptr+round(toneDuration/method.PSTHbinWidth) -1;
        CNHSRsaturated(levelNo,freqNo)=sum(PSTH(PSTHMidptr:PSTHEndptr));
        CNHSRsaturated(levelNo,freqNo)=CNHSRsaturated(levelNo,freqNo)...
            /((toneDuration/2)*method.McGMultinNeuronsPerBF);
        [cvMMHSR, cvTimes, allTimeStamps, allISIs]=  UTIL_CV(MacGregorMultiHSRspikes, method.MacGregorMultidt);
        
        if doReceptiveFields  % receptive field
            figure(15), set(gcf,'name','CN HSR input')
            plotReceptiveFields(method, PSTH, PSTHplotCount, levels, frequencies)
            ylim([0 method.McGMultinNeuronsPerBF])
            xlabel(['CV= ' num2str(max(cvMMHSR),'%4.2f')],'fontsize',8)
        end
        
        MacGregorAllData=method.MacGregorData;
        numLSRfibers=length(method.nonlinCF)* (nTaus-1);
        
        %LSR (same as HSR if no LSR fibers present)
        channelPtr1=targetChannelNo;
        MacGregorLSR=MacGregorAllData(channelPtr1,:);
        method.dt=method.MacGregordt;
        PSTH=UTIL_PSTHmaker(MacGregorLSR, method);
        PSTHStartptr=round(silenceDuration/method.PSTHbinWidth)+1;
        PSTHMidptr=PSTHStartptr+round(toneDuration/(2*method.PSTHbinWidth)) -1;
        PSTHEndptr=PSTHStartptr+round(toneDuration/method.PSTHbinWidth) -1;
        ICLSRsaturated(levelNo,freqNo)=sum(PSTH(PSTHMidptr:PSTHEndptr));
        ICLSRsaturated(levelNo,freqNo)=ICLSRsaturated(levelNo,freqNo)/(toneDuration/2);
        
        %LSR (same as HSR if no LSR fibers present)
        channelPtr1=numLSRfibers+targetChannelNo;
        MacGregorHSR=MacGregorAllData(channelPtr1,:);
        method.dt=method.MacGregordt;
        PSTH=UTIL_PSTHmaker(MacGregorHSR, method);
        PSTHStartptr=round(silenceDuration/method.PSTHbinWidth)+1;
        PSTHMidptr=PSTHStartptr+round(toneDuration/(2*method.PSTHbinWidth)) -1;
        PSTHEndptr=PSTHStartptr+round(toneDuration/method.PSTHbinWidth) -1;
        ICHSRsaturated(levelNo,freqNo)=sum(PSTH(PSTHMidptr:PSTHEndptr));
        ICHSRsaturated(levelNo,freqNo)=ICHSRsaturated(levelNo,freqNo)/(toneDuration/2);
        [cvICHSR, cvTimes, allTimeStamps, allISIs]=  UTIL_CV(MacGregorHSR, method.MacGregordt);
        
%         if doReceptiveFields  % receptive field
%             figure(16), set(gcf,'name','IC HSR input')
%             plotReceptiveFields(method, PSTH, PSTHplotCount, levels, frequencies)
%             ylim([0 method.McGMultinNeuronsPerBF])
%             xlabel(['CV= ' num2str(max(cvICHSR),'%4.2f')],'fontsize',8)
%         end
    end % frequency
end % level
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


% ---------------------------------------------------------- display parameters
disp(['parameter file was: ' experiment.name])
fprintf('\n')
UTIL_showStruct(IHC_VResp_VivoParams, 'IHC_cilia_RPParams')
UTIL_showStruct(IHCpreSynapseParams, 'IHCpreSynapseParams')
UTIL_showStruct(AN_IHCsynapseParams, 'AN_IHCsynapseParams')



function plotReceptiveFields(method, PSTH, PSTHplotCount,  levels, frequencies)

% show PSTH for each level/frequency combination
nLevels=length(levels);
nFrequencies=length(frequencies);

PSTHtime=method.PSTHbinWidth:method.PSTHbinWidth:method.PSTHbinWidth*length(PSTH);
subplot(nLevels,nFrequencies,PSTHplotCount)
bar(PSTHtime, PSTH)
xlim([0 max(PSTHtime)])
% write axis labels only at left and bottom
if PSTHplotCount< (nLevels-1) * nFrequencies+1
    set(gca,'xticklabel',[])
end
if ~isequal(mod(PSTHplotCount,nFrequencies),1)
    set(gca,'yticklabel',[])
else
    ylabel([num2str(levels(round(PSTHplotCount/nFrequencies) +1)) ' dB'])
end
% add titles only on top row
if PSTHplotCount<=nFrequencies
    title([num2str(frequencies(PSTHplotCount)) ' Hz'])
end
