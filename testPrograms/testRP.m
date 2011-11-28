function testRP(BFs,MAPparamsName,paramChanges)
% testIHC evaluates IHC I/O function
% multiple BFs can be used but only one is easier to interpret.
% e.g. testRP(1000,'Normal',{});

global experiment method inputStimulusParams
global stimulusParameters IHC_VResp_VivoParams IHC_cilia_RPParams
savePath=path;
addpath (['..' filesep 'utilities'],['..' filesep 'MAP'])
dbstop if error

figure(4), clf,
set (gcf, 'name', ['IHC'])
set(gcf,'position',[613   354   360   322])
drawColors='rgbkmcy';
drawnow

if nargin<3
    paramChanges=[]; 
end
if nargin<2
    MAPparamsName='Normal';
end
if nargin<3
    BFs=800; 
end

levels=-20:10:100;
nLevels=length(levels);
toneDuration=.05;
silenceDuration=.01;
sampleRate=50000;
dt=1/sampleRate;

allIHC_RP_peak=[];
allIHC_RP_dc=[];

for BFno=1:length(BFs)
    BF=BFs(BFno);
    targetFrequency=BF;
    % OR
    %Patuzzi and Sellick test (see ELP & AEM, 2006)
    % targetFrequency=100;

    IHC_RP_peak=zeros(nLevels,1);
    IHC_RP_min=zeros(nLevels,1);
    IHC_RP_dc=zeros(nLevels,1);

    time=dt:dt:toneDuration;

    rampDuration=0.004;
    rampTime=dt:dt:rampDuration;
    ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
        ones(1,length(time)-length(rampTime))];
    ramp=ramp.*fliplr(ramp);

    silence=zeros(1,round(silenceDuration/dt));

    toneStartptr=length(silence)+1;
    toneMidptr=toneStartptr+round(toneDuration/(2*dt)) -1;
    toneEndptr=toneStartptr+round(toneDuration/dt) -1;

    levelNo=0;
    for leveldB=levels
        levelNo=levelNo+1;
        % replicate at all levels
        amp=28e-6*10^(leveldB/20);

        %% create signal (leveldB/ frequency)
        inputSignal=amp*sin(2*pi*targetFrequency'*time);
        inputSignal= ramp.*inputSignal;
        inputSignal=[silence inputSignal silence];
        inputStimulusParams.sampleRate=1/dt;
%         global IHC_ciliaParams

        %% disable efferent for fast processing
        method.DRNLSave=1;
        method.IHC_cilia_RPSave=1;
        method.IHCpreSynapseSave=1;
        method.IHC_cilia_RPSave=1;
        method.segmentDuration=-1;
        moduleSequence=1:4;

        %% run the model
        global  DRNLoutput IHC_cilia_output IHCrestingCiliaCond IHCrestingV...
            IHCoutput
        AN_spikesOrProbability='probability';

        MAP1_14(inputSignal, sampleRate, BF, ...
            MAPparamsName, AN_spikesOrProbability, paramChanges);

        % DRNL
        DRNLoutput=DRNLoutput;
        DRNL_peak(levelNo,1)=max(DRNLoutput(toneMidptr:toneEndptr));
        DRNL_min(levelNo,1)=min(DRNLoutput(toneMidptr:toneEndptr));
        DRNL_dc(levelNo,1)=mean(DRNLoutput(toneMidptr:toneEndptr));

        % cilia
        IHC_ciliaData=IHC_cilia_output;
        IHC_ciliaData=IHC_ciliaData;
        IHC_cilia_peak(levelNo,1)=...
            max(IHC_ciliaData(toneMidptr:toneEndptr));
        IHC_cilia_min(levelNo,1)=...
            min(IHC_ciliaData(toneMidptr:toneEndptr));
        IHC_cilia_dc(levelNo,1)=...
            mean(IHC_ciliaData(toneMidptr:toneEndptr));

        % RP
        IHC_RPData=IHCoutput;
        IHC_RPData=IHC_RPData;
        IHC_RP_peak(levelNo,1)=...
            max(IHC_RPData(toneMidptr:toneEndptr));
        IHC_RP_min(levelNo,1)=...
            min(IHC_RPData(toneMidptr:toneEndptr));
        IHC_RP_dc(levelNo,1)=...
            mean(IHC_RPData(toneMidptr:toneEndptr));

    end % level


    disp(['parameter file was: ' MAPparamsName])
    fprintf('\n')

    %% plot DRNL
    subplot(2,2,1)
%     referenceDisp= 9e-9*1000/targetFrequency;
%     plot(levels,20*log10(DRNL_peak/referenceDisp), drawColors(BFno), ...
%         'linewidth',2), hold on
    referenceDisp=10e-9;
    plot(levels,20*log10(DRNL_peak/referenceDisp), drawColors(BFno), ...
        'linewidth',2), hold on
    title([' DRNL peak:  ' num2str(BFs) ' Hz'])
    ylabel ('log10DRNL(m)'), xlabel('dB SPL')
    xlim([min(levels) max(levels)]), ylim([-10 50])
    grid on

    %% plot cilia displacement
    figure(4)
    subplot(2,2,2)
    restingIHC_cilia=IHCrestingCiliaCond;
    plot(levels, IHC_cilia_peak,'k', 'linewidth',2), hold on
%     plot(levels, IHC_cilia_min,'r', 'linewidth',2)
%     hold on,
%     plot([min(levels) max(levels)], ...
%         [restingIHC_cilia restingIHC_cilia], 'g')
    title(' IHC apical cond.')
    ylabel ('IHCcilia(conductance)'), xlabel('dB SPL')
    xlim([min(levels) max(levels)])
    grid on

    %% plot receptor potentials
    figure(4)
    subplot(2,2,3)
    % RP I/O function min and max
    restingRP=IHC_RP_peak(1);
    toPlot= [fliplr(IHC_RP_min(:,1)') IHC_RP_peak(:,1)'];
    microPa=   28e-6*10.^(levels/20);
    microPa=[-fliplr(microPa) microPa];
    plot(microPa,toPlot, drawColors(BFno), 'linewidth',2)
    % ylim([0 300])

    %% Dallos and Harris data
    dallosx=[-0.9	-0.1	-0.001	0.001	0.01	0.9];
    dallosy=[-8	-7.8	-6.5	11	16.5	22]/1000 + restingRP;
    hold on, plot(dallosx,dallosy, 'o')
    plot([-1 1], [restingRP restingRP], 'r')
    title(' Dallos(86) data at 800 Hz')
    ylabel ('receptor potential(V)'), xlabel('Pa')
    ylim([-0.08 -0.02]), xlim([-1 1])
    grid on

    %% RP I/O function min and max
    figure(4)
    subplot(2,2,4)
    restingRP=IHC_RP_peak(1);
    peakRP=max(IHC_RP_peak);
    plot(levels, IHC_RP_peak,drawColors(BFno), 'linewidth',2)
    hold on
    plot(levels, IHC_RP_dc, [drawColors(BFno) ':'], 'linewidth',2)
    hold on,
    plot([min(levels) max(levels)], [restingRP restingRP], 'r')
    xlim([min(levels) max(levels)])
    % animal data
    sndLevel=[5	15	25	35	45	55	65	75];
RPanimal=restingRP+[0.5	2	4.6	5.8	6.4	7.2	8	10.2]/1000;
% could be misleading when restingRP changes
RPanimal=-0.060+[0.5	2	4.6	5.8	6.4	7.2	8	10.2]/1000;
hold on, plot(sndLevel,RPanimal,'o')

    grid on
    title(['Et= ' num2str(IHC_cilia_RPParams.Et) ':  RP data 7 kHz Patuzzi'])
    ylabel ('RP(V) peak and DC'), xlabel('dB SPL')
    ylim([-0.08 -0.04])
    allIHC_RP_peak=[allIHC_RP_peak IHC_RP_peak];
    allIHC_RP_dc=[allIHC_RP_dc IHC_RP_dc];

    fprintf('level\t peak\t DC\n')
    UTIL_printTabTable([levels' IHC_RP_peak IHC_RP_dc])
    % disp(['restingIHC_cilia= ' num2str(restingIHC_cilia)])
    fprintf('peakRP= \t%6.3f', peakRP)
    fprintf('\nrestingRP= \t%6.3f', restingRP)
    fprintf('\ndifference= \t%6.3f\n', (peakRP-restingRP))
    drawnow
end
% UTIL_showStruct(IHC_VResp_VivoParams, 'IHC_VResp_VivoParams')
UTIL_showStruct(IHC_cilia_RPParams, 'IHC_cilia_RPParams')
    fprintf('level\t peak\n')
    UTIL_printTabTable([levels' allIHC_RP_peak])
    fprintf('level\t DC\n')
    UTIL_printTabTable([levels' allIHC_RP_dc])

path(savePath);
disp(paramChanges)
