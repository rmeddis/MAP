function testRP2(MAPparamsName,paramChanges)
% testIHC evaluates IHC I/O function
% multiple BF can be used but only one is easier to interpret.
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
    BF=800;
end

levels=-20:10:100;
nLevels=length(levels);
toneDuration=.05;
silenceDuration=.01;
sampleRate=50000;
dt=1/sampleRate;

allIHC_RP_peak=[];
allIHC_RP_dc=[];

%% Ruggero
%%Ruggero data
RuggeroData=[
    0	2.00E-10;
    10	5.00E-10;
    20	1.50E-09;
    30	2.50E-09;
    40	5.30E-09;
    50	1.00E-08;
    60	1.70E-08;
    70	2.50E-08;
    80	4.00E-08;
    90	6.00E-08;
    100	1.50E-07;
    110	3.00E-07;
    ];


BF=10000;
targetFrequency=BF;

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
    
end
%% plot DRNL
subplot(2,2,1)
referenceDisp=10e-9;
semilogy(levels,DRNL_peak, 'linewidth',2), hold on
semilogy(RuggeroData(:,1),RuggeroData(:,2),'o')
title(['BM: Ruggero ' num2str(BF) ' Hz'])
ylabel ('displacement(m)'), xlabel('dB SPL')
xlim([min(levels) max(levels)]), ylim([1e-10 1e-7])
grid on


%% Dallos
BF=800;
targetFrequency=BF;

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
    
end

%% Dallos and Harris data
%% plot receptor potentials
figure(4)
subplot(2,2,3)
% RP I/O function min and max
restingRP=IHC_RP_peak(1);
toPlot= [fliplr(IHC_RP_min(:,1)') IHC_RP_peak(:,1)'];
microPa=   28e-6*10.^(levels/20);
microPa=[-fliplr(microPa) microPa];
plot(microPa,toPlot, 'linewidth',2)

dallosx=[-0.9	-0.1	-0.001	0.001	0.01	0.9];
dallosy=[-8	-7.8	-6.5	11	16.5	22]/1000 + restingRP;
subplot(2,2,3)
hold on, plot(dallosx,dallosy, 'o')
plot([-1 1], [restingRP restingRP], 'r')
title(' Dallos(86) data at 800 Hz')
ylabel ('receptor potential(V)'), xlabel('Pa')
ylim([-0.08 -0.02]), xlim([-1 1])
grid on


%% Patuzzi
BF=7000;
targetFrequency=BF;

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
    
    %% run the model
    global  DRNLoutput IHC_cilia_output IHCrestingCiliaCond IHCrestingV...
        IHCoutput
    AN_spikesOrProbability='probability';
    
    MAP1_14(inputSignal, sampleRate, BF, ...
        MAPparamsName, AN_spikesOrProbability, paramChanges);
    
    % RP
    IHC_RPData=IHCoutput;
    IHC_RP_peak(levelNo,1)=...
        max(IHC_RPData(toneMidptr:toneEndptr));
    IHC_RP_min(levelNo,1)=...
        min(IHC_RPData(toneMidptr:toneEndptr));
    IHC_RP_dc(levelNo,1)=...
        mean(IHC_RPData(toneMidptr:toneEndptr));
end


%% RP I/O function min and max
figure(4)
subplot(2,2,4)
restingRP=IHC_RP_peak(1);
peakRP=max(IHC_RP_peak);
plot(levels, IHC_RP_peak, 'linewidth',2)
hold on
plot(levels, IHC_RP_dc, ':', 'linewidth',2)
plot([min(levels) max(levels)], [restingRP restingRP], 'r')
xlim([min(levels) max(levels)])
% animal data
sndLevel=[5	15	25	35	45	55	65	75];
RPanimal=restingRP+[0.5	2	4.6	5.8	6.4	7.2	8	10.2]/1000;
% could be misleading when restingRP changes
RPanimal=-0.060+[0.5	2	4.6	5.8	6.4	7.2	8	10.2]/1000;
hold on, plot(sndLevel,RPanimal,'o')

grid on
title(' 7 kHz Patuzzi')
ylabel ('RP(V) peak and DC'), xlabel('dB SPL')
ylim([-0.07 -0.04])
path(savePath);
disp(paramChanges)
