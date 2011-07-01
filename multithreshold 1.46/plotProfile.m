function plotProfile(longTone,shortTone,gaps,BFs,TMC,offBFs,IFMCs)

%% plot profile
if nargin<1
    load profile
end

normLongTone=[	11.4	1.55	-13.5	-6.35	-6.4	7.45];
normShortTone=[	23.85	18.9	9.85	10.6	9.55	21.9];

normGaps=0.01:0.01:0.09;
normTMC=	[
28.5	35.0	49.3	70.1	80.5	85.5    NaN      NaN    NaN;			
31.1	44.3	48.4	59.5	56.4	76.7	70.2	82.4	76.3;
33.4	38.4	48.8	55.8	64.5	78.7	84.2	88.3	90.3;
25.4	37.0	49.2	49.7	58.2	69.6	87.7	95.8	93.0;
18.2	23.5	27.4	41.5	64.3	82.1	86.7	91.2	NaN;
32.5	35.8	43.5	52.1	69.1	78.6	86.6	86.0	NaN;
    ];
normTMC=normTMC';

normIFMC=[
50	42	34	35	34	33	37;
58	51	38	33	28	41	49;
57	41	27	20	28	37	66;
61	49	27	20	34	68	79;
67	45	27	22	46	74	87;
62	62	43	22	47	56	83;
];
normIFMC=normIFMC';

% absolute thresholds
figure(90), clf
subplot(2,1,2)
semilogx(BFs,longTone,'ko-','lineWidth',2); hold on
semilogx(BFs,shortTone,'bo-','lineWidth',2); hold on
semilogx(BFs,normLongTone,'ko:'); hold on
semilogx(BFs,normShortTone,'bo:'); hold on
ylim([0 100])

% TMC
for BFno=1:length(BFs)
    subplot(2,6,BFno)
    plot(gaps,TMC(:,BFno)-longTone(BFno),'r','lineWidth',3), hold on
    plot(normGaps,normTMC(:,BFno)-longTone(BFno),'k:')
    ylim([-10 90])
    xlim([0.01 0.1])
    grid on    
    if BFno==1
        ylabel('masker dB SL')
        xlabel('gap')
        text(0.02,80,' TMC','backgroundColor','w')
    end
        title([num2str(BFs(BFno)) ' Hz'])
        set(gca,'XTick',[ 0.1],'xTickLabel', { '0.1'})
end

% IFMCs
for BFno=1:length(BFs)
    BF=BFs(BFno);
    freq=offBFs'*BFs(BFno);
    subplot(2,1,2)
    semilogx(freq,IFMCs(:,BFno),'r','lineWidth',3), hold on
    semilogx(freq,normIFMC(:,BFno),'k:')
    ylim([0 100])
    xlim([100 12000])
    grid on
end
xlabel('frequency (Hz)')
ylabel('masker dB / probe dB')
set(gca,'XTick',BFs)
