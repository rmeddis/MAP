function plotProfile(longTone,shortTone,gaps,BFs,TMC,offBFs,IFMCs)

%% plot profile
if nargin<1
    load profile
end

% comparison data (e.g. participants)
% rows are BFs

% % -------------------------------------------JSan
% compareBFs=[250	500	1000	2000	3000
% ];
% compareLongTone=  [67.46485	56.95655	65.01985	61.46655	73.33265
% ];
% compareShortTone=[	72.3185	63.2818	69.0373	65.2853	76
% ];
% 
% compareGaps=[0.02 0.05 0.08];
% compareTMC=	[
%     84	69	77	75	93
% 88	73	81	79	97
% 95	79	85	83	98
% ];
% 
% compareMaskerFreqs=[0.7  0.9 1 1.1  1.3 ];
% compareIFMCs=[83.1698	77.3165	79.8474	82.9074	82.3294
% 80.9667	73.6653	80.9446	80.7005	79.0022
% 82.0135	71.2284	78.7345	74.3342	84.126
% 79.3348	70.3347	78.5769	79.7274	91.9849
% 79.2308	76.3164	83.6881	86.2105	NaN
%   ];
% 
% -------------------------------------------JE
compareBFs=[250 500 1000 2000 4000 8000];
compareLongTone=  [32 30	31	40	54 NaN];
compareShortTone=[	49	50	47	56	63	NaN];

compareGaps=0.01:0.01:0.09;
compareTMC=	[
    69	83	82	NaN	NaN	NaN	NaN	NaN	NaN
    61	68	79	88	93	NaN	NaN	NaN	NaN
    63	69	79	84	92	NaN	NaN	NaN	NaN
    67	71	75	80	82	84	88	93	NaN
    82	82	86	86	NaN	83	88	90	75
    NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
    ];
compareTMC=compareTMC';

compareMaskerFreqs=[0.5	0.7	0.9	1	1.1	1.3	1.6];
compareIFMCs=[
    64	60	58	58	57	60	64
    65	63	58	55	54	59	69
    68	64	60	59	62	73	79
    76	75	71	67	68	71	77
    79	71	68	69	73	75	77
    76	73	75	75	76	80	NaN
    ];
compareIFMCs=compareIFMCs';

% % -------------------------------------------CMR
% % CMR 
% compareBFs=[250 500 1000 2000 4000 8000];
% compareLongTone=[	11.4	1.55	-13.5	-6.35	-6.4	7.45];
% compareShortTone=[	23.85	18.9	9.85	10.6	9.55	21.9];
% 
% compareGaps=0.01:0.01:0.09;
% compareTMC=	[
%     28.5	35.0	49.3	70.1	80.5	85.5    NaN      NaN    NaN;
%     31.1	44.3	48.4	59.5	56.4	76.7	70.2	82.4	76.3;
%     33.4	38.4	48.8	55.8	64.5	78.7	84.2	88.3	90.3;
%     25.4	37.0	49.2	49.7	58.2	69.6	87.7	95.8	93.0;
%     18.2	23.5	27.4	41.5	64.3	82.1	86.7	91.2	NaN;
%     32.5	35.8	43.5	52.1	69.1	78.6	86.6	86.0	NaN;
%     ];
% compareTMC=compareTMC';
% 
% compareMaskerFreqs=[0.5	0.7	0.9	1	1.1	1.3	1.6];
% compareIFMCs=[
%     50	42	34	35	34	33	37;
%     58	51	38	33	28	41	49;
%     57	41	27	20	28	37	66;
%     61	49	27	20	34	68	79;
%     67	45	27	22	46	74	87;
%     62	62	43	22	47	56	83;
%     ];
% compareIFMCs=compareIFMCs';

% absolute thresholds
figure(90), clf
subplot(2,1,2)
semilogx(BFs,longTone,'ko-','lineWidth',2); hold on
semilogx(BFs,shortTone,'bo-','lineWidth',2); hold on
semilogx(compareBFs,compareLongTone,'ko:'); hold on
semilogx(compareBFs,compareShortTone,'bo:'); hold on
ylim([0 100])

% TMC
for BFno=1:length(BFs)
    subplot(2,6,BFno)
    plot(gaps,TMC(:,BFno)-longTone(BFno),'r','lineWidth',3), hold on
    plot(gaps,TMC(:,BFno),'b','lineWidth',3), hold on
    ylim([-10 110])
    xlim([0.01 0.1])
    grid on
    if BFno==1
        ylabel('masker dB SL')
        xlabel('gap')
%         text(0.02,80,' TMC','backgroundColor','w')
    end
    title([num2str(BFs(BFno)) ' Hz'])
    set(gca,'XTick',[ 0.1],'xTickLabel', { '0.1'})
end

% IFMCs
for BFno=1:length(BFs)
    freq=offBFs'*BFs(BFno);
    subplot(2,1,2)
    semilogx(freq,IFMCs(:,BFno),'r','lineWidth',3), hold on
    ylim([0 100])
    xlim([100 12000])
    grid on
end
xlabel('frequency (Hz)')
ylabel('masker dB / probe dB')
set(gca,'XTick',BFs)

for BFno=1:length(compareBFs)
    subplot(2,6,BFno)
    plot(compareGaps,compareTMC(:,BFno)-longTone(BFno),'k:')
    plot(compareGaps,compareTMC(:,BFno),'k:')
    ylim([-10 110])
    xlim([0.01 0.1])
    grid on
    if BFno==1
        ylabel('masker dB SL')
        xlabel('gap')
%         text(0.02,80,' TMC','backgroundColor','w')
    end
    title([num2str(BFs(BFno)) ' Hz'])
    set(gca,'XTick',[ 0.1],'xTickLabel', { '0.1'})
end

% IFMCs
for BFno=1:length(compareBFs)
    compareFreq=compareMaskerFreqs'*BFs(BFno);
    subplot(2,1,2)
    semilogx(compareFreq,compareIFMCs(:,BFno),'k:')
    ylim([0 100])
    xlim([100 12000])
    grid on
end
mydate=datestr(now); idx=findstr(':',mydate); mydate(idx)='_';

fileName= ['savedData/' mydate ];

save (fileName)
set(gcf,'name', mydate)
disp(fileName)
