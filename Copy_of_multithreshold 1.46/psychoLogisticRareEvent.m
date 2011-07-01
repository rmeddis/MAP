function psycho 

% DATA from AbsThresholdCM, LS

levelSequence=[16.0836	6.0836	11.0836	12.725	13.4995	11.6559	12.7168	13.9507	12.0632	11.3845	11.7754	12.8639	14.1061	14.5545	13.5589	13.4736	14.1554	13.5123	13.1537	11.5334	12.6921	11.4461	10.7112	11.9171	10.8412	12.0536	13.06	14.2345	13.8535	14.0575	14.0145	15.6385	13.7959	14.0397	12.2497	11.4383	12.7685	12.7317	11.6741	10.3955	9.0334	10.4787	11.1883	11.4687	12.0583	13.8982	13.0644	14.4722	14.0304	14.7569	13.8571	12.1521	14.0865	13.2258	11.6988	12.2186	12.2561	12.2711	13.0345	14.2948	13.7838	12.9155	12.6927	12.5757	13.4818	12.8397	14.1803	12.7524	12.3169	13.0973	14.9722	14.4422	15.8944	17.2728	16.8621	15.3074	13.7319	11.7716	12.1233	12.1576	13.0023	13.2442	12.1454	13.8298	13.3997	11.7799	12.3834	12.8387	14.3718	12.7063	14.4002	14.3674	14.7882	15.0271	13.5737	15.2095	14.1789	13.4982	14.3926	13.6277	15.0551	13.493	15.3368	16.9437	16.1645	15.9275	15.7577	14.2939	14.1069	13.6184	12.728	14.4283	15.6494	14.2216	15.8188	14.3906	12.7635	13.2973	11.9934	11.1053	12.2359	11.439	12.4735	12.5466	12.9186	13.6937	13.792	13.085	12.5016	13.8401	11.9688	13.0676	14.3708	13.4891	14.4394	15.8265	17.1478	15.2337	13.6512	12.7727	14.4818	13.7045	14.5207	13.8211	13.2874	14.9221	14.3822	14.3562	13.3132	13.5015	13.1372	14.4495	13.957	15.2124	16.6607	15.114	13.8978	12.6665	12.3469	13.3699	14.1991	15.2259	13.2326	11.9325	13.1954	12.2766	13.7996	15.5274	14.7623	13.21	12.5874	10.6431	12.2121	10.4043	11.9979	10.2153	10.7945	11.5215	12.7081	14.2288	14.4475	13.9777	13.8564	14.1923	15.393	14.5229	15.3377	13.9385	13.6328	14.2861	15.0021	15.5488	13.8067	13.0269	13.3644	13.6195	14.5834	15.6554	13.8914	13.3657	14.0566	12.9857	12.9343	14.0372	14.0156	15.6395	15.4225	14.6514	14.6498	13.6421	13.4263	13.227	14.0428	14.681	15.1661	13.5285	14.9074	13.5498	14.4265	14.0798	14.6391	12.8226	12.225	13.7883	12.7012	12.7932	13.5318	14.9207	14.064	12.9551	14.2066	14.3463	13.4565	13.9379	15.7189	14.3516	12.6425	12.8485	12.2879	13.1142	11.8093	12.0037	10.2154	11.8395	12.4643	13.0501	13.5452	13.2338	12.6246	13.292	11.5661	11.6015	13.2896	15.2506	15.5295	14.523	13.6199	13.083	12.1714	12.5977	14.2862	15.7924	14.0221	14.4042	16.0364	16.6967	14.7112	14.7594	13.9874	13.3428	13.5231	14.441	13.9983	13.6414	11.862	13.2632	13.1065	12.8701	12.2578	13.7517	14.3268	13.3849	13.9609	13.3239	14.6694	14.8172	13.654	13.8144	13.5015	14.34	14.0139	15.1029	13.3895	13.3701	13.4554	11.9646	12.2519	12.4145	13.1343	14.9372	13.2483
];


responseSequence=[1	0	0	0	1	0	0	1	1	0	0	0	0	1	1	0	1	1	1	0	1	1	0	1	0	0	0	1	0	1	0	1	0	1	1	0	1	1	1	1	0	0	0	0	0	1	0	1	0	1	1	0	1	1	0	0	0	0	0	1	1	1	1	0	1	0	1	1	0	0	1	0	0	1	1	1	1	0	0	0	0	1	0	1	1	0	0	0	1	0	1	0	0	1	0	1	1	0	1	0	1	0	0	1	1	1	1	1	1	1	0	0	1	0	1	1	0	1	1	0	1	0	0	0	0	0	1	1	0	1	0	0	1	0	0	0	1	1	1	0	1	0	1	1	0	1	1	1	0	1	0	1	0	0	1	1	1	1	0	0	0	1	1	0	1	0	0	1	1	1	1	0	1	0	1	0	0	0	0	0	1	1	0	0	1	0	1	1	0	0	0	1	1	0	0	0	0	1	1	0	1	1	0	1	0	1	1	1	1	1	1	0	0	0	1	0	1	0	1	0	1	1	0	1	0	0	0	1	1	0	0	1	0	0	1	1	0	1	0	1	0	1	0	0	0	0	1	1	0	1	0	0	0	0	1	1	1	1	0	0	0	1	0	0	0	1	0	1	1	0	0	1	1	1	0	1	1	1	0	0	1	0	1	0	0	1	0	1	0	1	0	1	1	0	1	0	0	0	0	1	1
];


duration=0.008;


% levelSequence=[23.3987	13.3987	18.9769	25.9114	18.6243	16.4689	23.7829	20.4015	19.3156	25.1232	21.4541	24.353	22.1732	24.83	21.8257	19.3798	17.885	26.0332	23.3135	20.7889	20.1446	22.3504	25.7491	25.8517	23.2701	25.2773	24.8727	26.416	23.9425	21.4652	25.3512	18.3484	26.0988	23.2303	17.0335	26.7833	17.1426	21.1627	21.8802	24.8192	18.2626
% ];    
% responseSequence=[1	0	0	1	0	0	1	0	0	1	0	1	1	1	0	0	0	1	1	0	0	1	1	1	1	1	1	1	0	0	1	0	1	1	0	1	0	0	0	1	0
% ];
% duration=0.016;


% levelSequence=[24.9709	14.9709	18.3185	20.893	19.1644	17.8633	12.8493	20.2354	19.4386	14.4227	14.3361	12.0587	18.6425	20.1803	18.1533	20.8734	19.0579	20.4904	18.3061	14.8561	16.6948	19.2928	17.8098	17.6181	21.7709	16.6804	15.2661	17.967	20.0085	20.7096	15.7485	20.0475	16.0825	11.8798	19.4383	13.838	18.3323	15.8729	12.066	15.4367	12.4179
% ];   
% responseSequence=[1	0	1	1	1	1	0	1	1	0	0	0	1	1	0	1	1	1	1	0	1	1	1	0	1	1	1	1	1	1	0	1	0	0	1	0	1	1	0	1	0
% ];
% duration=0.512;
% 
% 
% levelSequence=[20.9248	10.9248	17.6621	16.5301	12.8652	8.4223	7.2398	13.566	14.5773	15.7202	13.4907	15.5048	10.9308	11.9566	13.3548	18.6419	15.8288	13.6961	12.8613	11.1181	10.7736	14.6912	11.3496	18.4924	16.9867	14.6383	17.4752	19.0135	14.7333	14.0739	14.637	11.7261	10.1013	17.7782	18.4642	14.1431	18.7764	19.312	10.9071	17.2639	17.9951
% ];   
% responseSequence=[1	0	1	1	1	0	0	1	0	0	0	1	0	0	0	1	0	1	0	0	0	0	0	1	1	1	1	1	1	0	0	0	0	1	0	1	1	1	0	1	1
% ];
% duration=0.256;


functionType='logistic';
functionType='rareEvent';

[R C]=size(levelSequence);
levelSequence=reshape(levelSequence',1,R*C);

[R C]=size(responseSequence);
responseSequence=reshape(responseSequence',1,R*C);

figure(2),clf
figure(3),clf
figure(2), plot(levelSequence,responseSequence,'k.'),hold on
set(gca,'ytick',[0 1])
set(gca,'yticklabel',{'0'; '1'})
ylim([-.5 1.5])
figure(3), plot(levelSequence,responseSequence,'k.'),hold on

figure(2), xlabel('dB SPL'), ylabel('p(yes)')
figure(3), xlabel('dB SPL'), ylabel('p(yes)')

set(gca,'ytick',[0 1])
set(gca,'yticklabel',{'0'; '1'})
ylim([-.5 1.5])
binWidth=1;
minBinCenter=-10;
maxBinCenter=100;
binCenters=minBinCenter:binWidth:maxBinCenter;
noPointers=find(responseSequence==0);
noLevels=levelSequence(noPointers);
yesPointers=find(responseSequence==1);
yesLevels=levelSequence(yesPointers);
noHist=hist(noLevels,binCenters);
yesHist=hist(yesLevels,binCenters);
warning off MATLAB:divideByZero
yesHist=hist(yesLevels,binCenters);
proportions=yesHist./(yesHist+noHist);
warning on MATLAB:divideByZero

idx=find(noHist | yesHist);
noHist=noHist(idx);
yesHist=yesHist(idx);
proportions=proportions(idx);
binCenters=binCenters(idx);

for i=1:length(noHist)
    marker=10*(noHist(i)+yesHist(i))/max(noHist);
    figure(2), plot(binCenters(i),proportions(i),'ko','markersize',marker), hold on
    figure(3), plot(binCenters(i),proportions(i),'ko','markersize',marker), hold on
end







rareEventProfile= fitRareEvent2(levelSequence,responseSequence, duration)
figure(2), hold on
plot(rareEventProfile.predictionLevels, rareEventProfile.predictionsRE,'g')
xlim([min(binCenters) max(binCenters)])
figure(3), hold on
plot(rareEventProfile.predictionLevels, rareEventProfile.predictionsRE,'g')
xlim([20 35])



[bestVmin, bestG, smallestEuclid]=...
    rareEvent(levelSequence,responseSequence, duration);

[bestThreshold, bestSlope, smallestEuclid]=...
    logistic(levelSequence,responseSequence);

% -------------------------------------------- rareEvent
function [bestVmin, bestG, smallestEuclid]=...
    rareEvent(levelSequence,responseSequence, duration)
'rare event'
gs=.01:.01:4;

Vmins=1:10:5000;
P=28*10.^(levelSequence/20);
allSmallestEuclid=[ ];
allBestGs=[];

for Vmin=Vmins
    Euclids=[];
    x=[];y=[];
    for G=gs;
        predictions= 1-exp(-duration*(G*P-Vmin));
        idx=find(predictions<0); predictions(idx)=0;
        Euclid=mean((predictions - responseSequence).^2);
        Euclids=[Euclids Euclid];
    end
    
    [smallestEuclid idx]=min(Euclids);
    % disp(num2str([idx smallestEuclid]))
    allSmallestEuclid=[allSmallestEuclid smallestEuclid];
    bestG=gs(idx);
    allBestGs=[allBestGs bestG];
    
end
[smallestEuclid idx2]=min(allSmallestEuclid)
bestVmin=Vmins(idx2);
bestG=allBestGs(idx2);
% [bestVmin bestG smallestEuclid]

% predictions= 1-exp(-duration*(bestG*P-Vmin));
% errors=responseSequence-predictions;
% rho=corr([errors;  levelSequence]');
% disp(['correlation=' num2str(rho(1,2))])


levels=[min(levelSequence):max(levelSequence)];
P=28*10.^(levels/20);
predictions= 1-exp(-duration*(bestG*P-bestVmin));
idx=find(predictions<0); predictions(idx)=0;
figure(2), plot(levels,predictions,'r'), hold off
title(['g=' num2str(bestG) ':  Vmin=' num2str(bestVmin) ':  Euclid= ' num2str(smallestEuclid)])
figure(4), plot(levels(1:end-1), diff(predictions))

xlim([levels(1)-10 levels(end)+10])
ylim([-0.5 1.5])
pause(.1)

% -------------------------------------------- logistic
function [bestThreshold, bestSlope, smallestEuclid]=...
    logistic(levelSequence,responseSequence)
'logistic'
candidateThresholds=-20:1:50;
candidateSlopes=.01:.01:10;

allSmallestEuclid=[ ];
allThresholds=[];

for thisSlope=candidateSlopes;
    Euclids=[];
    for thisThreshold=candidateThresholds;
        predictions=1./(1+exp(-thisSlope.*(levelSequence-thisThreshold)));
        Euclid=mean((predictions - responseSequence).^2);
        Euclids=[Euclids Euclid];
    end
    [smallestEuclid idx]=min(Euclids);
    % disp(num2str([idx smallestEuclid]))
    allSmallestEuclid=[allSmallestEuclid smallestEuclid];
    bestThreshold=candidateThresholds(idx);
    allThresholds=[allThresholds bestThreshold];
end
[smallestEuclid idx2]=min(allSmallestEuclid)
bestSlope=candidateSlopes(idx2);
bestThreshold=allThresholds(idx2);
% [bestThreshold bestSlope smallestEuclid]

% predictions= 1-exp(-duration*(bestG*P-Vmin));
% errors=responseSequence-predictions;
% rho=corr([errors;  levelSequence]');
% disp(['correlation=' num2str(rho(1,2))])


levels=[min(levelSequence):max(levelSequence)];
bestLogistic=1./(1+exp(-bestSlope*(levels-bestThreshold)));
figure(3), plot(levels,bestLogistic,'r'), hold off
title(['k=' num2str(bestSlope) ':  threshold=' num2str(bestThreshold) ':  Euclid= ' num2str(smallestEuclid)])
% figure(2), hold on, plot(levels,bestLogistic,'g'), hold off
figure(4), hold on, plot(levels(1:end-1), diff(bestLogistic), 'r')


xlim([levels(1)-10 levels(end)+10])
ylim([-0.5 1.5])



% --------------------------------------------------- fitRareEvent
function rareEvent=fitRareEvent2(stimulusLevels, responses, duration, gains, Vmins)
% least squares estimate of *rare event* function
% p(event)=gain*levelmPa -Vmin
% 'responses' is a binary vector of subject's decision.
% 'stimulusLevels' are the corresponding signal levesl (values)
% duration is required to compute the expectation of an event occurring
% gains is an optional list of gains to be tried
% Vmins is an optional list of Vmins to be tried

global experiment
if nargin<5
    minVmin=1; maxVmin=100000; nVmins=30;
    Vmins=[0 logspace(log10(minVmin),log10(maxVmin),nVmins)];
    % disp(Vmins)
    nGains=10;
    dGain=1/nGains;
    gains=dGain:dGain:20;
end

nVmins=length(Vmins);
nGains=length(gains);

rareEvent.bestGain=NaN;
rareEvent.bestVMin=NaN;
rareEvent.thresholddB=0;
rareEvent.bestPaMindB=NaN;
rareEvent.predictionLevels=[];
rareEvent.predictionsRE=[];
rareEvent.Euclid=NaN;

if isempty(stimulusLevels), return, end

% expected slope is negative, rareEvent can not be used
% if experiment.psyFunSlope<0
%     return
% end

% NB calculations in microPascals!
stimulusLevelsAsPressure=28 * 10.^(stimulusLevels/20);

allGains=reshape(repmat(gains,nVmins,1), 1, nVmins*nGains);
allVmins=repmat(Vmins, 1, nGains);

for repeat=1:2
    predictions=NaN(1,length(stimulusLevels));
    gainCount=0; VminCount=0;
    Euclid=inf; bestVmin=0; bestGain=0;
    for gain= gains
        gainCount=gainCount+1;
        VminCount=0;
        for Vmin=Vmins
            VminCount=VminCount+1;
            % all levels are simultaneously assessed
            gP_Vmin=gain*stimulusLevelsAsPressure-Vmin;
            idx=(gP_Vmin>0);
            predictions(idx)= 1-exp(-duration*(gP_Vmin(idx)));
            predictions(~idx)=0;

            error=(predictions - responses).^2;
            error=mean(error(~isnan(error)));
            if error<Euclid
                Euclid=error;
                bestVmin=Vmin;
                bestVminCount=VminCount;
                bestGainCount=gainCount;
                bestGain=gain;
            end
        end
    end
    if repeat==1
        if bestVminCount>1 & bestVminCount<nVmins
            minVmin=Vmins(bestVminCount-1); maxVmin=Vmins(bestVminCount+1); nVmins=30;
            Vmins=[0 logspace(log10(minVmin),log10(maxVmin),nVmins)];
        elseif bestVminCount==1
            Vmins=0:Vmins(2)/30:Vmins(2);
        else
            disp('rareEvent: Vmin estimate may be out of range')
        end
    end
    % disp(Vmins)
end

if bestGainCount==1 | bestGainCount==nGains
    disp('gain estimate may be out of range')
end


[rareEvent.Euclid idx]=min(Euclid);
rareEvent.bestGain=bestGain;
rareEvent.bestVMin=round(bestVmin);
rareEvent.thresholdPa=(-log(0.5)/duration + rareEvent.bestVMin)/rareEvent.bestGain;
rareEvent.thresholddB=20*log10(rareEvent.thresholdPa/28);
rareEvent.bestPaMindB=20*log10((rareEvent.bestVMin/rareEvent.bestGain)/28);

predictionLevels= -50:1:120;
rareEvent.predictionLevels=predictionLevels;
rareEvent.predictionsRE=...
    1-exp(-duration*(rareEvent.bestGain*28 * 10.^(predictionLevels/20)-rareEvent.bestVMin));
rareEvent.predictionsRE(rareEvent.predictionsRE<0)=0;

