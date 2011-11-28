function testPhaseLocking(paramsName, paramChanges)

if nargin<2
    paramChanges=[];
end

if nargin<1
    paramsName='Normal';
end

testFrequencies=[250 500 1000 2000 4000 8000];
levels=0:10:100;

figure(14), clf
set(gcf,'position', [980    36   383   321])
set(gcf,'name', 'phase locking')

allStrengths=zeros(length(testFrequencies), length(levels));
peakVectorStrength=zeros(1,length(testFrequencies));

freqCount=0;
for targetFrequency=testFrequencies;
    %single test
    freqCount=freqCount+1;
    vectorStrength=...
        testAN(targetFrequency,targetFrequency, levels,...
        paramsName, paramChanges);
    allStrengths(freqCount,:)=vectorStrength';
    peakVectorStrength(freqCount)=max(vectorStrength');
end
%% plot results
figure(14)
subplot(2,1,2)
plot(levels,allStrengths, '+')
xlabel('levels')
ylabel('vector strength')
legend (num2str(testFrequencies'),'location','eastOutside')

subplot(2,1,1)
semilogx(testFrequencies,peakVectorStrength)
grid on
title ('peak vector strength')
xlabel('frequency')
ylabel('vector strength')

johnson=[250	0.79
500	0.82
1000	0.8
2000	0.7
4000	0.25
5500	0.05
];
hold on
plot(johnson(:,1),johnson(:,2),'o')
legend({'model','Johnson 80'},'location','eastOutside')
hold off


