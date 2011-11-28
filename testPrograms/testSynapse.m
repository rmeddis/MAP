function testSynapse(BF,paramsName, AN_spikesOrProbability, paramChanges)
% testSynapse tracks the quantity of available transmitter vesicles
%  the computations are single channel using the first frequency
%  in the targetFrequency box of the expGUI.
% This function uses only probability and HSR fibers.
% testSynapse(1000,'Normal',[])

global experiment  IHCpreSynapseParams
global  AN_IHCsynapseParams  stimulusParameters
global savePavailable saveNavailable

savePath=path;
addpath (['..' filesep 'utilities'],['..' filesep 'MAP'])

if nargin<3
    paramChanges=[];
end

if length(BF)>1
    error('Only one value allowed for BF')
end
% AN_spikesOrProbability='probability';
% AN_spikesOrProbability='spikes';
% showPlotsAndDetails=0;

figure(6),clf
plotColors='rbgkrbgkrbgkrbgkrbgkrbgk';
set(gcf,'position',[5    32   264   243])

sampleRate=5e4; dt=1/sampleRate;

switch AN_spikesOrProbability
    case 'probability'
        maskerLevels=-0:20:100;
    case 'spikes'
        maskerLevels=80;
end

targetFrequency=BF;

silenceDuration=0.015;
maskerDuration=0.1;
recoveryDuration=0.15;
rampDuration=0.004;

maskerTime=dt:dt:maskerDuration;

rampTime=dt:dt:rampDuration;
ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
    ones(1,length(maskerTime)-length(rampTime))];
ramp=ramp.*fliplr(ramp);

initialSilence=zeros(1,round(silenceDuration/dt));
recoverySilence=zeros(1,round(recoveryDuration/dt));

signal=sin(2*pi*targetFrequency'*maskerTime);
signal= ramp.*signal;
signal=[initialSilence signal  recoverySilence];

levelCount=0;
qtMatrix=[];
for leveldB=maskerLevels
    levelCount=levelCount+1;

    amp=28e-6*10^(leveldB/20);
    inputSignal=amp*signal;

    MAP1_14(inputSignal, 1/dt, targetFrequency, ...
        paramsName, AN_spikesOrProbability, paramChanges);

    % ignore LSR channels (if any) at the top of the matrix
    switch AN_spikesOrProbability
        case 'probability'
            qt=savePavailable(end, :);
        case 'spikes'
            qt=saveNavailable;
    end
    synapsedt=dt;
    time=synapsedt:synapsedt:synapsedt*length(qt);

    figure(6)
    qtMatrix=[qtMatrix; qt];
    plot(time,qt,  plotColors(levelCount))
    hold on
    xlim([0 max(time)])
    ylim([0 AN_IHCsynapseParams.M])
    xlabel ('time')
end

set(gcf,'name','pre-synaptic available transmitter')
title(['pre-synaptic transmitter:' num2str(BF) ' Hz'])
ylabel(['q - available vesicles'])
legend(strvcat(num2str(maskerLevels')),'location','southeast')
legend boxoff
grid on

switch AN_spikesOrProbability
    case 'probability'
        figure(88), [c,H]=contour(time, maskerLevels,qtMatrix,1:12); clabel(c, H);
        set(gcf,'position',[ 276    31   328   246])
        xlabel('time'), ylabel('maskerLevels')
        title('contour plot of available transmitter')
        grid on
end
path(savePath);
