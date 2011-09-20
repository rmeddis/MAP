function testSynapse(BFlist,paramsName, paramChanges)
% testSynapse tracks the quantity of available transmitter vesicles
%  the computations are single channel using the first frequency
%  in the targetFrequency box of the expGUI.
% For, speed this function uses only probability and HSR fibers.
% This cannot be changed because of the way AN_IHCsynapse is coded.This).

global experiment  method IHCpreSynapseParams
global  AN_IHCsynapseParams  stimulusParameters
savePath=path;
addpath (['..' filesep 'utilities'],['..' filesep 'MAP'])

if nargin<3
    paramChanges=[];
end

figure(6),clf
plotColors='rbgkrbgkrbgkrbgkrbgkrbgk';
set(gcf,'position',[5    32   264   243])

sampleRate=5e4; dt=1/sampleRate;

maskerLevels=-0:10:100;

targetFrequency=BFlist;

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

    AN_spikesOrProbability='probability';
    showPlotsAndDetails=0;

    global savePavailable
   
        MAP1_14(inputSignal, 1/dt, targetFrequency, ...
            paramsName, AN_spikesOrProbability, paramChanges);

    % ignore LSR channels (if any) at the top of the matrix
    qt=savePavailable(end, :);

    synapsedt=dt;
    time=synapsedt:synapsedt:synapsedt*length(qt);

    figure(6)
    qtMatrix=[qtMatrix; qt];
    plot(time,qt,  plotColors(levelCount))
    hold on
    xlim([0 max(time)])
    ylim([0 AN_IHCsynapseParams.M])
end

set(gcf,'name','pre-synaptic available transmitter')
title(['q - available vesicles:' num2str(BFlist) ' Hz'])
legend(strvcat(num2str(maskerLevels')),'location','southeast')
legend boxoff
grid on

figure(88), [c,H]=contour(time, maskerLevels,qtMatrix,1:12); clabel(c, H);
set(gcf,'position',[ 276    31   328   246])
xlabel('time'), ylabel('maskerLevels')
title('contour plot of available transmitter')
grid on
path(savePath);
