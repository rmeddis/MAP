function testEfferent(probeFrequency,BFlist, levels, ...
    paramsName,paramChanges)
% generates rate/level functions for AAR and MOC
%
% e.g.
% testEfferent(1000,1000, -10:10:80,'Normal',[]);

global dtSpikes MOCattenuation ANtauCas

tic
dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'], ['..' filesep 'utilities'], ...
    ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'],...
    ['..' filesep 'testPrograms'])

if nargin<5, paramChanges=[]; end
if nargin<4, paramsName='Normal'; end
if nargin<3, levels=-10:10:100; end
if nargin==0,
    probeFrequency=1000;
    probeFrequency=100:100:8000;
        lowestBF=250; 	highestBF= 8000; 	numChannels=21;
    % 21 chs (250-8k)includes BFs at 250 500 1000 2000 4000 8000
    BFlist=round(logspace(log10(lowestBF),log10(highestBF),numChannels));
    keyChannel=round(numChannels/2);
%     BFlist=1000;
end
nLevels=length(levels);

toneDuration=.2;   rampDuration=0.002;   silenceDuration=.02;
localPSTHbinwidth=0.001;


sampleRate=64000; dt=1/sampleRate;

%% delare 'showMap' options to control graphical output
showMapOptions.printModelParameters=0;   % prints all parameters
showMapOptions.showModelOutput=0;       % plot of all stages
showMapOptions.printFiringRates=1;      % prints stage activity levels
showMapOptions.showACF=0;               % shows SACF (probability only)
showMapOptions.showEfferent=1;          % tracks of AR and MOC
showMapOptions.surfAN=0;       % 2D plot of HSR response
showMapOptions.surfSpikes=0;            % 2D plot of spikes histogram
showMapOptions.ICrates=0;               % IC rates by CNtauGk


%% pre-allocate storage
AN_HSRonset=zeros(nLevels,1);
AN_HSRsaturated=zeros(nLevels,1);
AN_LSRonset=zeros(nLevels,1);
AN_LSRsaturated=zeros(nLevels,1);
CNLSRrate=zeros(nLevels,1);
CNHSRsaturated=zeros(nLevels,1);
ICHSRsaturated=zeros(nLevels,1);
ICLSRsaturated=zeros(nLevels,1);
vectorStrength=zeros(nLevels,1);

AR=zeros(nLevels,1);
MOC=zeros(nLevels,1);
maxMOC=[];

%% main computational loop (vary level)
levelNo=0;
for leveldB=levels
    levelNo=levelNo+1;
    amp=28e-6*10^(leveldB/20);
    fprintf('%4.0f\t', leveldB)

    %% generate tone and silences
    time=dt:dt:toneDuration;
    rampTime=dt:dt:rampDuration;
    ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
        ones(1,length(time)-length(rampTime))];
    ramp=ramp.*fliplr(ramp);

    silence=zeros(1,round(silenceDuration/dt));

    inputSignal=amp*sin(2*pi*probeFrequency'*time);
    inputSignal=sum(inputSignal);
    inputSignal= ramp.*inputSignal;
    inputSignal=[silence inputSignal];

    %% run the model
    AN_spikesOrProbability='spikes';
%     nExistingParamChanges=length(paramChanges);
%     paramChanges{nExistingParamChanges+1}=...
%         ['AN_IHCsynapseParams.spikesTargetSampleRate=' ...
%         num2str(spikesSampleRate) ';'];

    MAP1_14(inputSignal, 1/dt, BFlist, ...
        paramsName, AN_spikesOrProbability, paramChanges);
    
maxMOC=[maxMOC min(MOCattenuation(keyChannel,:))];
    UTIL_showMAP(showMapOptions)
    pause(0.1)

    %% Auditory nerve evaluate and display (Fig. 5)
    %LSR (same as HSR if no LSR fibers present)
    nTaus=length(ANtauCas);


end % level

%% MOC atten/ level function
figure(21), subplot(2,1,2)
plot(levels, 20*log10(maxMOC), 'k'), hold off
title(' MOC dB attenuation'), ylabel('dB attenuation')
ylim([-30 0])
figure(21), subplot(2,1,1)
plot(levels, maxMOC, 'k'), hold off
title(' MOC attenuation (scalar)'), ylabel('attenuation (scalar)')
ylim([0 1])

set(gcf,'name','MOC atten/level')

path(restorePath)
toc
