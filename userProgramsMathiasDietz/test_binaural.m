function test_binaural
% test_binaural is a first attempt to produce a binaural model
%  incorporating MSO and IC models.
% The monaural response is computed first for left and right stimuli
%  before using the CN response as input to the binaural MSO model
%  that, in turn, feeds a single cell IC model.
%
% The function has no arguments and everything is set up internally
%
% #1
% Identify the file (in 'MAPparamsName') containing the model parameters
%  the default is 'PL' which uses primary like neurons in the CN to
%  simulate spherical bushy cells
%
% #2
% Set AN_spikesOrProbability'). to 'spikes'
%
% #3
% Choose between a tone signal or file input (in 'signalType')
%
% #4
% Set the signal rms level (in leveldBSPL)
%
% #5
% Identify the channels in terms of their best frequencies in the vector
%  BFlist. This is currently a single-channel model, so only one BF needed
%
% #6
% Last minute changes to the parameters fetched earlier can be made using
%  the cell array of strings 'paramChanges'.
%  Each string must have the same format as the corresponding line in the
%  file identified in 'MAPparamsName'
% Currently this is used to specify that only HSR fibers are used and
%  for changing the current per AN spike at the CN dendrite
%
% #7
% specify the parameters of the MSO cells in the MSOParams structure
%
% #8
% specify the parameters of the IC cells in the ICMSOParams structure
%
% #9
% identify the plots required from MAP1_14 (i.e. before the bonaural model)
%
% #10
% Specify ITDs. The program cycles through different ITDs
%

global CNoutput dtSpikes
dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])

%%  #1 monaural model parameter file name
MAPparamsName='PL';


%% #2 'spikes' are mandatory for this model
AN_spikesOrProbability='spikes';


%% #3 pure tone, harmonic sequence or speech file input
signalType= 'tones';
sampleRate= 50000;
duration=0.050;                 % seconds
rampDuration=.005;              % raised cosine ramp (seconds)
beginSilence=0.050;
endSilence=0.050;
toneFrequency= 750;            % or a pure tone (Hz)

%   or
% harmonic sequence (Hz)
% F0=210;
% toneFrequency= F0:F0:8000;

%   or
% signalType= 'file';
% fileName='twister_44kHz';

if strcmp(signalType, 'file')
    % needed for labeling plot
    showMapOptions.fileName=fileName;
else
    showMapOptions.fileName=[];
end

%% #4 rms level
leveldBSPL= 70;


%% #5 number of channels in the model
BFlist=toneFrequency;

%% #6 change model parameters

paramChanges={'IHCpreSynapseParams.tauCa=80e-6;',...
    'MacGregorMultiParams.currentPerSpike=0.800e-6;'};

% Parameter changes can be used to change one or more model parameters
%  *after* the MAPparams file has been read

%% #7 MSO parameters
MSOParams.nNeuronsPerBF=	10;   % N neurons per BF
MSOParams.type = 'primary-like cell';
MSOParams.fibersPerNeuron=4;   % N input fibers
MSOParams.dendriteLPfreq=2000;  % dendritic filter
MSOParams.currentPerSpike=0.11e-6; % (A) per spike
MSOParams.currentPerSpike=0.5e-6; % (A) per spike
MSOParams.Cap=4.55e-9;   % cell capacitance (Siemens)
MSOParams.tauM=5e-4;     % membrane time constant (s)
MSOParams.Ek=-0.01;      % K+ eq. potential (V)
MSOParams.dGkSpike=3.64e-5; % K+ cond.shift on spike,S
MSOParams.tauGk=	0.0012; % K+ conductance tau (s)
MSOParams.Th0=	0.01;   % equilibrium threshold (V)
MSOParams.c=	0.01;       % threshold shift on spike, (V)
MSOParams.tauTh=	0.015;  % variable threshold tau
MSOParams.Er=-0.06;      % resting potential (V)
MSOParams.Eb=0.06;       % spike height (V)
MSOParams.debugging=0;        % (special)
MSOParams.wideband=0;         % special for wideband units

%% #8 IC parameters
ICchopperParams.nNeuronsPerBF=	10;   % N neurons per BF
ICchopperParams.type = 'chopper cell';
ICchopperParams.fibersPerNeuron=10;  % N input fibers
ICchopperParams.dendriteLPfreq=50;   % dendritic filter
ICchopperParams.currentPerSpike=50e-9; % *per spike
ICchopperParams.currentPerSpike=100e-9; % *per spike
ICchopperParams.Cap=1.67e-8; % ??cell capacitance (Siemens)
ICchopperParams.tauM=0.002;  % membrane time constant (s)
ICchopperParams.Ek=-0.01;    % K+ eq. potential (V)
ICchopperParams.dGkSpike=1.33e-4; % K+ cond.shift on spike,S
ICchopperParams.tauGk=	0.0005;% K+ conductance tau (s)
ICchopperParams.Th0=	0.01; % equilibrium threshold (V)
ICchopperParams.c=	0;        % threshold shift on spike, (V)
ICchopperParams.tauTh=	0.02; % variable threshold tau
ICchopperParams.Er=-0.06;    % resting potential (V)
ICchopperParams.Eb=0.06;     % spike height (V)
ICchopperParams.PSTHbinWidth=	1e-4;

%% #9 delare 'showMap' options to control graphical output
% this applies to the monaural input model only
showMapOptions.printModelParameters=0;   % prints all parameters
showMapOptions.showModelOutput=1;   % plot all stages if monaural input

%% #10 ITDs
% the program cycles through a range of stimulus ITDs
ITDs=0e-6:100e-6:2000e-6;
% ITDs=0; % single shot

%% Now start computing!
figure(98), clf, set(gcf, 'name', 'binaural demo')
MSOcounts=zeros(1,length(ITDs));
ICcounts=zeros(1,length(ITDs));
ITDcount=0;
for ITD=ITDs
    ITDcount=ITDcount+1;
    delaySamples=round(ITD* sampleRate);
    %% Generate stimuli
    switch signalType
        case 'tones'
            % Create pure tone stimulus
            dt=1/sampleRate; % seconds
            time=dt: dt: duration;
            inputSignal=sum(sin(2*pi*toneFrequency'*time), 1);
            amp=10^(leveldBSPL/20)*28e-6;   % converts to Pascals (peak)
            inputSignal=amp*inputSignal;
            % apply ramps
            % catch rampTime error
            if rampDuration>0.5*duration, rampDuration=duration/2; end
            rampTime=dt:dt:rampDuration;
            ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
                ones(1,length(time)-length(rampTime))];
            inputSignal=inputSignal.*ramp;
            ramp=fliplr(ramp);
            inputSignal=inputSignal.*ramp;
            % add silence
            intialSilence= zeros(1,round(beginSilence/dt));
            finalSilence= zeros(1,round(endSilence/dt));
            inputSignal= [intialSilence inputSignal finalSilence];

        case 'file'
            %% file input simple or mixed
            [inputSignal sampleRate]=wavread(fileName);
            dt=1/sampleRate;
            inputSignal=inputSignal(:,1);
            targetRMS=20e-6*10^(leveldBSPL/20);
            rms=(mean(inputSignal.^2))^0.5;
            amp=targetRMS/rms;
            inputSignal=inputSignal*amp;
            intialSilence= zeros(1,round(0.1/dt));
            finalSilence= zeros(1,round(0.2/dt));
            inputSignal= [intialSilence inputSignal' finalSilence];
    end

    %% run the monaural model twice
    t=dt:dt:dt*length(inputSignal);
    for ear={'left','right'}
        figure(98), subplot(4,1,1), colour='b'; hold off
        switch ear{1}
            case 'right'
                inputSignal=circshift(inputSignal', delaySamples)';
                hold on
                colour='r';
        end
        plot(t, inputSignal, colour)
        title ('binaural inputs signals')
        ylabel('Pa'), xlabel('time')
        xlim([0 max(t)])

        % call to monaural model
        MAP1_14(inputSignal, sampleRate, BFlist, ...
            MAPparamsName, AN_spikesOrProbability, paramChanges);

        % the model run is now complete. Now display the results
        UTIL_showMAP(showMapOptions, paramChanges)

        % copy the CN inputspiking pattern to the binaural display
        figure(98)
        switch ear{1}
            case 'left'
                CNoutputLeft=CNoutput;
                subplot(4,2,3)
            case 'right'
                CNoutputRight=CNoutput;
                subplot(4,2,4)
        end
        plotInstructions=[];
        plotInstructions.axes=gca;
        plotInstructions.displaydt=dtSpikes;
        plotInstructions.title= ['CN spikes ' ear{1}];
        plotInstructions.rasterDotSize=2;
        if sum(sum(CNoutput))<100
            plotInstructions.rasterDotSize=3;
        end
        UTIL_plotMatrix(CNoutput, plotInstructions);

    end % left/ right ear

    %%  MSO
    % run MSO model using left and right CN spikes
    MSOspikes=MSO(CNoutputLeft,CNoutputRight, dtSpikes, MSOParams);
    
    sumspikes=sum(sum(MSOspikes));
    disp(['ITD/ MSO spikes count= ' num2str([ITD sumspikes])])
    MSOcounts(ITDcount)=sumspikes;
    figure(98), subplot(4,2, 8), cla, hold off, plot(ITDs*1e6, MSOcounts)

    %% IC chopper
    % run IC model using all MSO spikes as input to a single IC cell
    ICMSOspikes=ICchopper(MSOspikes, dtSpikes, ICchopperParams);
    
    sumspikes=sum(sum(ICMSOspikes));
    disp(['ITD/ ICMSO spikes count= ' num2str([ITD sumspikes])])
    ICcounts(ITDcount)=sumspikes;
    figure(98), subplot(4,2,8),hold on, plot(ITDs*1e6, ICcounts,'r')
    xlabel('ITD'), ylabel(' spike count')
    title('MSO (blue)/ IC (red) spike counts')
    legend({'MSO','IC'})

end % ITDs

path(restorePath)


function MSOspikes=MSO(CNoutputLeft,CNoutputRight, dtSpikes, MSOparams)
%%
[nMSOcells nEpochs]=size(CNoutputLeft);
inputCurrent=zeros(nMSOcells, nEpochs);
MSOmembranePotential=inputCurrent;

MSO_tauM=MSOparams.tauM;
MSO_tauGk=MSOparams.tauGk;
MSO_tauTh=MSOparams.tauTh;
MSO_cap=MSOparams.Cap;
MSO_c=MSOparams.c;
MSO_b=MSOparams.dGkSpike;
MSO_Th0=MSOparams.Th0;
MSO_Ek=MSOparams.Ek;
MSO_Eb= MSOparams.Eb;
MSO_Er=MSOparams.Er;

MSO_E=zeros(nMSOcells,1);
MSO_Gk=zeros(nMSOcells,1);
MSO_Th=MSO_Th0*ones(nMSOcells,1);

% Dendritic filtering, all spikes are replaced by CNalphaFunction functions
MSOdendriteLPfreq= MSOparams.dendriteLPfreq;
MSOcurrentPerSpike=MSOparams.currentPerSpike;
MSOspikeToCurrentTau=1/(2*pi*MSOdendriteLPfreq);
t=dtSpikes:dtSpikes:5*MSOspikeToCurrentTau;
MSO_CNalphaFunction= (MSOcurrentPerSpike / ...
    MSOspikeToCurrentTau)*t.*exp(-t / MSOspikeToCurrentTau);
% show alpha function
% figure(84), subplot(4,2,2), plot(t,MSO_CNalphaFunction)
% title(['LP cutoff ' num2str(MSOdendriteLPfreq)])

% convert CN spikes to post-dendritic current
CN_spikes=CNoutputLeft+CNoutputRight;
for i=1:nMSOcells
    x= conv2(CN_spikes(i,:),  MSO_CNalphaFunction);
    inputCurrent(i,:)=x(1:nEpochs);
end

if MSO_c==0
    % faster computation when threshold is stable (c==0)
    for t=1:nEpochs
        s=MSO_E>MSO_Th0;
        dE = (-MSO_E/MSO_tauM + inputCurrent(:,t)/MSO_cap +...
            (MSO_Gk/MSO_cap).*(MSO_Ek-MSO_E))*dtSpikes;
        dGk=-MSO_Gk*dtSpikes/MSO_tauGk +MSO_b*s;
        MSO_E=MSO_E+dE;
        MSO_Gk=MSO_Gk+dGk;
        MSOmembranePotential(:,t)=MSO_E+s.*(MSO_Eb-MSO_E)+MSO_Er;
    end
else
    %  threshold is changing (MSO_c>0; e.g. bushy cell)
    for t=1:nEpochs
        dE = (-MSO_E/MSO_tauM + ...
            inputCurrent(:,t)/MSO_cap + (MSO_Gk/MSO_cap)...
            .*(MSO_Ek-MSO_E))*dtSpikes;
        MSO_E=MSO_E+dE;
        s=MSO_E>MSO_Th;
        MSOmembranePotential(:,t)=MSO_E+s.*(MSO_Eb-MSO_E)+MSO_Er;
        dGk=-MSO_Gk*dtSpikes/MSO_tauGk +MSO_b*s;
        MSO_Gk=MSO_Gk+dGk;

        % After a spike, the threshold is raised
        % otherwise it settles to its baseline
        dTh=-(MSO_Th-MSO_Th0)*dtSpikes/MSO_tauTh +s*MSO_c;
        MSO_Th=MSO_Th+dTh;
    end
end

figure(98),subplot(4,1,3)
hold off, imagesc(MSOmembranePotential)
title ('MSO (V)')

MSOspikes=MSOmembranePotential> -0.01;
% Remove any spike that is immediately followed by a spike
%  NB 'find' works on columns (whence the transposing)
MSOspikes=MSOspikes';
idx=find(MSOspikes);
idx=idx(1:end-1);
MSOspikes(idx+1)=0;
MSOspikes=MSOspikes';


function ICMSOspikes=ICchopper(ICMSOspikes, dtSpikes, ICMSOParams)
%%
ICMSOspikes=sum(ICMSOspikes);
[nICMSOcells nEpochs]=size(ICMSOspikes);
inputCurrent=zeros(nICMSOcells, nEpochs);
ICMSOmembranePotential=inputCurrent;

ICMSO_tauM=ICMSOParams.tauM;
ICMSO_tauGk=ICMSOParams.tauGk;
ICMSO_tauTh=ICMSOParams.tauTh;
ICMSO_cap=ICMSOParams.Cap;
ICMSO_c=ICMSOParams.c;
ICMSO_b=ICMSOParams.dGkSpike;
ICMSO_Th0=ICMSOParams.Th0;
ICMSO_Ek=ICMSOParams.Ek;
ICMSO_Eb= ICMSOParams.Eb;
ICMSO_Er=ICMSOParams.Er;

ICMSO_E=zeros(nICMSOcells,1);
ICMSO_Gk=zeros(nICMSOcells,1);
ICMSO_Th=ICMSO_Th0*ones(nICMSOcells,1);

% Dendritic filtering, all spikes are replaced by CNalphaFunction functions
ICMSOdendriteLPfreq= ICMSOParams.dendriteLPfreq;
ICMSOcurrentPerSpike=ICMSOParams.currentPerSpike;
ICMSOspikeToCurrentTau=1/(2*pi*ICMSOdendriteLPfreq);
t=dtSpikes:dtSpikes:5*ICMSOspikeToCurrentTau;
ICMSOalphaFunction= (ICMSOcurrentPerSpike / ...
    ICMSOspikeToCurrentTau)*t.*exp(-t / ICMSOspikeToCurrentTau);
% show alpha function
% figure(84), subplot(4,2,5), plot(t,ICMSOalphaFunction)
% title(['IC MSO LP cutoff ' num2str(ICMSOdendriteLPfreq)])

% post-dendritic current
for i=1:nICMSOcells
    x= conv2(1*ICMSOspikes(i,:),  ICMSOalphaFunction);
    inputCurrent(i,:)=x(1:nEpochs);
end

if ICMSO_c==0
    % faster computation when threshold is stable (c==0)
    for t=1:nEpochs
        s=ICMSO_E>ICMSO_Th0;
        dE = (-ICMSO_E/ICMSO_tauM + inputCurrent(:,t)/ICMSO_cap +...
            (ICMSO_Gk/ICMSO_cap).*(ICMSO_Ek-ICMSO_E))*dtSpikes;
        dGk=-ICMSO_Gk*dtSpikes/ICMSO_tauGk +ICMSO_b*s;
        ICMSO_E=ICMSO_E+dE;
        ICMSO_Gk=ICMSO_Gk+dGk;
        ICMSOmembranePotential(:,t)=ICMSO_E+s.*(ICMSO_Eb-ICMSO_E)+ICMSO_Er;
    end
else
    %  threshold is changing (ICMSO_c>0; e.g. bushy cell)
    for t=1:nEpochs
        dE = (-ICMSO_E/ICMSO_tauM + ...
            inputCurrent(:,t)/ICMSO_cap + (ICMSO_Gk/ICMSO_cap)...
            .*(ICMSO_Ek-ICMSO_E))*dtSpikes;
        ICMSO_E=ICMSO_E+dE;
        s=ICMSO_E>ICMSO_Th;
        ICMSOmembranePotential(:,t)=ICMSO_E+s.*(ICMSO_Eb-ICMSO_E)+ICMSO_Er;
        dGk=-ICMSO_Gk*dtSpikes/ICMSO_tauGk +ICMSO_b*s;
        ICMSO_Gk=ICMSO_Gk+dGk;

        % After a spike, the threshold is raised
        % otherwise it settles to its baseline
        dTh=-(ICMSO_Th-ICMSO_Th0)*dtSpikes/ICMSO_tauTh +s*ICMSO_c;
        ICMSO_Th=ICMSO_Th+dTh;
    end
end

t=dtSpikes:dtSpikes:dtSpikes*length(ICMSOmembranePotential);
figure(98),subplot(4,2,7)
plot(t, ICMSOmembranePotential)
ylim([-0.07 0]), xlim([0 max(t)])
title('IC (V)')

ICMSOspikes=ICMSOmembranePotential> -0.01;
% now remove any spike that is immediately followed by a spike
% NB 'find' works on columns (whence the transposing)
ICMSOspikes=ICMSOspikes';
idx=find(ICMSOspikes);
idx=idx(1:end-1);
ICMSOspikes(idx+1)=0;
ICMSOspikes=ICMSOspikes';

