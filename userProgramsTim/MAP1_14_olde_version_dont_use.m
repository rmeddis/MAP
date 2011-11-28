function  MAP1_14(inputSignal, sampleRate, BFlist, MAPparamsName, ...
    AN_spikesOrProbability, paramChanges)
% To test this function use test_MAP1_14 in this folder
%
% example:
% <navigate to 'MAP1_14\MAP'>
%  [inputSignal FS] = wavread('../wavFileStore/twister_44kHz');
%  MAP1_14(inputSignal, FS, -1, 'Normal', 'probability', [])
%
% All arguments are mandatory.
%
%  BFlist is a vector of BFs but can be '-1' to allow MAPparams to choose
%  MAPparamsName='Normal';          % source of model parameters
%  AN_spikesOrProbability='spikes'; % or 'probability'
%  paramChanges is a cell array of strings that can be used to make last
%   minute parameter changes, e.g., to simulate OHC loss
%  e.g.  paramChanges{1}= 'DRNLParams.a=0;'; % disable OHCs
%  e.g.  paramchanges={};                    % no changes
% The model parameters are established in the MAPparams<***> file
%  and stored as global

restorePath=path;
addpath (['..' filesep 'parameterStore'])

global OMEParams DRNLParams IHC_cilia_RPParams IHCpreSynapseParams
global AN_IHCsynapseParams MacGregorParams MacGregorMultiParams

% All of the results of this function are stored as global
global dt ANdt  savedBFlist saveAN_spikesOrProbability saveMAPparamsName...
    savedInputSignal OMEextEarPressure TMoutput OMEoutput ARattenuation ...
    DRNLoutput IHC_cilia_output IHCrestingCiliaCond IHCrestingV...
    IHCoutput ANprobRateOutput ANoutput savePavailable ANtauCas  ...
    CNtauGk CNoutput  ICoutput ICmembraneOutput ICfiberTypeRates ...
    MOCattenuation

% Normally only ICoutput(logical spike matrix) or ANprobRateOutput will be
% needed by the user; so the following will suffice
%   global ANdt ICoutput ANprobRateOutput

% Note that sampleRate has not changed from the original function call and
%  ANprobRateOutput is sampled at this rate
% However ANoutput, CNoutput and IC output are stored as logical
%  'spike' matrices using a lower sample rate (see ANdt).

% When AN_spikesOrProbability is set to probability,
%  no spike matrices are computed.
% When AN_spikesOrProbability is set to 'spikes',
%  no probability output is computed

% Efferent control variables are ARattenuation and MOCattenuation
%  These are scalars between 1 (no attenuation) and 0.
%  They are represented with dt=1/sampleRate (not ANdt)
%  They are computed using either AN probability rate output
%   or IC (spikes) output as approrpriate.
% AR is computed using across channel activity
% MOC is computed on a within-channel basis.

if nargin<1
    error(' MAP1_14 is not a script but a function that must be called')
end

if nargin<6
    paramChanges=[];
end
% Read parameters from MAPparams<***> file in 'parameterStore' folder
% Beware, 'BFlist=-1' is a legitimate argument for MAPparams<>
%  It means that the calling program allows MAPparams to specify the list
cmd=['method=MAPparams' MAPparamsName ...
    '(BFlist, sampleRate, 0, paramChanges);'];
eval(cmd);
BFlist=DRNLParams.nonlinCFs;

% save as global for later plotting if required
savedBFlist=BFlist;
saveAN_spikesOrProbability=AN_spikesOrProbability;
saveMAPparamsName=MAPparamsName;

dt=1/sampleRate;
duration=length(inputSignal)/sampleRate;
% segmentDuration is specified in parameter file (must be >efferent delay)
segmentDuration=method.segmentDuration;
segmentLength=round(segmentDuration/ dt);
segmentTime=dt*(1:segmentLength); % used in debugging plots

% all spiking activity is computed using longer epochs
ANspeedUpFactor=AN_IHCsynapseParams.ANspeedUpFactor;  % e.g.5 times

% inputSignal must be  row vector
[r c]=size(inputSignal);
if r>c, inputSignal=inputSignal'; end       % transpose
% ignore stereo signals
inputSignal=inputSignal(1,:);               % drop any second channel
savedInputSignal=inputSignal;

% Segment the signal
% The sgment length is given but the signal length must be adjusted to be a
% multiple of both the segment length and the reduced segmentlength
[nSignalRows signalLength]=size(inputSignal);
segmentLength=ceil(segmentLength/ANspeedUpFactor)*ANspeedUpFactor;
% Make the signal length a whole multiple of the segment length
nSignalSegments=ceil(signalLength/segmentLength);
padSize=nSignalSegments*segmentLength-signalLength;
pad=zeros(nSignalRows,padSize);
inputSignal=[inputSignal pad];
[ignore signalLength]=size(inputSignal);

% AN (spikes) is computed at a lower sample rate when spikes required
%  so it has a reduced segment length (see 'ANspeeUpFactor' above)
% AN CN and IC all use this sample interval
ANdt=dt*ANspeedUpFactor;
reducedSegmentLength=round(segmentLength/ANspeedUpFactor);
reducedSignalLength= round(signalLength/ANspeedUpFactor);

%% Initialise with respect to each stage before computing
%  by allocating memory,
%  by computing constants
%  by establishing easy to read variable names
% The computations are made in segments and boundary conditions must
%  be established and stored. These are found in variables with
%  'boundary' or 'bndry' in the name

%% OME ---
% external ear resonances
OMEexternalResonanceFilters=OMEParams.externalResonanceFilters;
[nOMEExtFilters c]=size(OMEexternalResonanceFilters);
% details of external (outer ear) resonances
OMEgaindBs=OMEexternalResonanceFilters(:,1);
OMEgainScalars=10.^(OMEgaindBs/20);
OMEfilterOrder=OMEexternalResonanceFilters(:,2);
OMElowerCutOff=OMEexternalResonanceFilters(:,3);
OMEupperCutOff=OMEexternalResonanceFilters(:,4);
% external resonance coefficients
ExtFilter_b=cell(nOMEExtFilters,1);
ExtFilter_a=cell(nOMEExtFilters,1);
for idx=1:nOMEExtFilters
    Nyquist=sampleRate/2;
    [b, a] = butter(OMEfilterOrder(idx), ...
        [OMElowerCutOff(idx) OMEupperCutOff(idx)]...
        /Nyquist);
    ExtFilter_b{idx}=b;
    ExtFilter_a{idx}=a;
end
OMEExtFilterBndry=cell(2,1);
OMEextEarPressure=zeros(1,signalLength); % pressure at tympanic membrane

% pressure to velocity conversion using smoothing filter (50 Hz cutoff)
tau=1/(2*pi*50);
a1=dt/tau-1; a0=1;
b0=1+ a1;
TMdisp_b=b0; TMdisp_a=[a0 a1];
% figure(9), freqz(TMdisp_b, TMdisp_a)
OME_TMdisplacementBndry=[];

% OME high pass (simulates poor low frequency stapes response)
OMEhighPassHighCutOff=OMEParams.OMEstapesLPcutoff;
Nyquist=sampleRate/2;
[stapesDisp_b,stapesDisp_a] = butter(1, OMEhighPassHighCutOff/Nyquist, 'high');
% figure(10), freqz(stapesDisp_b, stapesDisp_a)

OMEhighPassBndry=[];

% OMEampStapes might be reducdant (use OMEParams.stapesScalar)
stapesScalar= OMEParams.stapesScalar;

% Acoustic reflex
efferentDelayPts=round(OMEParams.ARdelay/dt);
% smoothing filter
a1=dt/OMEParams.ARtau-1; a0=1;
b0=1+ a1;
ARfilt_b=b0; ARfilt_a=[a0 a1];

ARattenuation=ones(1,signalLength);
ARrateThreshold=OMEParams.ARrateThreshold; % may not be used
ARrateToAttenuationFactor=OMEParams.rateToAttenuationFactor;
ARrateToAttenuationFactorProb=OMEParams.rateToAttenuationFactorProb;
ARboundary=[];
ARboundaryProb=0;

% save complete OME record (stapes displacement)
OMEoutput=zeros(1,signalLength);
TMoutput=zeros(1,signalLength);

%% BM ---
% BM is represented as a list of locations identified by BF
DRNL_BFs=BFlist;
nBFs= length(DRNL_BFs);

% DRNLchannelParameters=DRNLParams.channelParameters;
DRNLresponse= zeros(nBFs, segmentLength);

MOCrateToAttenuationFactor=DRNLParams.rateToAttenuationFactor;
rateToAttenuationFactorProb=DRNLParams.rateToAttenuationFactorProb;
MOCrateThresholdProb=DRNLParams.MOCrateThresholdProb;

% smoothing filter for MOC
a1=dt/DRNLParams.MOCtau-1; a0=1;
b0=1+ a1;
MOCfilt_b=b0; MOCfilt_a=[a0 a1];
% figure(9), freqz(stapesDisp_b, stapesDisp_a)
MOCboundary=cell(nBFs,1);
MOCprobBoundary=cell(nBFs,1);

MOCattSegment=zeros(nBFs,reducedSegmentLength);
MOCattenuation=ones(nBFs,signalLength);

% if DRNLParams.a>0
%     DRNLcompressionThreshold=10^((1/(1-DRNLParams.c))* ...
%     log10(DRNLParams.b/DRNLParams.a));
% else
%     DRNLcompressionThreshold=inf;
% end
DRNLcompressionThreshold=DRNLParams.cTh;
DRNLlinearOrder= DRNLParams.linOrder;
DRNLnonlinearOrder= DRNLParams.nonlinOrder;

DRNLa=DRNLParams.a;
DRNLb=DRNLParams.b;
DRNLc=DRNLParams.c;
linGAIN=DRNLParams.g;
%
% gammatone filter coefficients for linear pathway
bw=DRNLParams.linBWs';
phi = 2 * pi * bw * dt;
cf=DRNLParams.linCFs';
theta = 2 * pi * cf * dt;
cos_theta = cos(theta);
sin_theta = sin(theta);
alpha = -exp(-phi).* cos_theta;
b0 = ones(nBFs,1);
b1 = 2 * alpha;
b2 = exp(-2 * phi);
z1 = (1 + alpha .* cos_theta) - (alpha .* sin_theta) * i;
z2 = (1 + b1 .* cos_theta) - (b1 .* sin_theta) * i;
z3 = (b2 .* cos(2 * theta)) - (b2 .* sin(2 * theta)) * i;
tf = (z2 + z3) ./ z1;
a0 = abs(tf);
a1 = alpha .* a0;
GTlin_a = [b0, b1, b2];
GTlin_b = [a0, a1];
GTlinOrder=DRNLlinearOrder;
GTlinBdry=cell(nBFs,GTlinOrder);

% nonlinear gammatone filter coefficients 
bw=DRNLParams.nlBWs';
phi = 2 * pi * bw * dt;
cf=DRNLParams.nonlinCFs';
theta = 2 * pi * cf * dt;
cos_theta = cos(theta);
sin_theta = sin(theta);
alpha = -exp(-phi).* cos_theta;
b0 = ones(nBFs,1);
b1 = 2 * alpha;
b2 = exp(-2 * phi);
z1 = (1 + alpha .* cos_theta) - (alpha .* sin_theta) * i;
z2 = (1 + b1 .* cos_theta) - (b1 .* sin_theta) * i;
z3 = (b2 .* cos(2 * theta)) - (b2 .* sin(2 * theta)) * i;
tf = (z2 + z3) ./ z1;
a0 = abs(tf);
a1 = alpha .* a0;
GTnonlin_a = [b0, b1, b2];
GTnonlin_b = [a0, a1];
GTnonlinOrder=DRNLnonlinearOrder;
GTnonlinBdry1=cell(nBFs, GTnonlinOrder);
GTnonlinBdry2=cell(nBFs, GTnonlinOrder);

% complete BM record (BM displacement)
DRNLoutput=zeros(nBFs, signalLength);


%% IHC ---
% IHC cilia activity and receptor potential
% viscous coupling between BM and stereocilia displacement
% Nyquist=sampleRate/2;
% IHCcutoff=1/(2*pi*IHC_cilia_RPParams.tc);
% [IHCciliaFilter_b,IHCciliaFilter_a]=...
%     butter(1, IHCcutoff/Nyquist, 'high');
a1=dt/IHC_cilia_RPParams.tc-1; a0=1;
b0=1+ a1;
% high pass (i.e. low pass reversed)
IHCciliaFilter_b=[a0 a1]; IHCciliaFilter_a=b0;
% figure(9), freqz(IHCciliaFilter_b, IHCciliaFilter_a)

IHCciliaBndry=cell(nBFs,1);

% IHC apical conductance (Boltzman function)
IHC_C= IHC_cilia_RPParams.C;
IHCu0= IHC_cilia_RPParams.u0;
IHCu1= IHC_cilia_RPParams.u1;
IHCs0= IHC_cilia_RPParams.s0;
IHCs1= IHC_cilia_RPParams.s1;
IHCGmax= IHC_cilia_RPParams.Gmax;
IHCGa= IHC_cilia_RPParams.Ga; % (leakage)

IHCGu0 = IHCGa+IHCGmax./(1+exp(IHCu0/IHCs0).*(1+exp(IHCu1/IHCs1)));
IHCrestingCiliaCond=IHCGu0;

% Receptor potential
IHC_Cab= IHC_cilia_RPParams.Cab;
IHC_Gk= IHC_cilia_RPParams.Gk;
IHC_Et= IHC_cilia_RPParams.Et;
IHC_Ek= IHC_cilia_RPParams.Ek;
IHC_Ekp= IHC_Ek+IHC_Et*IHC_cilia_RPParams.Rpc;

IHCrestingV= (IHC_Gk*IHC_Ekp+IHCGu0*IHC_Et)/(IHCGu0+IHC_Gk);

IHC_Vnow= IHCrestingV*ones(nBFs,1); % initial voltage
IHC_RP= zeros(nBFs,segmentLength);

% complete record of IHC receptor potential (V)
IHCciliaDisplacement= zeros(nBFs,segmentLength);
IHCoutput= zeros(nBFs,signalLength);
IHC_cilia_output= zeros(nBFs,signalLength);


%% pre-synapse ---
% Each BF is replicated using a different fiber type to make a 'channel'
% The number of channels is nBFs x nANfiberTypes
% Fiber types are specified in terms of tauCa
nANfiberTypes= length(IHCpreSynapseParams.tauCa);
ANtauCas= IHCpreSynapseParams.tauCa;
nANchannels= nANfiberTypes*nBFs;
synapticCa= zeros(nANchannels,segmentLength);

% Calcium control (more calcium, greater release rate)
ECa=IHCpreSynapseParams.ECa;
gamma=IHCpreSynapseParams.gamma;
beta=IHCpreSynapseParams.beta;
tauM=IHCpreSynapseParams.tauM;
mICa=zeros(nANchannels,segmentLength);
GmaxCa=IHCpreSynapseParams.GmaxCa;
synapse_z= IHCpreSynapseParams.z;
synapse_power=IHCpreSynapseParams.power;

% tauCa vector is established across channels to allow vectorization
%  (one tauCa per channel). Do not confuse with ANtauCas (one pre fiber type)
tauCa=repmat(ANtauCas, nBFs,1);
tauCa=reshape(tauCa, nANchannels, 1);

% presynapse startup values (vectors, length:nANchannels)
% proportion (0 - 1) of Ca channels open at IHCrestingV
mICaCurrent=((1+beta^-1 * exp(-gamma*IHCrestingV))^-1)...
    *ones(nBFs*nANfiberTypes,1);
% corresponding startup currents
ICaCurrent= (GmaxCa*mICaCurrent.^3) * (IHCrestingV-ECa);
CaCurrent= ICaCurrent.*tauCa;

% vesicle release rate at startup (one per channel)
% kt0 is used only at initialisation
kt0= -synapse_z * CaCurrent.^synapse_power;


%% AN ---
% each row of the AN matrices represents one AN fiber
% The results computed either for probabiities *or* for spikes (not both)
% Spikes are necessary if CN and IC are to be computed
nFibersPerChannel= AN_IHCsynapseParams.numFibers;
nANfibers= nANchannels*nFibersPerChannel;
AN_refractory_period= AN_IHCsynapseParams.refractory_period;

y=AN_IHCsynapseParams.y;
l=AN_IHCsynapseParams.l;
x=AN_IHCsynapseParams.x;
r=AN_IHCsynapseParams.r;
M=round(AN_IHCsynapseParams.M);

% probability            (NB initial 'P' on everything)
PAN_ydt = repmat(AN_IHCsynapseParams.y*dt, nANchannels,1);
PAN_ldt = repmat(AN_IHCsynapseParams.l*dt, nANchannels,1);
PAN_xdt = repmat(AN_IHCsynapseParams.x*dt, nANchannels,1);
PAN_rdt = repmat(AN_IHCsynapseParams.r*dt, nANchannels,1);
PAN_rdt_plus_ldt = PAN_rdt + PAN_ldt;
PAN_M=round(AN_IHCsynapseParams.M);

% compute starting values
Pcleft    = kt0* y* M ./ (y*(l+r)+ kt0* l);
Pavailable    = Pcleft*(l+r)./kt0;
Preprocess    = Pcleft*r/x; % canbe fractional

ANprobability=zeros(nANchannels,segmentLength);
ANprobRateOutput=zeros(nANchannels,signalLength);
lengthAbsRefractoryP= round(AN_refractory_period/dt);
cumANnotFireProb=ones(nANchannels,signalLength);
% special variables for monitoring synaptic cleft (specialists only)
savePavailableSeg=zeros(nANchannels,segmentLength);
savePavailable=zeros(nANchannels,signalLength);

% spikes     % !  !  !    ! !        !   !  !
lengthAbsRefractory= round(AN_refractory_period/ANdt);

AN_ydt= repmat(AN_IHCsynapseParams.y*ANdt, nANfibers,1);
AN_ldt= repmat(AN_IHCsynapseParams.l*ANdt, nANfibers,1);
AN_xdt= repmat(AN_IHCsynapseParams.x*ANdt, nANfibers,1);
AN_rdt= repmat(AN_IHCsynapseParams.r*ANdt, nANfibers,1);
AN_rdt_plus_ldt= AN_rdt + AN_ldt;
AN_M= round(AN_IHCsynapseParams.M);

% kt0  is initial release rate
% Establish as a vector (length=channel x number of fibers)
kt0= repmat(kt0', nFibersPerChannel, 1);
kt0=reshape(kt0, nANfibers,1);

% starting values for reservoirs
AN_cleft    = kt0* y* M ./ (y*(l+r)+ kt0* l);
AN_available    = round(AN_cleft*(l+r)./kt0); %must be integer
AN_reprocess    = AN_cleft*r/x;

% output is in a logical array spikes = 1/ 0.
ANspikes= false(nANfibers,reducedSegmentLength);
ANoutput= false(nANfibers,reducedSignalLength);


%% CN (first brain stem nucleus - could be any subdivision of CN)
% Input to a CN neuorn is a random selection of AN fibers within a channel
%  The number of AN fibers used is ANfibersFanInToCN
% CNtauGk (Potassium time constant) determines the rate of firing of
%  the unit when driven hard by a DC input (not normally >350 sp/s)
% If there is more than one value, everything is replicated accordingly

ANavailableFibersPerChan=AN_IHCsynapseParams.numFibers;
ANfibersFanInToCN=MacGregorMultiParams.fibersPerNeuron;

CNtauGk=MacGregorMultiParams.tauGk; % row vector of CN types (by tauGk)
nCNtauGk=length(CNtauGk);

% the total number of 'channels' is now greater
% 'channel' is defined as collections of units with the same parameters
%  i.e. same BF, same ANtau, same CNtauGk
nCNchannels=nANchannels*nCNtauGk;

nCNneuronsPerChannel=MacGregorMultiParams.nNeuronsPerBF;
tauGk=repmat(CNtauGk, nCNneuronsPerChannel,1);
tauGk=reshape(tauGk,nCNneuronsPerChannel*nCNtauGk,1);

% Now the number of neurons has been increased
nCNneurons=nCNneuronsPerChannel*nCNchannels;
CNmembranePotential=zeros(nCNneurons,reducedSegmentLength);

% establish which ANfibers (by name) feed into which CN nuerons
CNinputfiberLists=zeros(nANchannels*nCNneuronsPerChannel, ANfibersFanInToCN);
unitNo=1;
for ch=1:nANchannels
    % Each channel contains a number of units =length(listOfFanInValues)
    for idx=1:nCNneuronsPerChannel
        for idx2=1:nCNtauGk
            fibersUsed=(ch-1)*ANavailableFibersPerChan + ...
                ceil(rand(1,ANfibersFanInToCN)* ANavailableFibersPerChan);
            CNinputfiberLists(unitNo,:)=fibersUsed;
            unitNo=unitNo+1;
        end
    end
end

% input to CN units
AN_PSTH=zeros(nCNneurons,reducedSegmentLength);

% Generate CNalphaFunction function
%  by which spikes are converted to post-synaptic currents
CNdendriteLPfreq= MacGregorMultiParams.dendriteLPfreq;
CNcurrentPerSpike=MacGregorMultiParams.currentPerSpike;
CNspikeToCurrentTau=1/(2*pi*CNdendriteLPfreq);
t=ANdt:ANdt:5*CNspikeToCurrentTau;
CNalphaFunction= (1 / ...
    CNspikeToCurrentTau)*t.*exp(-t /CNspikeToCurrentTau);
CNalphaFunction=CNalphaFunction*CNcurrentPerSpike;

% figure(98), plot(t,CNalphaFunction)
% working memory for implementing convolution

CNcurrentTemp=...
    zeros(nCNneurons,reducedSegmentLength+length(CNalphaFunction)-1);
% trailing alphas are parts of humps carried forward to the next segment
CNtrailingAlphas=zeros(nCNneurons,length(CNalphaFunction));

CN_tauM=MacGregorMultiParams.tauM;
CN_tauTh=MacGregorMultiParams.tauTh;
CN_cap=MacGregorMultiParams.Cap;
CN_c=MacGregorMultiParams.c;
CN_b=MacGregorMultiParams.dGkSpike;
CN_Ek=MacGregorMultiParams.Ek;
CN_Eb= MacGregorMultiParams.Eb;
CN_Er=MacGregorMultiParams.Er;
CN_Th0= MacGregorMultiParams.Th0;
CN_E= zeros(nCNneurons,1);
CN_Gk= zeros(nCNneurons,1);
CN_Th= MacGregorMultiParams.Th0*ones(nCNneurons,1);
CN_Eb=CN_Eb.*ones(nCNneurons,1);
CN_Er=CN_Er.*ones(nCNneurons,1);
CNtimeSinceLastSpike=zeros(nCNneurons,1);
% tauGk is the main distinction between neurons
%  in fact they are all the same in the standard model
tauGk=repmat(tauGk,nANchannels,1);

CNoutput=false(nCNneurons,reducedSignalLength);


%% MacGregor (IC - second nucleus) --------
nICcells=nANchannels*nCNtauGk;  % one cell per channel
CN_PSTH=zeros(nICcells ,reducedSegmentLength);

% ICspikeWidth=0.00015;   % this may need revisiting
% epochsPerSpike=round(ICspikeWidth/ANdt);
% if epochsPerSpike<1
%     error(['MacGregorMulti: sample rate too low to support ' ...
%         num2str(ICspikeWidth*1e6) '  microsec spikes']);
% end

% short names
IC_tauM=MacGregorParams.tauM;
IC_tauGk=MacGregorParams.tauGk;
IC_tauTh=MacGregorParams.tauTh;
IC_cap=MacGregorParams.Cap;
IC_c=MacGregorParams.c;
IC_b=MacGregorParams.dGkSpike;
IC_Th0=MacGregorParams.Th0;
IC_Ek=MacGregorParams.Ek;
IC_Eb= MacGregorParams.Eb;
IC_Er=MacGregorParams.Er;

IC_E=zeros(nICcells,1);
IC_Gk=zeros(nICcells,1);
IC_Th=IC_Th0*ones(nICcells,1);

% Dendritic filtering, all spikes are replaced by CNalphaFunction functions
ICdendriteLPfreq= MacGregorParams.dendriteLPfreq;
ICcurrentPerSpike=MacGregorParams.currentPerSpike;
ICspikeToCurrentTau=1/(2*pi*ICdendriteLPfreq);
t=ANdt:ANdt:3*ICspikeToCurrentTau;
IC_CNalphaFunction= (ICcurrentPerSpike / ...
    ICspikeToCurrentTau)*t.*exp(-t / ICspikeToCurrentTau);
% figure(98), plot(t,IC_CNalphaFunction)

% working space for implementing alpha function
ICcurrentTemp=...
    zeros(nICcells,reducedSegmentLength+length(IC_CNalphaFunction)-1);
ICtrailingAlphas=zeros(nICcells, length(IC_CNalphaFunction));

ICfiberTypeRates=zeros(nANfiberTypes,reducedSignalLength);
ICoutput=false(nICcells,reducedSignalLength);

ICmembranePotential=zeros(nICcells,reducedSegmentLength);
ICmembraneOutput=zeros(nICcells,signalLength);


%% Main program %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%

%  Compute the entire model for each segment
segmentStartPTR=1;
reducedSegmentPTR=1; % when sampling rate is reduced
while segmentStartPTR<signalLength
    segmentEndPTR=segmentStartPTR+segmentLength-1;
    % shorter segments after speed up.
    shorterSegmentEndPTR=reducedSegmentPTR+reducedSegmentLength-1;

    inputPressureSegment=inputSignal...
        (:,segmentStartPTR:segmentStartPTR+segmentLength-1);

    % segment debugging plots
    % figure(98)
    % plot(segmentTime,inputPressureSegment), title('signalSegment')


    % OME ----------------------

    % OME Stage 1: external resonances. Add to inputSignal pressure wave
    y=inputPressureSegment;
    for n=1:nOMEExtFilters
        % any number of resonances can be used
        [x  OMEExtFilterBndry{n}] = ...
            filter(ExtFilter_b{n},ExtFilter_a{n},...
            inputPressureSegment, OMEExtFilterBndry{n});
        x= x* OMEgainScalars(n);
        % This is a parallel resonance so add it
        y=y+x;
    end
    inputPressureSegment=y;
    OMEextEarPressure(segmentStartPTR:segmentEndPTR)= inputPressureSegment;
    
    % OME stage 2: convert input pressure (velocity) to
    %  tympanic membrane(TM) displacement using low pass filter
    [TMdisplacementSegment  OME_TMdisplacementBndry] = ...
        filter(TMdisp_b,TMdisp_a,inputPressureSegment, ...
        OME_TMdisplacementBndry);
    % and save it
    TMoutput(segmentStartPTR:segmentEndPTR)= TMdisplacementSegment;

    % OME stage 3: middle ear high pass effect to simulate stapes inertia
    [stapesDisplacement  OMEhighPassBndry] = ...
        filter(stapesDisp_b,stapesDisp_a,TMdisplacementSegment, ...
        OMEhighPassBndry);

    % OME stage 4:  apply stapes scalar
    stapesDisplacement=stapesDisplacement*stapesScalar;

    % OME stage 5:    acoustic reflex stapes attenuation
    %  Attenuate the TM response using feedback from LSR fiber activity
    if segmentStartPTR>efferentDelayPts
        stapesDisplacement= stapesDisplacement.*...
            ARattenuation(segmentStartPTR-efferentDelayPts:...
            segmentEndPTR-efferentDelayPts);
    end

    % segment debugging plots
    % figure(98)
    % plot(segmentTime, stapesDisplacement), title ('stapesDisplacement')

    % and save
    OMEoutput(segmentStartPTR:segmentEndPTR)= stapesDisplacement;


    %% BM ------------------------------
    % Each location is computed separately
    for BFno=1:nBFs

        %            *linear* path
        linOutput = stapesDisplacement * linGAIN;  % linear gain
       
        for order = 1 : GTlinOrder
            [linOutput GTlinBdry{BFno,order}] = ...
                filter(GTlin_b(BFno,:), GTlin_a(BFno,:), linOutput, ...
                GTlinBdry{BFno,order});
        end

        %           *nonLinear* path
        % efferent attenuation (0 <> 1)
        if segmentStartPTR>efferentDelayPts
            MOC=MOCattenuation(BFno, segmentStartPTR-efferentDelayPts:...
                segmentEndPTR-efferentDelayPts);
        else    % no MOC available yet
            MOC=ones(1, segmentLength);
        end
        % apply MOC to nonlinear input function       
        nonlinOutput=stapesDisplacement.* MOC;

        %       first gammatone filter (nonlin path)
        for order = 1 : GTnonlinOrder
            [nonlinOutput GTnonlinBdry1{BFno,order}] = ...
                filter(GTnonlin_b(BFno,:), GTnonlin_a(BFno,:), ...
                nonlinOutput, GTnonlinBdry1{BFno,order});
        end
        
%         %    original   broken stick instantaneous compression
%         y= nonlinOutput.* DRNLa;  % linear section.
%         % compress parts of the signal above the compression threshold
%         abs_x = abs(nonlinOutput);
%         idx=find(abs_x>DRNLcompressionThreshold);
%         if ~isempty(idx)>0
%             y(idx)=sign(y(idx)).* (DRNLb*abs_x(idx).^DRNLc);
%         end
%         nonlinOutput=y;

        
        %   new broken stick instantaneous compression
        y= nonlinOutput.* DRNLa;  % linear section attenuation/gain.
        % compress parts of the signal above the compression threshold
%         holdY=y;
        abs_y = abs(y);
        idx=find(abs_y>DRNLcompressionThreshold);
        if ~isempty(idx)>0
%             y(idx)=sign(y(idx)).* (DRNLcompressionThreshold + ...
%                 (abs_y(idx)-DRNLcompressionThreshold).^DRNLc);
            y(idx)=sign(y(idx)).* (DRNLcompressionThreshold + ...
                (abs_y(idx)-DRNLcompressionThreshold)*DRNLc);
        end
        nonlinOutput=y;
       
% % Boltzmann compression
% y=(nonlinOutput * DRNLa*100000);
% holdY=y;
% y=abs(y);
% s=10; u=0.0;
% x=1./(1+exp(-(y-u)/s))-0.5;
% nonlinOutput=sign(nonlinOutput).*x/10000;


%         if segmentStartPTR==10*segmentLength+1
%             figure(90)
%         plot(holdY,'b'), hold on
%         plot(nonlinOutput, 'r'), hold off
%         ylim([-1e-5 1e-5])
%         pause(1)
%         end

        %       second filter removes distortion products
        for order = 1 : GTnonlinOrder
            [ nonlinOutput GTnonlinBdry2{BFno,order}] = ...
                filter(GTnonlin_b(BFno,:), GTnonlin_a(BFno,:), ...
                nonlinOutput, GTnonlinBdry2{BFno,order});
        end

        %  combine the two paths to give the DRNL displacement
        DRNLresponse(BFno,:)=linOutput+nonlinOutput;
%         disp(num2str(max(linOutput)))
    end % BF

    % segment debugging plots
    % figure(98)
    %     if size(DRNLresponse,1)>3
    %         imagesc(DRNLresponse)  % matrix display
    %         title('DRNLresponse'); % single or double channel response
    %     else
    %         plot(segmentTime, DRNLresponse)
    %     end

    % and save it
    DRNLoutput(:, segmentStartPTR:segmentEndPTR)= DRNLresponse;


    %% IHC ------------------------------------
    %  BM displacement to IHCciliaDisplacement is  a high-pass filter
    %   because of viscous coupling
    for idx=1:nBFs
        [IHCciliaDisplacement(idx,:)  IHCciliaBndry{idx}] = ...
            filter(IHCciliaFilter_b,IHCciliaFilter_a, ...
            DRNLresponse(idx,:), IHCciliaBndry{idx});
    end
    
    % apply scalar
    IHCciliaDisplacement=IHCciliaDisplacement* IHC_C;

    % and save it
    IHC_cilia_output(:,segmentStartPTR:segmentStartPTR+segmentLength-1)=...
        IHCciliaDisplacement;

    % compute apical conductance
    G=IHCGmax./(1+exp(-(IHCciliaDisplacement-IHCu0)/IHCs0).*...
        (1+exp(-(IHCciliaDisplacement-IHCu1)/IHCs1)));
    Gu=G + IHCGa;

    % Compute receptor potential
    for idx=1:segmentLength
        IHC_Vnow=IHC_Vnow+ (-Gu(:, idx).*(IHC_Vnow-IHC_Et)-...
            IHC_Gk*(IHC_Vnow-IHC_Ekp))*  dt/IHC_Cab;
        IHC_RP(:,idx)=IHC_Vnow;
    end

    % segment debugging plots
    %     if size(IHC_RP,1)>3
    %         surf(IHC_RP), shading interp, title('IHC_RP')
    %     else
    %         plot(segmentTime, IHC_RP)
    %     end

    % and save it
    IHCoutput(:, segmentStartPTR:segmentStartPTR+segmentLength-1)=IHC_RP;


    %% synapse -----------------------------
    % Compute the vesicle release rate for each fiber type at each BF
    % replicate IHC_RP for each fiber type
    Vsynapse=repmat(IHC_RP, nANfiberTypes,1);

    % look-up table of target fraction channels open for a given IHC_RP
    mICaINF=    1./( 1 + exp(-gamma  * Vsynapse)  /beta);
    % fraction of channel open - apply time constant
    for idx=1:segmentLength
        % mICaINF is the current 'target' value of mICa
        mICaCurrent=mICaCurrent+(mICaINF(:,idx)-mICaCurrent)*dt./tauM;
        mICa(:,idx)=mICaCurrent;
    end

    ICa=   (GmaxCa* mICa.^3) .* (Vsynapse- ECa);

    for idx=1:segmentLength
        CaCurrent=CaCurrent +  ICa(:,idx)*dt - CaCurrent*dt./tauCa;
        synapticCa(:,idx)=CaCurrent;
    end
    synapticCa=-synapticCa; % treat synapticCa as positive substance

    % NB vesicleReleaseRate is /s and is independent of dt
    vesicleReleaseRate = synapse_z * synapticCa.^synapse_power; % rate

    % segment debugging plots
    %     if size(vesicleReleaseRate,1)>3
    %         surf(vesicleReleaseRate), shading interp, title('vesicleReleaseRate')
    %     else
    %         plot(segmentTime, vesicleReleaseRate)
    %     end


    %% AN
    switch AN_spikesOrProbability
        case 'probability'
            % No refractory effect is applied
            for t = 1:segmentLength;
                M_Pq=PAN_M-Pavailable;
                M_Pq(M_Pq<0)=0;
                Preplenish = M_Pq .* PAN_ydt;
                Pejected = Pavailable.* vesicleReleaseRate(:,t)*dt;
                Preprocessed = M_Pq.*Preprocess.* PAN_xdt;

                ANprobability(:,t)= min(Pejected,1);
                reuptakeandlost= PAN_rdt_plus_ldt .* Pcleft;
                reuptake= PAN_rdt.* Pcleft;

                Pavailable= Pavailable+ Preplenish- Pejected+ Preprocessed;
                Pcleft= Pcleft + Pejected - reuptakeandlost;
                Preprocess= Preprocess + reuptake - Preprocessed;
                Pavailable(Pavailable<0)=0;
                savePavailableSeg(:,t)=Pavailable;    % synapse tracking
                
            end

            % and save it as *rate*
            ANrate=ANprobability/dt;
            ANprobRateOutput(:, segmentStartPTR:...
                segmentStartPTR+segmentLength-1)=  ANrate;
            % monitor synapse contents (only sometimes used)
            savePavailable(:, segmentStartPTR:segmentStartPTR+segmentLength-1)=...
                savePavailableSeg;
            
            %% Apply refractory effect
               % the probability of a spike's occurring in the preceding
               %  refractory window (t= tnow-refractory period to tnow-)
               %    pFired= 1 - II(1-p(t)),
               % we need a running account of cumProb=II(1-p(t))
               %   cumProb(t)= cumProb(t-1)*(1-p(t))/(1-p(t-refracPeriod))
               %   cumProb(0)=0
               %   pFired(t)= 1-cumProb(t)
               % This gives the fraction of firing events that must be 
               %  discounted because of a firing event in the refractory
               %  period
               %   p(t)= ANprobOutput(t) * pFired(t)
               % where ANprobOutput is the uncorrected firing probability
               %  based on vesicle release rate
               % NB this covers only the absoute refractory period
               % not the relative refractory period. To approximate this it
               % is necessary to extend the refractory period by 50%


                for t = segmentStartPTR:segmentEndPTR;
                    if t>1
                    ANprobRateOutput(:,t)= ANprobRateOutput(:,t)...
                        .* cumANnotFireProb(:,t-1);
                    end
                    % add recent and remove distant probabilities
                    refrac=round(lengthAbsRefractoryP * 1.5);
                    if t>refrac
                        cumANnotFireProb(:,t)= cumANnotFireProb(:,t-1)...
                            .*(1-ANprobRateOutput(:,t)*dt)...
                            ./(1-ANprobRateOutput(:,t-refrac)*dt);
                    end
                end
%                 figure(88), plot(cumANnotFireProb'), title('cumNotFire')
%                 figure(89), plot(ANprobRateOutput'), title('ANprobRateOutput')

            %% Estimate efferent effects. ARattenuation (0 <> 1)
            %  acoustic reflex
            [r c]=size(ANrate);
            if r>nBFs % Only if LSR fibers are computed
                ARAttSeg=mean(ANrate(1:nBFs,:),1); %LSR channels are 1:nBF
                % smooth
                [ARAttSeg, ARboundaryProb] = ...
                    filter(ARfilt_b, ARfilt_a, ARAttSeg, ARboundaryProb);
                ARAttSeg=ARAttSeg-ARrateThreshold;
                ARAttSeg(ARAttSeg<0)=0;   % prevent negative strengths
                ARattenuation(segmentStartPTR:segmentEndPTR)=...
                    (1-ARrateToAttenuationFactorProb.* ARAttSeg);
            end
            %             plot(ARattenuation)

            % MOC attenuation based on within-channel HSR fiber activity
            HSRbegins=nBFs*(nANfiberTypes-1)+1;
            rates=ANrate(HSRbegins:end,:);
            if rateToAttenuationFactorProb<0
                % negative factor implies a fixed attenuation
                MOCattenuation(:,segmentStartPTR:segmentEndPTR)= ...
                    ones(size(rates))* -rateToAttenuationFactorProb;
            else
                for idx=1:nBFs
                    [smoothedRates, MOCprobBoundary{idx}] = ...
                        filter(MOCfilt_b, MOCfilt_a, rates(idx,:), ...
                        MOCprobBoundary{idx});
                    smoothedRates=smoothedRates-MOCrateThresholdProb;
                    smoothedRates=max(smoothedRates, 0);
                    
                    x=(1- smoothedRates* rateToAttenuationFactorProb);
                    x=max(x, 10^(-30/20));
                    MOCattenuation(idx,segmentStartPTR:segmentEndPTR)= x;
                end
            end
            MOCattenuation(MOCattenuation<0)=0.001;

            %             plot(MOCattenuation)


        case 'spikes'
            ANtimeCount=0;
            % implement speed upt
            for t = ANspeedUpFactor:ANspeedUpFactor:segmentLength;
                ANtimeCount=ANtimeCount+1;
                % convert release rate to probabilities
                releaseProb=vesicleReleaseRate(:,t)*ANdt;
                % releaseProb is the release probability per channel
                %  but each channel has many synapses
                releaseProb=repmat(releaseProb',nFibersPerChannel,1);
                releaseProb=reshape(releaseProb, nFibersPerChannel*nANchannels,1);

                % AN_available=round(AN_available); % vesicles must be integer, (?needed)
                M_q=AN_M- AN_available;     % number of missing vesicles
                M_q(M_q<0)= 0;              % cannot be less than 0

                % AN_N1 converts probability to discrete events
                %   it considers each event that might occur
                %   (how many vesicles might be released)
                %   and returns a count of how many were released

                % slow line
%                 probabilities= 1-(1-releaseProb).^AN_available;
                probabilities= 1-intpow((1-releaseProb), AN_available);
                ejected= probabilities> rand(length(AN_available),1);

                reuptakeandlost = AN_rdt_plus_ldt .* AN_cleft;
                reuptake = AN_rdt.* AN_cleft;

                % slow line
%                 probabilities= 1-(1-AN_reprocess.*AN_xdt).^M_q;
                probabilities= 1-intpow((1-AN_reprocess.*AN_xdt), M_q);
                reprocessed= probabilities>rand(length(M_q),1);

                % slow line
%                 probabilities= 1-(1-AN_ydt).^M_q;
                 probabilities= 1-intpow((1-AN_ydt), M_q);

                replenish= probabilities>rand(length(M_q),1);

                AN_available = AN_available + replenish - ejected ...
                    + reprocessed;
                AN_cleft = AN_cleft + ejected - reuptakeandlost;
                AN_reprocess = AN_reprocess + reuptake - reprocessed;

                % ANspikes is logical record of vesicle release events>0
                ANspikes(:, ANtimeCount)= ejected;
            end % t

            % zero any events that are preceded by release events ...
            %  within the refractory period
            % The refractory period consist of two periods
            %  1) the absolute period where no spikes occur
            %  2) a relative period where a spike may occur. This relative
            %     period is realised as a variable length interval
            %     where the length is chosen at random
            %     (esentially a linear ramp up)

            % Andreas has a fix for this
            for t = 1:ANtimeCount-2*lengthAbsRefractory;
                % identify all spikes across fiber array at time (t)
                % idx is a list of channels where spikes occurred
                % ?? try sparse matrices?
                idx=find(ANspikes(:,t));
                for j=idx  % consider each spike
                    % specify variable refractory period
                    %  between abs and 2*abs refractory period
                    nPointsRefractory=lengthAbsRefractory+...
                        round(rand*lengthAbsRefractory);
                    % disable spike potential for refractory period
                    % set all values in this range to 0
                    ANspikes(j,t+1:t+nPointsRefractory)=0;
                end
            end  %t

            % segment debugging
            % plotInstructions.figureNo=98;
            % plotInstructions.displaydt=ANdt;
            %  plotInstructions.numPlots=1;
            %  plotInstructions.subPlotNo=1;
            % UTIL_plotMatrix(ANspikes, plotInstructions);

            % and save it. NB, AN is now on 'speedUp' time
            ANoutput(:, reducedSegmentPTR: shorterSegmentEndPTR)=ANspikes;


            %% CN Macgregor first neucleus -------------------------------
            % input is from AN so ANdt is used throughout
            % Each CNneuron has a unique set of input fibers selected
            %  at random from the available AN fibers (CNinputfiberLists)

            % Create the dendritic current for that neuron
            % First get input spikes to this neuron
            synapseNo=1;
            for ch=1:nCNchannels
                for idx=1:nCNneuronsPerChannel
                    % determine candidate fibers for this unit
                    fibersUsed=CNinputfiberLists(synapseNo,:);
                    % ANpsth has a bin width of ANdt
                    %  (just a simple sum across fibers)
                    AN_PSTH(synapseNo,:) = ...
                        sum(ANspikes(fibersUsed,:), 1);
                    synapseNo=synapseNo+1;
                end
            end

            % One alpha function per spike
            [alphaRows alphaCols]=size(CNtrailingAlphas);

            for unitNo=1:nCNneurons
                CNcurrentTemp(unitNo,:)= ...
                    conv2(AN_PSTH(unitNo,:),CNalphaFunction);
            end
%             disp(['sum(AN_PSTH)= ' num2str(sum(AN_PSTH(1,:)))])
            % add post-synaptic current  left over from previous segment
            CNcurrentTemp(:,1:alphaCols)=...
                CNcurrentTemp(:,1:alphaCols)+ CNtrailingAlphas;

            % take post-synaptic current for this segment
            CNcurrentInput= CNcurrentTemp(:, 1:reducedSegmentLength);
%                 disp(['mean(CNcurrentInput)= ' num2str(mean(CNcurrentInput(1,:)))])

            % trailingalphas are the ends of the alpha functions that
            % spill over into the next segment
            CNtrailingAlphas= ...
                CNcurrentTemp(:, reducedSegmentLength+1:end);

            if CN_c>0
                % variable threshold condition (slow)
                for t=1:reducedSegmentLength
                    CNtimeSinceLastSpike=CNtimeSinceLastSpike-ANdt;
                    s=CN_E>CN_Th & CNtimeSinceLastSpike<0 ;
                    CNtimeSinceLastSpike(s)=0.0005;         % 0.5 ms for sodium spike
                    dE =(-CN_E/CN_tauM + ...
                        CNcurrentInput(:,t)/CN_cap+(...
                        CN_Gk/CN_cap).*(CN_Ek-CN_E))*ANdt;
                    dGk=-CN_Gk*ANdt./tauGk + CN_b*s;
                    dTh=-(CN_Th-CN_Th0)*ANdt/CN_tauTh + CN_c*s;
                    CN_E=CN_E+dE;
                    CN_Gk=CN_Gk+dGk;
                    CN_Th=CN_Th+dTh;
                    CNmembranePotential(:,t)=CN_E+s.*(CN_Eb-CN_E)+CN_Er;
                end
            else
                % static threshold (faster)
                E=zeros(1,reducedSegmentLength);
                Gk=zeros(1,reducedSegmentLength);
                ss=zeros(1,reducedSegmentLength);
                for t=1:reducedSegmentLength
                    % time of previous spike moves back in time
                    CNtimeSinceLastSpike=CNtimeSinceLastSpike-ANdt;
                    % action potential if E>threshold
                    %  allow time for s to reset between events
                    s=CN_E>CN_Th0 & CNtimeSinceLastSpike<0 ;  
                    ss(t)=s(1);
                    CNtimeSinceLastSpike(s)=0.0005; % 0.5 ms for sodium spike
                    dE = (-CN_E/CN_tauM + ...
                        CNcurrentInput(:,t)/CN_cap +...
                        (CN_Gk/CN_cap).*(CN_Ek-CN_E))*ANdt;
                    dGk=-CN_Gk*ANdt./tauGk +CN_b*s;
                    CN_E=CN_E+dE;
                    CN_Gk=CN_Gk+dGk;
                    E(t)=CN_E(1);
                    Gk(t)=CN_Gk(1);
                    % add spike to CN_E and add resting potential (-60 mV)
                    CNmembranePotential(:,t)=CN_E +s.*(CN_Eb-CN_E)+CN_Er;
                end
            end
%             disp(['CN_E= ' num2str(sum(CN_E(1,:)))])
%             disp(['CN_Gk= ' num2str(sum(CN_Gk(1,:)))])
%             disp(['CNmembranePotential= ' num2str(sum(CNmembranePotential(1,:)))])
%             plot(CNmembranePotential(1,:))


            % extract spikes.  A spike is a substantial upswing in voltage
            CN_spikes=CNmembranePotential> -0.02;
%             disp(['CNspikesbefore= ' num2str(sum(sum(CN_spikes)))])

            % now remove any spike that is immediately followed by a spike
            % NB 'find' works on columns (whence the transposing)
            % for each spike put a zero in the next epoch
            CN_spikes=CN_spikes';
            idx=find(CN_spikes);
            idx=idx(1:end-1);
            CN_spikes(idx+1)=0;
            CN_spikes=CN_spikes';
%             disp(['CNspikes= ' num2str(sum(sum(CN_spikes)))])

            % segment debugging
            % plotInstructions.figureNo=98;
            % plotInstructions.displaydt=ANdt;
            %  plotInstructions.numPlots=1;
            %  plotInstructions.subPlotNo=1;
            % UTIL_plotMatrix(CN_spikes, plotInstructions);

            % and save it
            CNoutput(:, reducedSegmentPTR:shorterSegmentEndPTR)=...
                CN_spikes;


            %% IC ----------------------------------------------
                %  MacGregor or some other second order neurons

                % combine CN neurons in same channel (BF x AN tau x CNtau) 
                %  i.e. same BF, same tauCa, same CNtau
                %  to generate inputs to single IC unit
                channelNo=0;
                for idx=1:nCNneuronsPerChannel: ...
                        nCNneurons-nCNneuronsPerChannel+1;
                    channelNo=channelNo+1;
                    CN_PSTH(channelNo,:)=...
                        sum(CN_spikes(idx:idx+nCNneuronsPerChannel-1,:));
                end

                [alphaRows alphaCols]=size(ICtrailingAlphas);
                for ICneuronNo=1:nICcells
                    ICcurrentTemp(ICneuronNo,:)= ...
                        conv2(CN_PSTH(ICneuronNo,:),  IC_CNalphaFunction);
                end

                % add the unused current from the previous convolution
                ICcurrentTemp(:,1:alphaCols)=ICcurrentTemp(:,1:alphaCols)...
                    + ICtrailingAlphas;
                % take what is required and keep the trailing part for next time
                inputCurrent=ICcurrentTemp(:, 1:reducedSegmentLength);
                ICtrailingAlphas=ICcurrentTemp(:, reducedSegmentLength+1:end);

                if IC_c==0
                    % faster computation when threshold is stable (c==0)
                    for t=1:reducedSegmentLength
                        s=IC_E>IC_Th0;
                        dE = (-IC_E/IC_tauM + inputCurrent(:,t)/IC_cap +...
                            (IC_Gk/IC_cap).*(IC_Ek-IC_E))*ANdt;
                        dGk=-IC_Gk*ANdt/IC_tauGk +IC_b*s;
                        IC_E=IC_E+dE;
                        IC_Gk=IC_Gk+dGk;
                        ICmembranePotential(:,t)=IC_E+s.*(IC_Eb-IC_E)+IC_Er;
                    end
                else
                    %  threshold is changing (IC_c>0; e.g. bushy cell)
                    for t=1:reducedSegmentLength
                        dE = (-IC_E/IC_tauM + ...
                            inputCurrent(:,t)/IC_cap + (IC_Gk/IC_cap)...
                            .*(IC_Ek-IC_E))*ANdt;
                        IC_E=IC_E+dE;
                        s=IC_E>IC_Th;
                        ICmembranePotential(:,t)=IC_E+s.*(IC_Eb-IC_E)+IC_Er;
                        dGk=-IC_Gk*ANdt/IC_tauGk +IC_b*s;
                        IC_Gk=IC_Gk+dGk;

                        % After a spike, the threshold is raised
                        % otherwise it settles to its baseline
                        dTh=-(IC_Th-Th0)*ANdt/IC_tauTh +s*IC_c;
                        IC_Th=IC_Th+dTh;
                    end
                end

                ICspikes=ICmembranePotential> -0.01;
                % now remove any spike that is immediately followed by a spike
                % NB 'find' works on columns (whence the transposing)
                ICspikes=ICspikes';
                idx=find(ICspikes);
                idx=idx(1:end-1);
                ICspikes(idx+1)=0;
                ICspikes=ICspikes';

                nCellsPerTau= nICcells/nANfiberTypes;
                firstCell=1;
                lastCell=nCellsPerTau;
                for tauCount=1:nANfiberTypes
                    % separate rates according to fiber types
                    % currently only the last segment is saved
                    ICfiberTypeRates(tauCount, ...
                        reducedSegmentPTR:shorterSegmentEndPTR)=...
                        sum(ICspikes(firstCell:lastCell, :))...
                        /(nCellsPerTau*ANdt);
                    firstCell=firstCell+nCellsPerTau;
                    lastCell=lastCell+nCellsPerTau;
                end
                
                ICoutput(:,reducedSegmentPTR:shorterSegmentEndPTR)=ICspikes;
                
                % store membrane output on original dt scale
                % do this for single channel models only 
                % and only for the HSR-driven IC cells
                if round(nICcells/length(ANtauCas))==1  % single channel
                    % select HSR driven cells
                    x= ICmembranePotential(length(ANtauCas),:);
                    % restore original dt
                    x= repmat(x, ANspeedUpFactor,1);
                    x= reshape(x,1,segmentLength);
                    if nANfiberTypes>1  % save HSR and LSR
                        y=repmat(ICmembranePotential(end,:),...
                            ANspeedUpFactor,1);
                        y= reshape(y,1,segmentLength);
                        x=[x; y];
                    end
                    ICmembraneOutput(:, segmentStartPTR:segmentEndPTR)= x;
                end

                % estimate efferent effects.
                % ARis based on LSR units. LSR channels are 1:nBF
                if nANfiberTypes>1  % AR is multi-channel only
                    ARAttSeg=sum(ICspikes(1:nBFs,:),1)/ANdt;
                    [ARAttSeg, ARboundary] = ...
                        filter(ARfilt_b, ARfilt_a, ARAttSeg, ARboundary);
                    ARAttSeg=ARAttSeg-ARrateThreshold;
                    ARAttSeg(ARAttSeg<0)=0;   % prevent negative strengths
                    % scale up to dt from ANdt
                    x=    repmat(ARAttSeg, ANspeedUpFactor,1);
                    x=reshape(x,1,segmentLength);
                    ARattenuation(segmentStartPTR:segmentEndPTR)=...
                        (1-ARrateToAttenuationFactor* x);
                    ARattenuation(ARattenuation<0)=0.001;
                else
                    % single channel model; disable AR
                    ARattenuation(segmentStartPTR:segmentEndPTR)=...
                        ones(1,segmentLength);
                end

                % MOC attenuation using HSR response only
                % Separate MOC effect for each BF
                HSRbegins=nBFs*(nANfiberTypes-1)+1;
                rates=ICspikes(HSRbegins:end,:)/ANdt;
                for idx=1:nBFs
                    [smoothedRates, MOCboundary{idx}] = ...
                        filter(MOCfilt_b, MOCfilt_a, rates(idx,:), ...
                        MOCboundary{idx});
                    % spont 'rates' is zero for IC
                    MOCattSegment(idx,:)=smoothedRates;
                    % expand timescale back to model dt from ANdt
                    x= repmat(MOCattSegment(idx,:), ANspeedUpFactor,1);
                    x= reshape(x,1,segmentLength);
                    MOCattenuation(idx,segmentStartPTR:segmentEndPTR)= ...
                        (1- MOCrateToAttenuationFactor*  x);
                end
                MOCattenuation(MOCattenuation<0)=0.04;
                % segment debugging
                % plotInstructions.figureNo=98;
                % plotInstructions.displaydt=ANdt;
                %  plotInstructions.numPlots=1;
                %  plotInstructions.subPlotNo=1;
                % UTIL_plotMatrix(ICspikes, plotInstructions);

    end     % AN_spikesOrProbability
    segmentStartPTR=segmentStartPTR+segmentLength;
    reducedSegmentPTR=reducedSegmentPTR+reducedSegmentLength;


end  % segment

%% apply refractory correction to spike probabilities

% the probability of a spike's having occurred in the preceding 
%  refractory window 
%    pFired= 1 - II(1-p(t)), t= tnow-refractory period: tnow-1
% we need a running account of cumProb=II(1-p(t))
%   cumProb(t)= cumProb(t-1)*(1-p(t-1))/(1-p(t-refracPeriod))
%   cumProb(0)=0
%   pFired(t)= 1-cumProb(t)
% whence
%   p(t)= ANprobOutput(t) * pFired(t)
% where ANprobOutput is the uncorrected probability


% switch AN_spikesOrProbability
%     case 'probability'
%         ANprobOutput=ANprobRateOutput*dt;
%         [r nEpochs]=size(ANprobOutput);
%         % find probability of no spikes in refractory period
%         pNoSpikesInRefrac=ones(size(ANprobOutput));
%         pSpike=zeros(size(ANprobOutput));
%         for epochNo=lengthAbsRefractoryP+2:nEpochs
%             pNoSpikesInRefrac(:,epochNo)=...
%                 pNoSpikesInRefrac(:,epochNo-2)...
%                 .*(1-pSpike(:,epochNo-1))...
%                 ./(1-pSpike(:,epochNo-lengthAbsRefractoryP-1));
%             pSpike(:,epochNo)= ANprobOutput(:,epochNo)...
%                 .*pNoSpikesInRefrac(:,epochNo);
%         end
%         ANprobRateOutput=pSpike/dt;
% end

path(restorePath)
