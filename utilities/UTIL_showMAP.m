function UTIL_showMAP (showMapOptions)
% UTIL_showMAP produces summaries of the output from MAP's mmost recent run
%  All MAP_14 outputs are stored in global variables and UTIL_showMAP
%   simply assumes that they are in place. It uses whatever it there.
%
% showMapOptions:
% showMapOptions.printModelParameters=1; % print model parameters
% showMapOptions.showModelOutput=1;      % plot all stages output
% showMapOptions.printFiringRates=1;     % mean activity at all stages
% showMapOptions.showACF=1;              % SACF (probabilities only)
% showMapOptions.showEfferent=1;         % plot of efferent activity
% showMapOptions.surfAN=0;               % AN output surf plot
% showMapOptions.fileName='';            % parameter filename for plot title
% showMapOptions.PSTHbinwidth=0.001      % binwidth for surface plots
% showMapOptions.colorbar=1;             % request color bar if appropriate
% showMapOptions.view=[0 90];
% dbstop if warning

% Discover results left behind by MAP1_14
global  savedParamChanges savedBFlist saveAN_spikesOrProbability ...
    saveMAPparamsName savedInputSignal dt dtSpikes ...
    OMEextEarPressure TMoutput OMEoutput DRNLoutput...
    IHC_cilia_output IHCrestingCiliaCond IHCrestingV...
    IHCoutput ANprobRateOutput ANoutput savePavailable saveNavailable ...
    ANtauCas CNtauGk CNoutput ICoutput ICmembraneOutput ICfiberTypeRates...
    MOCattenuation ARattenuation 

% Find parameters created in MAPparams<name>
global inputStimulusParams OMEParams DRNLParams IHC_cilia_RPParams
global IHCpreSynapseParams  AN_IHCsynapseParams
global MacGregorParams MacGregorMultiParams  filteredSACFParams
global ICrate


restorePath=path;
addpath ( ['..' filesep 'utilities'], ['..' filesep 'parameterStore'])

if nargin<1,  showMapOptions=[]; end
paramChanges=savedParamChanges; 

% defaults (plot staged outputs and print rates only) if options omitted
if ~isfield(showMapOptions,'printModelParameters')
    showMapOptions.printModelParameters=0; end
if ~isfield(showMapOptions,'showModelOutput'),showMapOptions.showModelOutput=0;end
if ~isfield(showMapOptions,'printFiringRates'),showMapOptions.printFiringRates=0;end
if ~isfield(showMapOptions,'surfAN'),showMapOptions.surfAN=0;end
if ~isfield(showMapOptions,'ICrates'),showMapOptions.ICrates=0;end
if ~isfield(showMapOptions,'showACF'),showMapOptions.showACF=0;end
if ~isfield(showMapOptions,'showEfferent'),showMapOptions.showEfferent=0;end
if ~isfield(showMapOptions,'PSTHbinwidth'),showMapOptions.PSTHbinwidth=0.001;end
if ~isfield(showMapOptions,'fileName'),showMapOptions.fileName=[];end
if ~isfield(showMapOptions,'colorbar'),showMapOptions.colorbar=1;end
if ~isfield(showMapOptions,'view'),showMapOptions.view=[-25 80];end
if ~isfield(showMapOptions,'ICrasterPlot'),showMapOptions.ICrasterPlot=0;end

%% send all model parameters to command window
if showMapOptions.printModelParameters
    % Read parameters from MAPparams<***> file in 'parameterStore' folder
    %  and print out all parameters
        cmd=['MAPparams' saveMAPparamsName '(-1, 1/dt, 1, paramChanges);'];
        eval(cmd);
end

% ignore zero elements in signal
signalRMS=mean(savedInputSignal(find(savedInputSignal)).^2)^0.5;
signalRMSdb=20*log10(signalRMS/20e-6);
nANfiberTypes=length(ANtauCas);

%% summarise firing rates in command window
if showMapOptions.printFiringRates
    %% print summary firing rates
    fprintf('\n')
    disp('summary')
    disp(['AR: ' num2str(min(ARattenuation))])
    disp(['MOC: ' num2str(min(min(MOCattenuation)))])
    if strcmp(saveAN_spikesOrProbability, 'spikes')
        nANfibers=size(ANoutput,1);
        nHSRfibers=nANfibers/nANfiberTypes;
        duration=size(TMoutput,2)*dt;
        disp(['AN(HSR): ' num2str(sum(sum(ANoutput(end-nHSRfibers+1:end,:)))/...
            (nHSRfibers*duration))])

        nCNneurons=size(CNoutput,1);
        nHSRCNneuronss=nCNneurons/nANfiberTypes;
        disp(['CN(HSR): ' num2str(sum(sum(CNoutput(end-nHSRCNneuronss+1:end,:)))...
            /(nHSRCNneuronss*duration))])
        nICneurons=size(ICoutput,1);
        nHSRICneurons= round(nICneurons/nANfiberTypes);
        ICrate=sum(sum(ICoutput(end-nHSRICneurons+1:end,:)))/duration/nHSRICneurons;
        disp(['IC(HSR): ' num2str(ICrate)])
        %         disp(['IC by type: ' num2str(mean(ICfiberTypeRates,2)')])
    else
        HSRprobOutput= ANprobRateOutput(end-length(savedBFlist)+1:end,:);
        disp(['AN(HSR): ' num2str(mean(mean(HSRprobOutput)))])
        PSTH= UTIL_PSTHmakerb(HSRprobOutput, dt, 0.001);
        disp(['max max AN: ' num2str(max(max(PSTH)))])
    end
end


%% figure (99): display output from all model stages
if showMapOptions.showModelOutput
    plotInstructions.figureNo=99;

    % plot signal (1)
    plotInstructions.displaydt=dt;
    plotInstructions.numPlots=6;
    plotInstructions.subPlotNo=1;
    plotInstructions.title=...
        ['stimulus (Pa).  ' num2str(signalRMSdb, '%4.0f') ' dB SPL'];
    UTIL_plotMatrix(savedInputSignal, plotInstructions);

    % stapes (2)
    plotInstructions.subPlotNo=2;
    plotInstructions.title= ['stapes displacement (m)'];
    UTIL_plotMatrix(OMEoutput, plotInstructions);

    % DRNL (3)
    plotInstructions.subPlotNo=3;
    plotInstructions.yValues= savedBFlist;
    [r c]=size(DRNLoutput);
    if r>1
        plotInstructions.title= ['BM displacement'];
        UTIL_plotMatrix(abs(DRNLoutput), plotInstructions);
    else
        plotInstructions.title= ['BM displacement. BF=' ...
            num2str(savedBFlist) ' Hz'];
        UTIL_plotMatrix(DRNLoutput, plotInstructions);
    end

    switch saveAN_spikesOrProbability
        case 'spikes'
            % AN (4)
            plotInstructions.displaydt=dtSpikes;
            plotInstructions.title='AN';
            plotInstructions.subPlotNo=4;
            plotInstructions.yLabel='BF';
            plotInstructions.yValues= savedBFlist;
            plotInstructions.rasterDotSize=1;
            if length(ANtauCas)==2
                plotInstructions.plotDivider=1;
            else
                plotInstructions.plotDivider=0;
            end
            if sum(sum(ANoutput))<100
                plotInstructions.rasterDotSize=3;
            end
            UTIL_plotMatrix(ANoutput, plotInstructions);

            % CN (5)
            plotInstructions.displaydt=dtSpikes;
            plotInstructions.subPlotNo=5;
            plotInstructions.title='CN spikes';
            if sum(sum(CNoutput))<100
                plotInstructions.rasterDotSize=3;
            end
            UTIL_plotMatrix(CNoutput, plotInstructions);

            % IC (6)
            plotInstructions.displaydt=dtSpikes;
            plotInstructions.subPlotNo=6;
            plotInstructions.title='Brainstem 2nd level';
            if size(ICoutput,1)>1
                if sum(sum(ICoutput))<100
                    plotInstructions.rasterDotSize=3;
                end
                UTIL_plotMatrix(ICoutput, plotInstructions);
            else
                plotInstructions.title='IC (HSR) membrane potential';
                plotInstructions.displaydt=dt;
                plotInstructions.yLabel='V';
                plotInstructions.zValuesRange= [-.1 0];
                UTIL_plotMatrix(ICmembraneOutput, plotInstructions);
            end

        otherwise % AN rate based on probability of firing
            PSTHbinWidth=0.0005;
            PSTH= UTIL_PSTHmakerb(ANprobRateOutput, dt, PSTHbinWidth);
            plotInstructions=[];
            plotInstructions.displaydt=PSTHbinWidth;
            plotInstructions.plotDivider=0;
            plotInstructions.numPlots=2;
            plotInstructions.subPlotNo=2;
            plotInstructions.yLabel='BF';
            plotInstructions.yValues=savedBFlist;
            plotInstructions.xLabel='time (s)';
            plotInstructions.zValuesRange= [0 300];

            if nANfiberTypes>1,
                plotInstructions.yLabel='LSR                   HSR';
                plotInstructions.plotDivider=1;
            end
            plotInstructions.title='AN - spike rate';
            UTIL_plotMatrix(PSTH, plotInstructions);
            shading interp
            if showMapOptions.colorbar, colorbar('southOutside'), end
    end
    set(gcf,'name','MAP output')
end


%% surface plot of AN response
%   probability
if showMapOptions.surfAN &&...
        strcmp(saveAN_spikesOrProbability,'probability') && ...
        length(savedBFlist)>2
    % select only HSR fibers
    figure(97), clf
    PSTHbinWidth=showMapOptions.PSTHbinwidth;
    PSTH= UTIL_PSTHmakerb(...
        ANprobRateOutput(end-length(savedBFlist)+1:end,:), ...
        dt, PSTHbinWidth);
    [nY nX]=size(PSTH);
    time=PSTHbinWidth*(1:nX);
    surf(time, savedBFlist, PSTH)
    caxis([0 500])
    shading interp
    set(gca, 'yScale','log')
    xlim([0 max(time)])
    ylim([0 max(savedBFlist)])
    zlim([0 500])
    myFontSize=10;
    xlabel('time (s)', 'fontsize',myFontSize)
    ylabel('BF (Hz)', 'fontsize',myFontSize)
    zlabel('spike rate')
    view(showMapOptions.view)
    set(gca,'ytick',[500 1000 2000 8000],'fontSize',myFontSize)
    title (['AN response. Level= ' ...
        num2str(signalRMSdb,'% 3.0f') ' dB'...
        '  binwidth= ' num2str(1000*PSTHbinWidth) ' s'])
    if showMapOptions.colorbar, colorbar('southOutside'), end
end

%% surfAN ('spikes')
if showMapOptions.surfAN ...
        && strcmp(saveAN_spikesOrProbability, 'spikes')...
        && length(savedBFlist)>2
    figure(97), clf
    % combine fibers across channels
    nFibersPerChannel=AN_IHCsynapseParams.numFibers;
    [r nEpochs]=size(ANoutput);
    nChannels=round(r/nFibersPerChannel);
    x=reshape(ANoutput,nChannels,nFibersPerChannel,nEpochs);
    x=squeeze(sum(x,2));
    % select only HSR fibers at the bottom of the matrix
    HSRoutput= x(end-length(savedBFlist)+1:end,:);
    PSTH=HSRoutput;
    PSTHbinWidth=showMapOptions.PSTHbinwidth;
    PSTH=UTIL_makePSTH(HSRoutput, dtSpikes, PSTHbinWidth);
    [nY nX]=size(PSTH);
    time=PSTHbinWidth*(1:nX);
    surf(time, savedBFlist, PSTH)
    shading interp
    set(gca, 'yScale','log')
    xlim([0 max(time)])
    ylim([0 max(savedBFlist)])
    %     zlim([0 1000])
    xlabel('time (s)')
    ylabel('best frequency (Hz)')
    zlabel('spike rate')
    view(showMapOptions.view)
    title ([showMapOptions.fileName ':   ' num2str(signalRMSdb,'% 3.0f') ' dB'])
    set(97,'name', 'spikes surface plot')
end

%% IC raster plot
if showMapOptions.ICrasterPlot &&...
        strcmp(saveAN_spikesOrProbability,'spikes') && ...
        length(savedBFlist)>2
    figure(91), clf
    plotInstructions=[];
    plotInstructions.numPlots=2;
    plotInstructions.subPlotNo=2;
    plotInstructions.title=...
        ['IC raster plot'];
    plotInstructions.figureNo=91;
    plotInstructions.displaydt=dtSpikes;
    plotInstructions.title='Brainstem 2nd level';
    plotInstructions.yLabel='BF';
    plotInstructions.yValues= savedBFlist;

    if size(ICoutput,1)>1
        if sum(sum(ICoutput))<100
            plotInstructions.rasterDotSize=3;
        end
        UTIL_plotMatrix(ICoutput, plotInstructions);
    end

    % plot signal (1)
    plotInstructions.displaydt=dt;
    plotInstructions.subPlotNo=1;
    plotInstructions.title=...
        ['stimulus (Pa).  ' num2str(signalRMSdb, '%4.0f') ' dB SPL'];
    UTIL_plotMatrix(savedInputSignal, plotInstructions);

end

%% figure(98) plot efferent control values as dB
if showMapOptions.showEfferent
    plotInstructions=[];
    plotInstructions.figureNo=98;
    figure(98), clf
    plotInstructions.displaydt=dt;
    plotInstructions.numPlots=4;
    plotInstructions.subPlotNo=1;
    plotInstructions.zValuesRange= [-1 1];
    plotInstructions.title= ['RMS level='...
        num2str(signalRMSdb, '%4.0f') ' dB SPL'];
    UTIL_plotMatrix(savedInputSignal, plotInstructions);


    plotInstructions.subPlotNo=2;
    plotInstructions.zValuesRange=[ -25 0];
    plotInstructions.title= ['AR stapes attenuation (dB); tau='...
        num2str(OMEParams.ARtau, '%4.3f') ' s'];
    UTIL_plotMatrix(20*log10(ARattenuation), plotInstructions);

    % MOCattenuation
    plotInstructions.numPlots=2;
    plotInstructions.subPlotNo=2;
    plotInstructions.yValues= savedBFlist;
    plotInstructions.yLabel= 'dB';
    if strcmp(saveAN_spikesOrProbability,'spikes')
        rate2atten=DRNLParams.rateToAttenuationFactor;
        plotInstructions.title= ['MOC atten; tau=' ...
            num2str(DRNLParams.MOCtau) '; factor='...
            num2str(rate2atten, '%6.4f')];
    else
        rate2atten=DRNLParams.rateToAttenuationFactorProb;
        plotInstructions.title= ['MOC atten; tauProb=' ...
            num2str(DRNLParams.MOCtauProb) '; factor='...
            num2str(rate2atten, '%6.4f')];
    end
    plotInstructions.zValuesRange=[0 -DRNLParams.minMOCattenuationdB+5];
    UTIL_plotMatrix(-20*log10(MOCattenuation), plotInstructions);
    hold on
    [r c]=size(MOCattenuation);
    if r>2 && showMapOptions.colorbar
        colorbar('southOutside')
    end
    set(plotInstructions.figureNo, 'name', 'AR/ MOC')

    binwidth=0.1;
    [PSTH ]=UTIL_PSTHmaker(20*log10(MOCattenuation), dt, binwidth);
    PSTH=PSTH*length(PSTH)/length(MOCattenuation);
    t=binwidth:binwidth:binwidth*length(PSTH);
    fprintf('\n\n')
    %     UTIL_printTabTable([t' PSTH'])
    %     fprintf('\n\n')

end

%% ACF plot
if showMapOptions.showACF
    tic
    if filteredSACFParams.plotACFs
        % plot original waveform on ACF plot
        figure(13), clf
        subplot(4,1,1)
        t=dt*(1:length(savedInputSignal));
        plot(t,savedInputSignal)
        xlim([0 t(end)])
        title(['stimulus: ' num2str(signalRMSdb, '%4.0f') ' dB SPL']);
    end

    % compute ACF
    switch saveAN_spikesOrProbability
        case 'probability'
            inputToACF=ANprobRateOutput(end-length(savedBFlist)+1:end,:);
        otherwise
            inputToACF=ANoutput;
    end

    disp ('computing ACF...')
    [P, BFlist, sacf, boundaryValue] = ...
        filteredSACF(inputToACF, dt, savedBFlist, filteredSACFParams);
    disp(' ACF done.')

    % plot original waveform on summary/smoothed ACF plot
    figure(96), clf
    subplot(2,1,1)
    t=dt*(1:length(savedInputSignal));
    plot(t,savedInputSignal)
    xlim([0 t(end)])
    title(['stimulus: ' num2str(signalRMSdb, '%4.0f') ' dB SPL']);

    % plot SACF
    figure(96)
    subplot(2,1,2)
    imagesc(P)
%     surf(filteredSACFParams.lags, t, P)
    ylabel('periodicities (Hz)')
    xlabel('time (s)')
    title(['running smoothed SACF. ' saveAN_spikesOrProbability ' input'])
    pt=[1 get(gca,'ytick')]; % force top ytick to show
    set(gca,'ytick',pt)
    pitches=1./filteredSACFParams.lags;     % autocorrelation lags vector
    set(gca,'ytickLabel', round(pitches(pt)))
    [nCH nTimes]=size(P);
    t=dt:dt:dt*nTimes;
    tt=get(gca,'xtick');
    set(gca,'xtickLabel', round(100*t(tt))/100)
end

path(restorePath)


%% IC chopper analysis
% global ICrate
% if showMapOptions.ICrates
% [r nEpochs]=size(ICoutput);
% ICrate=zeros(1,length(CNtauGk));
% % convert ICoutput to a 4-D matrix (time, CNtau, BF, fiberType)
% %  NB only one IC unit for any combination.
% y=reshape(ICoutput', ...
%     nEpochs, length(CNtauGk),length(savedBFlist),length(ANtauCas));
% for i=1:length(CNtauGk)
%     ICrate(i)=sum(sum(sum(y(:,i,:,:))))/duration;
%     fprintf('%10.5f\t%6.0f\n', CNtauGk(i), ICrate(i))
% end
% figure(95), plot(CNtauGk,ICrate)
% title ('ICrate'), xlabel('CNtauGk'), ylabel('ICrate')
% end

function ANsmooth = makeANsmooth(ANresponse, sampleRate, winSize, hopSize)
if nargin < 3
    winSize = 25; %default 25 ms window
end
if nargin < 4
    hopSize = 10; %default 10 ms jump between windows
end

winSizeSamples = round(winSize*sampleRate/1000);
hopSizeSamples = round(hopSize*sampleRate/1000);

% smooth
hann = hanning(winSizeSamples);

ANsmooth = [];%Cannot pre-allocate a size as it is unknown until the enframing
for chan = 1:size(ANresponse,1)
    f = enframe(ANresponse(chan,:), hann, hopSizeSamples);
    ANsmooth(chan,:) = mean(f,2)'; %#ok<AGROW> see above comment
end
%         end% ------ OF makeANsmooth
