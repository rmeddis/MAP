function showMAP (options)
% defaults
% options.showModelParameters=1;
% options.showModelOutput=1;
% options.printFiringRates=1;
% options.showACF=1;
% options.showEfferent=1;
% options.surfProbability=0;
% options.fileName=[];

dbstop if warning

global dt ANdt saveAN_spikesOrProbability savedBFlist saveMAPparamsName...
    savedInputSignal TMoutput OMEoutput ARattenuation ...
    DRNLoutput IHC_cilia_output IHCrestingCiliaCond IHCrestingV...
    IHCoutput ANprobRateOutput ANoutput savePavailable tauCas  ...
    CNoutput  ICoutput ICmembraneOutput ICfiberTypeRates MOCattenuation
global OMEParams DRNLParams IHC_cilia_RPParams IHCpreSynapseParams
global AN_IHCsynapseParams MacGregorParams MacGregorMultiParams


restorePath=path;
addpath ( ['..' filesep 'utilities'], ['..' filesep 'parameterStore'])

if nargin<1
    options.showModelParameters=1;
    options.showModelOutput=1;
    options.printFiringRates=1;
    options.showACF=0;
    options.showEfferent=1;
    options.surfProbability=0;
    options.fileName=[];
end

if options.showModelParameters
    % Read parameters from MAPparams<***> file in 'parameterStore' folder
    %  and print out all parameters
    cmd=['MAPparams' saveMAPparamsName ...
        '(-1, 1/dt, 1);'];
    eval(cmd);
end

if options.printFiringRates
    %% print summary firing rates
    fprintf('\n\n')
    disp('summary')
    disp(['AR: ' num2str(min(ARattenuation))])
    disp(['MOC: ' num2str(min(min(MOCattenuation)))])
    nANfiberTypes=length(tauCas);
    if strcmp(saveAN_spikesOrProbability, 'spikes')
        nANfibers=size(ANoutput,1);
        nHSRfibers=nANfibers/nANfiberTypes;
        duration=size(TMoutput,2)*dt;
        disp(['AN: ' num2str(sum(sum(ANoutput(end-nHSRfibers+1:end,:)))/...
            (nHSRfibers*duration))])
        
        nCNneurons=size(CNoutput,1);
        nHSRCNneuronss=nCNneurons/nANfiberTypes;
        disp(['CN: ' num2str(sum(sum(CNoutput(end-nHSRCNneuronss+1:end,:)))...
            /(nHSRCNneuronss*duration))])
        disp(['IC: ' num2str(sum(sum(ICoutput))/duration)])
        %         disp(['IC by type: ' num2str(mean(ICfiberTypeRates,2)')])
    else
        disp(['AN: ' num2str(mean(mean(ANprobRateOutput)))])
        [PSTH pointsPerBin]= UTIL_makePSTH(ANprobRateOutput, dt, 0.001);
        disp(['max max AN: ' num2str(max(max(...
          PSTH/pointsPerBin )))])
    end
end


%% figure (99) summarises main model output
if options.showModelOutput
    plotInstructions.figureNo=99;
    signalRMS=mean(savedInputSignal.^2)^0.5;
    signalRMSdb=20*log10(signalRMS/20e-6);
    
    % plot signal (1)
    plotInstructions.displaydt=dt;
    plotInstructions.numPlots=6;
    plotInstructions.subPlotNo=1;
    plotInstructions.title=...
        ['stimulus: ' num2str(signalRMSdb, '%4.0f') ' dB SPL'];
    r=size(savedInputSignal,1);
    if r==1, savedInputSignal=savedInputSignal'; end
    UTIL_plotMatrix(savedInputSignal', plotInstructions);
    
    % stapes (2)
    plotInstructions.subPlotNo=2;
    plotInstructions.title= ['stapes displacement'];
    UTIL_plotMatrix(OMEoutput, plotInstructions);
    
    % DRNL (3)
    plotInstructions.subPlotNo=3;
    plotInstructions.title= ['BM displacement'];
    plotInstructions.yValues= savedBFlist;
    UTIL_plotMatrix(DRNLoutput, plotInstructions);
    
    switch saveAN_spikesOrProbability
        case 'spikes'
            % AN (4)
            plotInstructions.displaydt=ANdt;
            plotInstructions.title='AN';
            plotInstructions.subPlotNo=4;
            plotInstructions.yLabel='BF';
            plotInstructions.yValues= savedBFlist;
            plotInstructions.rasterDotSize=1;
            plotInstructions.plotDivider=1;
            if sum(sum(ANoutput))<100
                plotInstructions.rasterDotSize=3;
            end
            UTIL_plotMatrix(ANoutput, plotInstructions);
            
            % CN (5)
            plotInstructions.displaydt=ANdt;
            plotInstructions.subPlotNo=5;
            plotInstructions.title='CN spikes';
            if sum(sum(CNoutput))<100
                plotInstructions.rasterDotSize=3;
            end
            UTIL_plotMatrix(CNoutput, plotInstructions);
            
            % IC (6)
            plotInstructions.displaydt=ANdt;
            plotInstructions.subPlotNo=6;
            plotInstructions.title='IC';
            if size(ICoutput,1)>3
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
            
        otherwise % probability (4-6)
            plotInstructions.displaydt=dt;
            plotInstructions.numPlots=2;
            plotInstructions.subPlotNo=2;
            plotInstructions.yLabel='BF';
            if nANfiberTypes>1,
                plotInstructions.yLabel='LSR    HSR';
                plotInstructions.plotDivider=1;
            end
            plotInstructions.title='AN - spike probability';
            UTIL_plotMatrix(ANprobRateOutput, plotInstructions);
    end
end

%% surface plot of probability
if options.surfProbability
    figure(97), clf
    % select only HSR fibers at the bottom of the matrix
    ANprobRateOutput= ANprobRateOutput(end-length(savedBFlist)+1:end,:);
    [nY nX]=size(ANprobRateOutput);
    if nY>2
        time=dt*(1:nX);
        surf(time, savedBFlist, ANprobRateOutput)
        shading interp
        set(gca, 'yScale','log')
        xlim([0 max(time)]), ylim([0 max(savedBFlist)]), zlim([0 1000])
        xlabel('time (s)')
        ylabel('best frequency (Hz)')
        zlabel('spike rate')
        view([-20 60])
        if isfield(options, 'fileName')
            title ([options.fileName ':   ' num2str(signalRMSdb,'% 3.0f') ' dB'])
        else
            title ([ num2str(signalRMSdb,'% 3.0f') ' dB'])
        end
        
    end
end


%% plot efferent control values as dB
if options.showEfferent
    plotInstructions=[];
    plotInstructions.figureNo=98;
    figure(98), clf
    plotInstructions.displaydt=dt;
    plotInstructions.numPlots=2;
    plotInstructions.subPlotNo=1;
    plotInstructions.zValuesRange=[ -25 0];
    plotInstructions.title= ['AR strength.  Signal level= ' ...
        num2str(signalRMSdb,'%4.0f') ' dB SPL'];
    UTIL_plotMatrix(20*log10(ARattenuation), plotInstructions);
    
    plotInstructions.subPlotNo=2;
    plotInstructions.yValues= savedBFlist;
    plotInstructions.yLabel= 'BF';
    plotInstructions.title= ['MOC strength'];
    plotInstructions.zValuesRange=[ -25 0];
    subplot(2,1,2)
    % imagesc(MOCattenuation)
    UTIL_plotMatrix(20*log10(MOCattenuation), plotInstructions);
    colorbar
end
    
    %% ACF plot if required
    if options.showACF
        tic
        method.dt=dt;
        method.segmentNo=1;
        method.nonlinCF=savedBFlist;
        
        minPitch=	80; maxPitch=	4000; numPitches=100;    % specify lags
        pitches=10.^ linspace(log10(minPitch), log10(maxPitch),numPitches);
        pitches=fliplr(pitches);
        filteredSACFParams.lags=1./pitches;     % autocorrelation lags vector
        filteredSACFParams.acfTau=	.003;       % time constant of running ACF
        filteredSACFParams.lambda=	0.12;       % slower filter to smooth ACF
        filteredSACFParams.lambda=	0.01;       % slower filter to smooth ACF
        
        filteredSACFParams.plotACFs=0;          % special plot (see code)
        filteredSACFParams.plotFilteredSACF=0;  % 0 plots unfiltered ACFs
        filteredSACFParams.plotMoviePauses=.3;          % special plot (see code)
        
        filteredSACFParams.usePressnitzer=0; % attenuates ACF at  long lags
        filteredSACFParams.lagsProcedure=  'useAllLags';
        % filteredSACFParams.lagsProcedure=  'useBernsteinLagWeights';
        % filteredSACFParams.lagsProcedure=  'omitShortLags';
        filteredSACFParams.criterionForOmittingLags=3;
        filteredSACFParams.plotACFsInterval=200;
        
        if filteredSACFParams.plotACFs
            % plot original waveform on ACF plot
            figure(13), clf
            subplot(4,1,1)
            t=dt*(1:length(savedInputSignal));
            plot(t,savedInputSignal)
            xlim([0 t(end)])
            title(['stimulus: ' num2str(signalRMSdb, '%4.0f') ' dB SPL']);
        end
        
        % plot original waveform on summary/smoothed ACF plot
        figure(96), clf
        subplot(2,1,1)
        t=dt*(1:length(savedInputSignal));
        plot(t,savedInputSignal)
        xlim([0 t(end)])
        title(['stimulus: ' num2str(signalRMSdb, '%4.0f') ' dB SPL']);
        
        
        % compute ACF
        switch saveAN_spikesOrProbability
            case 'probability'
                inputToACF=ANprobRateOutput.^0.5;
            otherwise
                inputToACF=ANoutput;
        end
        
        disp ('computing ACF...')
        [P, BFlist, sacf, boundaryValue] = ...
            filteredSACF(inputToACF, method, filteredSACFParams);
        disp(' ACF done.')
        
        % SACF
        subplot(2,1,2)
        imagesc(P)
        ylabel('periodicities (Hz)')
        xlabel('time (s)')
        title(['running smoothed (root) SACF. ' saveAN_spikesOrProbability ' input'])
        pt=[1 get(gca,'ytick')]; % force top xtick to show
        set(gca,'ytick',pt)
        set(gca,'ytickLabel', round(pitches(pt)))
        tt=get(gca,'xtick');
        set(gca,'xtickLabel', round(100*t(tt))/100)
    end
    
    path(restorePath)
