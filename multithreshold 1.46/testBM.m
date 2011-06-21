function testBM (BMlocations, paramsName)
% testBM generates input output functions for DRNL model for any number
% of locations.
% Computations are bast on a single channel model (channelBFs=BF)
% peak displacement (peakAmp) is measured.
%  if DRNLParams.useMOC is chosen, the full model is run (slow)
%  otherwise only DRNL is computed.
% Tuning curves are generated based on a range of frequencies reletove to
% the BF of the location.
%

global    DRNLParams

savePath=path;

addpath (['..' filesep 'utilities'],['..' filesep 'MAP'])

levels=-10:10:90;   nLevels=length(levels);
% levels= 50;   nLevels=length(levels);

relativeFrequencies=[0.25    .5   .75  1  1.25 1.5    2];
relativeFrequencies=1;

% refBMdisplacement is the displacement of the BM at threshold
% 1 nm disp at  threshold (9 kHz, Ruggero)
refBMdisplacement= 1e-8;            % adjusted for 10 nm at 1 kHz

toneDuration=.200;
rampDuration=0.01;
silenceDuration=0.01;

sampleRate=30000;

dbstop if error
figure(3), clf
set(gcf,'position',[280   350   327   326])
set(gcf,'name','DRNL - BM')

finalSummary=[];
nBFs=length(BMlocations);
BFno=0; plotCount=0;
for BF=BMlocations
    BFno=BFno+1;
    plotCount=plotCount+nBFs;
    stimulusFrequencies=BF* relativeFrequencies;
    nFrequencies=length(stimulusFrequencies);

    peakAmpBM=zeros(nLevels,nFrequencies);
    peakAmpBMdB=NaN(nLevels,nFrequencies);
    peakEfferent=NaN(nLevels,nFrequencies);
    peakAREfferent=NaN(nLevels,nFrequencies);


    levelNo=0;
    for leveldB=levels
        disp(['level= ' num2str(leveldB)])
        levelNo=levelNo+1;

        freqNo=0;
        for frequency=stimulusFrequencies
            freqNo=freqNo+1;

            % Generate stimuli
            globalStimParams.FS=sampleRate;
            globalStimParams.overallDuration=...
                toneDuration+silenceDuration;  % s
            stim.phases='sin';
            stim.type='tone';
            stim.toneDuration=toneDuration;
            stim.frequencies=frequency;
            stim.amplitudesdB=leveldB;
            stim.beginSilence=silenceDuration;
            stim.rampOnDur=rampDuration;
            % no offset ramp
            stim.rampOffDur=rampDuration;
            doPlot=0;
            inputSignal=stimulusCreate(globalStimParams, stim, doPlot);
            inputSignal=inputSignal(:,1)';

            %% run the model
            MAPparamsName=paramsName;
            AN_spikesOrProbability='probability';
            % spikes are slow but can be used to study MOC using IC units
            AN_spikesOrProbability='spikes';

            global DRNLoutput MOCattenuation ARattenuation
            MAP1_14(inputSignal, sampleRate, BF, ...
                MAPparamsName, AN_spikesOrProbability);

            DRNLresponse=DRNLoutput;
            peakAmp=max(max(...
                DRNLresponse(:, end-round(length(DRNLresponse)/2):end)));
            peakAmpBM(levelNo,freqNo)=peakAmp;
            peakAmpBMdB(levelNo,freqNo)=...
                20*log10(peakAmp/refBMdisplacement);
            peakEfferent(levelNo,freqNo)=min(min(MOCattenuation));
            peakAREfferent(levelNo,freqNo)=min(min(ARattenuation));

        end  % tone frequency
    end  % level

    %% analyses results and plot

    % BM I/O plot (top panel)
    figure(3)
    subplot(3,nBFs,BFno), cla
    plot(levels,peakAmpBMdB, 'linewidth',2)
    hold on, plot(levels, repmat(refBMdisplacement,1,length(levels)))
    hold off
    title(['BF=' num2str(BF,'%5.0f') ' - ' paramsName])
    xlabel('level')
    %     set(gca,'xtick',levels),  grid on
    if length(levels)>1,xlim([min(levels) max(levels)]), end
    ylabel(['dB re:' num2str(refBMdisplacement,'%6.1e') 'm'])
    ylim([-20 50])
    set(gca,'ytick',[-10 0 10 20 40])
    grid on
    %     legend({num2str(stimulusFrequencies')}, 'location', 'EastOutside')
    UTIL_printTabTable([levels' peakAmpBMdB], ...
        num2str([0 stimulusFrequencies]','%5.0f'), '%5.0f')
    finalSummary=[finalSummary peakAmpBMdB];

    % Tuning curve
    if length(relativeFrequencies)>2
        figure(3), subplot(3,nBFs, 2*nBFs+BFno)
        %         contour(stimulusFrequencies,levels,peakAmpBM,...
        %             [refBMdisplacement refBMdisplacement],'r')
        contour(stimulusFrequencies,levels,peakAmpBM,...
            refBMdisplacement.*[1 5 10 50 100])
        title(['tuning curve at ' num2str(refBMdisplacement) 'm']);
        ylabel('level (dB) at reference')
        xlim([100 10000])
        grid on
        hold on
        set(gca,'xscale','log')
    end


    % MOC contribution
    figure(3)
    subplot(3,nBFs,nBFs+BFno), cla
    plot(levels,20*log10(peakEfferent), 'linewidth',2)
    ylabel('MOC (dB attenuation)'), xlabel('level')
    title(['peak MOC: model= ' AN_spikesOrProbability])
    grid on
    if length(levels)>1, xlim([min(levels) max(levels)]), end

    % AR contribution
    hold on
    plot(levels,20*log10(peakAREfferent), 'r')
    hold off

end % best frequency

UTIL_showStructureSummary(DRNLParams, 'DRNLParams', 10)

    UTIL_printTabTable([levels' finalSummary], ...
        num2str([0 stimulusFrequencies]','%5.0f'), '%5.0f')

path(savePath);
