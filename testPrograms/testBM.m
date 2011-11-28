function testBM (BFlist, paramsName,...
    relativeFrequencies, AN_spikesOrProbability, paramChanges)
% testBM generates input output functions for DRNL model for any number
% of locations (BFlist).
% Each BF is evaluated using a single channel model
%
% Peak displacement as a function of pure tone level(peakAmp) is displayed
%
% If relative Frequencies is set to a range of values, tuning curves will
%   be computed using these stimulus frequencie.
%
% If AN_spikesOrProbability is set to 'spikes' the full model is run (slow)
%  otherwise efferent activity is based on AN pspiking probabilities.
%
% testBM (1000, 'Normal', 1, 'probability', [])
%   for tuning curves at 6 locations:
% testBM ([250 500 1000 2000 4000 8000], 'Normal', [.5 .75 .9 1 1.1 1.25 1.5], 'probability', [])

global    DRNLParams

if nargin<5, paramChanges=[]; end
if nargin<4, AN_spikesOrProbability='spikes'; end
if nargin==0, BFlist=1000; paramsName='Normal'; 
    relativeFrequencies=1; end

savePath=path;
addpath (['..' filesep 'utilities'],['..' filesep 'MAP'])
tic

levels=-10:10:100;   nLevels=length(levels);
% levels= 50;   nLevels=length(levels);

% refBMdisplacement is the displacement of the BM at threshold
% 1 nm disp at  threshold (9 kHz, Ruggero)
% ? adjust for frequency
refBMdisplacement= 1e-8; % adjusted for 10 nm at 1 kHz 

toneDuration=.5;
% toneDuration=.050;
rampDuration=0.01;
silenceDuration=0.01;

sampleRate=30000;

dbstop if error
figure(3), clf
set(gcf,'position',[280   350   327   326])
set(gcf,'name','DRNL - BM')
pause(0.1)

finalSummary=[];
nBFs=length(BFlist);
BFno=0; plotCount=0;
for BF=BFlist
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

            global DRNLoutput MOCattenuation ARattenuation
            MAP1_14(inputSignal, sampleRate, BF, ...
                MAPparamsName, AN_spikesOrProbability, paramChanges);

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
if length(relativeFrequencies)>2
    maxRows=3;
else
    maxRows=2;
end

    % BM I/O plot (top panel)
    figure(3)
    subplot(maxRows,nBFs,BFno), cla
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
    UTIL_printTabTable([levels' peakAmpBMdB peakAmpBM*1e9], ...
        num2str([0 stimulusFrequencies]','%5.0f'), '%5.0f')
    finalSummary=[finalSummary peakAmpBMdB];

    % Tuning curve
    if length(relativeFrequencies)>2
        figure(3), subplot(maxRows,nBFs, 2*nBFs+BFno)
                contour(stimulusFrequencies,levels,peakAmpBM,...
                    [refBMdisplacement refBMdisplacement],'r','linewidth',4)
                ylim([-10 40])
                
%         contour(stimulusFrequencies,levels,peakAmpBM,...
%             refBMdisplacement.*[1 5 10 50 100])
%         ylim([-10 90])

        title(['tuning curve at ' num2str(refBMdisplacement) 'm']);
        ylabel('level (dB) at reference')
        xlim([100 10000])
        grid on
        hold on
        set(gca,'xscale','log')
    end


    % MOC contribution
    figure(3)
    subplot(maxRows,nBFs,nBFs+BFno), cla
    plot(levels,20*log10(peakEfferent), 'linewidth',2)
    ylabel('MOC (dB attenuation)'), xlabel('level (dB SPL)')
    title(['MOC: (' AN_spikesOrProbability ') duration= ' ...
        num2str(1000*toneDuration,'%5.0f') ' ms'])
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
    diff(finalSummary)
    if ~isempty(paramChanges)
        disp(paramChanges)
    end
    
toc
path(savePath);
