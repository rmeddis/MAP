% function MAPtwoToneDemo
%
% Demonstration of two-tone suppression in the AN using a
%   F1 tone suppressed by a 1.5*F1 kHz tone
%

global  ANprobRateOutput DRNLoutput
dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])
figure(5), clf

disp('')
primaryToneFrequency=2000;

probedBs=[-100 20 :20: 90];
nProbeLevels=length(probedBs);

% disp(['F1 F1 level= ' num2str([primaryToneFrequency probedB])])

suppressorLevels=0:10:80;
suppressorLevels=0:10:80;

lowestBF=250; 	highestBF= 8000; 	numChannels=21;
BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));

suppressorFrequencies=BFlist;
BFchannel=find(BFlist==primaryToneFrequency);

sampleRate= 40000; % Hz (higher sample rate needed for BF>8000 Hz)
dt=1/sampleRate; % seconds

duration=.050;		      % seconds
tone2Duration=duration/2; % s
startSilenceDuration=0.010;
startSilence=...
    zeros(1,startSilenceDuration*sampleRate);
suppressorStartsPTR=...
    round((startSilenceDuration+duration/2)*sampleRate);

probedBCount=0;
for probedB=probedBs
    probedBCount=probedBCount+1;
    
    peakResponse= zeros(length(suppressorLevels),length(suppressorFrequencies));
    suppFreqCount=0;
    for suppressorFrequency=suppressorFrequencies
        suppFreqCount=suppFreqCount+1;
        
        suppressorLevelCount=0;
        for suppressorDB=suppressorLevels
            suppressorLevelCount=suppressorLevelCount+1;
            
            % primary + suppressor
            primaryLevelDB=probedB;
            suppressorLevelDB=suppressorDB;
            
            
            % primary BF tone
            time1=dt: dt: duration;
            amp=10^(primaryLevelDB/20)*28e-6;
            inputSignal=amp*sin(2*pi*primaryToneFrequency*time1);
            rampDuration=.005; rampTime=dt:dt:rampDuration;
            ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(time1)-length(rampTime))];
            inputSignal=inputSignal.*ramp;
            inputSignal=inputSignal.*fliplr(ramp);
            
            % suppressor
            time2= dt: dt: tone2Duration;
            % B: tone on
            amp=10^(suppressorLevelDB/20)*28e-6;
            inputSignal2=amp*sin(2*pi*suppressorFrequency*time2);
            rampDuration=.005; rampTime=dt:dt:rampDuration;
            ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(time2)-length(rampTime))];
            inputSignal2=inputSignal2.*ramp;
            inputSignal2=inputSignal2.*fliplr(ramp);
            % A: initial silence (delay to suppressor)
            silence=zeros(1,length(time2));
            inputSignal2=[silence inputSignal2];
            
            % add tone and suppressor components
            inputSignal=inputSignal+inputSignal2;
            
            inputSignal=...
                [startSilence inputSignal ];
            
            % run MAP
            MAPparamsName='Normal';
            AN_spikesOrProbability='probability';
            paramChanges={'IHCpreSynapseParams.tauCa=80e-6;'};
            MAP1_14(inputSignal, sampleRate, BFlist, ...
                MAPparamsName, AN_spikesOrProbability, paramChanges);
            
            response=ANprobRateOutput;
            peakResponse(suppressorLevelCount,suppFreqCount)=...
                mean(response(BFchannel,suppressorStartsPTR:end));
            disp(['F2 level= ', num2str([suppressorFrequency suppressorDB ...
                peakResponse(suppressorLevelCount,suppFreqCount)])])
            
            
            %% make movie
            figure(6), subplot(4,1,2)
            axis tight
            set(gca,'nextplot','replacechildren');
            plot(dt:dt:dt*length(inputSignal), inputSignal, 'k')
            title('probe                 -                        suppressor')
            ylim([-.5 .5])
            
            figure(6), subplot(2,1,2)
            PSTHbinWidth=0.010;
            PSTH= UTIL_PSTHmakerb(...
                response, dt, PSTHbinWidth);
            [nY nX]=size(PSTH);
            time=PSTHbinWidth*(0:nX-1);
            surf(time, BFlist, PSTH)
            zlim([0 500]),
            caxis([0 500])
            shading interp
            colormap(jet)
            set(gca, 'yScale','log')
            xlim([0 max(time)+dt])
            ylim([0 max(BFlist)])
            view([0 90]) % view([-8 68])
            title('probe                 -                        suppressor')
            pause(0.05)
        end
    end
    
    %% plot matrix
    
    peakResponse=peakResponse/peakResponse(1,1);
    
    figure(5), subplot(2,nProbeLevels, probedBCount)
    surf(suppressorFrequencies,suppressorLevels,peakResponse)
    shading interp
    view([0 90])
    set(gca,'Xscale','log')
    xlim([min(suppressorFrequencies) max(suppressorFrequencies)])
    set(gca,'Xtick', [1000  4000],'xticklabel',{'1000', '4000'})
    ylim([min(suppressorLevels) max(suppressorLevels)])
    
    subplot(2,nProbeLevels,nProbeLevels+probedBCount)
    contour(suppressorFrequencies,suppressorLevels,peakResponse, ...
        [.1:.1:.9 1.1] )
    set(gca,'Xscale','log')
    set(gca,'Xtick', [1000  4000],'xticklabel',{'1000', '4000'})
    set(gcf, 'name',['primaryToneFrequency= ' num2str(primaryToneFrequency)])
    title([num2str(probedB) ' dB'])
end

path=restorePath;
