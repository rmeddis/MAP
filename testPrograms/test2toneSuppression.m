% function test2toneSuppression
%
% Demonstration of two-tone suppression
%
% A probe tone is played at a fixed level and a second tone is introduced
%  half way through the presentation. The response to the combined signal
%  is recorded and analysed.
%
% The second tone called the seeep tone is presented at a reange of
%  frequencies and levels. In a linear system we might expect the addition
%  of a second tone to increase the magnitude of the output. 
%  In fact, it often results in a reduction in the response. 
%  This is called two-tone suppression.
%
% The effect of the sweep tone is represented in a countour plot showing
%  both reductions and increases in response.
% The background colour in this plot is the response to the 
%  fixed tone alone


% Ruggero et al 1992 (Fig 2)
fixedToneFrequency=8600;
fixedToneLevelsdB=[45 :5: 90];
% fixedToneLevelsdB=[30];
fixedToneDuration=.050;		      % seconds
sweeptoneLevelsdB=[ 0 : 10:90];
nSweepToneFrequencies=21;
sampleSweepFrequency=10600;

global  ANprobRateOutput DRNLoutput
dbstop if error
restorePath=path;
addpath (['..' filesep 'MAP'],    ['..' filesep 'wavFileStore'], ...
    ['..' filesep 'utilities'])
figure(5), clf
figure(87), clf

nFixedToneLevels=length(fixedToneLevelsdB);

lowestSweepFrequency=fixedToneFrequency/6;
highestSweepFrequency=fixedToneFrequency*3;
sweepToneFrequencies=round(logspace(log10(lowestSweepFrequency), ...
    log10(highestSweepFrequency), nSweepToneFrequencies));

% key channels for snapshots
[a BFchannel]=min((sweepToneFrequencies-fixedToneFrequency).^2);
[a sampleChannel]=min((sweepToneFrequencies-sampleSweepFrequency).^2);

sampleRate= max(44100, 4*highestSweepFrequency);
dt=1/sampleRate; % seconds

startSilenceDuration=0.010;
startSilence= zeros(1,startSilenceDuration*sampleRate);

sweepStartPTR=...
    round((startSilenceDuration + fixedToneDuration/2)*sampleRate);

BF_BMresponse=zeros(length(sweepToneFrequencies), ...
    length(fixedToneLevelsdB), length(sweeptoneLevelsdB));

fixedTonedBCount=0;
for fixedTonedB=fixedToneLevelsdB
    fixedTonedBCount=fixedTonedBCount+1;

    BMpeakResponse= zeros(length(sweeptoneLevelsdB),length(sweepToneFrequencies));
    ANpeakResponse= zeros(length(sweeptoneLevelsdB),length(sweepToneFrequencies));
    sweepToneLevelCount=0;
    for sweepToneDB=sweeptoneLevelsdB
        sweepToneLevelCount=sweepToneLevelCount+1;
        suppFreqCount=0;
        for sweepToneFrequency=sweepToneFrequencies
            suppFreqCount=suppFreqCount+1;

            % fixedTone tone
            time1=dt: dt: fixedToneDuration;
            amp=10^(fixedTonedB/20)*28e-6;
            inputSignal=amp*sin(2*pi*fixedToneFrequency*time1);
            rampDuration=.005; rampTime=dt:dt:rampDuration;
            ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(time1)-length(rampTime))];
            inputSignal=inputSignal.*ramp;
            inputSignal=inputSignal.*fliplr(ramp);
            nsignalPoints=length(inputSignal);
            sweepStart=round(nsignalPoints/2);
            nSweepPoints=nsignalPoints-sweepStart;

            % sweepTone
            time2= dt: dt: dt*nSweepPoints;
            % B: tone on
            amp=10^(sweepToneDB/20)*28e-6;
            inputSignal2=amp*sin(2*pi*sweepToneFrequency*time2);
            rampDuration=.005; rampTime=dt:dt:rampDuration;
            ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(time2)-length(rampTime))];
            inputSignal2=inputSignal2.*ramp;
            inputSignal2=inputSignal2.*fliplr(ramp);
            % add tone and sweepTone components
            inputSignal(sweepStart+1:end)= inputSignal(sweepStart+1:end)+inputSignal2;

            inputSignal=...
                [startSilence inputSignal ];

            %% run MAP
                MAPparamsName='Normal';
                AN_spikesOrProbability='probability';
                % only use HSR fibers (NB no acoustic reflex)
                paramChanges={'IHCpreSynapseParams.tauCa=80e-6;'};
                paramChanges={'IHCpreSynapseParams.tauCa=80e-6;',...
                                'DRNLParams.g=0;'};
                MAP1_14(inputSignal, sampleRate, fixedToneFrequency, ...
                    MAPparamsName, AN_spikesOrProbability, paramChanges);

                % find toneAlone response
                if sweepToneLevelCount==1 && suppFreqCount==1
                    BMfixedToneAloneRate=...
                        mean(abs(DRNLoutput(sweepStartPTR:end)));
                    ANfixedToneAloneRate=...
                        mean(abs(ANprobRateOutput(sweepStartPTR:end)));
                end
                
                   BF_BMresponse(suppFreqCount,fixedTonedBCount, ...
                       sweepToneLevelCount)=...
                        mean(abs(DRNLoutput(sweepStartPTR:end)));

                BMpeakResponse(sweepToneLevelCount,suppFreqCount)=...
                    mean(abs(DRNLoutput(sweepStartPTR:end)))...
                    /BMfixedToneAloneRate;
                ANpeakResponse(sweepToneLevelCount,suppFreqCount)=...
                    mean(abs(ANprobRateOutput(sweepStartPTR:end)))...
                    /ANfixedToneAloneRate;
                disp(['F2: ', num2str([sweepToneFrequency sweepToneDB ...
                    BMpeakResponse(sweepToneLevelCount,suppFreqCount)...
                    ANpeakResponse(sweepToneLevelCount,suppFreqCount)])...
                    ' dB'])

                figure(5)
                time=dt:dt:dt*length(inputSignal);
                subplot(3,1,1), plot(time, inputSignal)
                title(['stimulus: Suppressor=' ...
                    num2str([sweepToneFrequency, sweepToneDB]) ' Hz/ dB'])

                time=dt:dt:dt*length(DRNLoutput);
                subplot(3,1,2), plot(time, DRNLoutput)
                xlim([0 fixedToneDuration])
                        ylim([0 inf])

                time=dt:dt:dt*length(ANprobRateOutput);
                subplot(3,1,2), plot(time, ANprobRateOutput)
                xlim([0 fixedToneDuration])
                        ylim([0 500])
                title(['ANresponse: fixedTone' num2str([fixedToneFrequency, fixedTonedB]) ' Hz/ dB'])

                subplot(3,2,5)
                contourf(sweepToneFrequencies,sweeptoneLevelsdB,BMpeakResponse, ...
                    [.1:.1:.9 1.05] )
                set(gca,'Xscale','log')
                set(gca,'Xtick', [1000  4000],'xticklabel',{'1000', '4000'})
                set(gcf, 'name',['fixedToneFrequency= ' num2str(fixedToneFrequency)])
                title(['BM' num2str(fixedTonedB) ' dB'])

                subplot(3,2,6)
                contourf(sweepToneFrequencies,sweeptoneLevelsdB,ANpeakResponse, ...
                    [.1:.1:.9 1.05] )
                set(gca,'Xscale','log')
                set(gca,'Xtick', [1000  4000],'xticklabel',{'1000', '4000'})
                set(gcf, 'name',['fixedToneFrequency= ' num2str(fixedToneFrequency)])
                title(['AN:' num2str(fixedTonedB) ' dB'])
                drawnow
        end
    end

    %% plot matrix

    figure (87)
    subplot(3, nFixedToneLevels, fixedTonedBCount)
    surf(sweepToneFrequencies,sweeptoneLevelsdB,BMpeakResponse)
    set(gca,'Xscale','log')
    zlabel('gain')
    xlim([lowestSweepFrequency highestSweepFrequency])
    ylim([min(sweeptoneLevelsdB) max(sweeptoneLevelsdB)])
    title('BM response')
    view([-11 52])

    subplot(3, nFixedToneLevels, nFixedToneLevels+fixedTonedBCount)
    contourf(sweepToneFrequencies,sweeptoneLevelsdB,BMpeakResponse, ...
        [.1:.5:.9 0.99 1.05] )
    hold on
    plot(fixedToneFrequency, fixedTonedB, 'or', 'markerfacecolor','w')
    set(gca,'Xscale','log')
    set(gca,'Xtick', [1000  5000],'xticklabel',{'1000', '5000'})
    ylabel('(BM) sweep level')
    xlabel('(BM) sweep freq')
    title(['fixed tone level=' num2str(fixedTonedB) ' dB'])
%     colorbar

    subplot(3, nFixedToneLevels, 2*nFixedToneLevels+fixedTonedBCount)
    contourf(sweepToneFrequencies,sweeptoneLevelsdB,ANpeakResponse, ...
        [.1:.5:.9 0.99 1.05] )
    hold on
    plot(fixedToneFrequency, fixedTonedB, 'or', 'markerfacecolor','w')
    set(gca,'Xscale','log')
    set(gca,'Xtick', [1000  5000],'xticklabel',{'1000', '5000'})
    ylabel('(AN) sweep level')
    xlabel('(AN) sweep freq')
    title(['fixed tone level=' num2str(fixedTonedB) ' dB'])
%     colorbar
    
end

% Ruggero fig 2 (probe tone level is x-axis, sweep tone level is y-axis
figure(1),semilogy(fixedToneLevelsdB,squeeze(BF_BMresponse(sampleChannel,:,:)))
ylim([-inf inf])
legend(num2str(sweeptoneLevelsdB'),'location','southeast')
xlabel('probe SPL')
ylabel ('displacement (m)')
title(['Probe ' num2str(fixedToneFrequency) ' Hz. Sweep ' ...
    num2str(sweepToneFrequencies(sampleChannel)) ' Hz.'])
path=restorePath;
