% function MAPtwoToneDemo
% Not for distribution but code worth keeping
% Demonstration of two-tone suppression in the AN using a
%   single channel
%

dbstop if error
% create access to all MAP 1_8 facilities
% addpath ('..\modules', '..\utilities',  '..\parameterStore',  '..\wavFileStore' , '..\testPrograms')
addpath (['..' filesep 'modules'], ['..' filesep 'utilities'],  ['..' filesep 'parameterStore'],  ['..' filesep 'wavFileStore'] , ['..' filesep 'testPrograms'])

moduleSequence= 1:7;  	% up to the AN
figure(3), clf
primaryToneFrequency=2000;
suppressorFrequency=primaryToneFrequency*1.5;

primaryDB=30;
suppressors=20:10:80;
frameCount=0;
for suppressorDB=suppressors
    frameCount=frameCount+1;
    lowestBF=1000; 	highestBF= 5000; 	numChannels=30;
    BFlist=round(logspace(log10(lowestBF), log10(highestBF), numChannels));
    BFlist=2000;
    duration=.020;		      % seconds
    sampleRate= 40000; % Hz (higher sample rate needed for BF>8000 Hz)
    dt=1/sampleRate; % seconds
    
    % for conditionNo=1:3  % probe alone/ suppressor alone/ combined
    for conditionNo=1:3  % probe alone/ suppressor alone/ combined
        switch conditionNo
            case 1
                primaryLevelDB=primaryDB;
                suppressorLevelDB= -100;
                plotGuide.subPlotNo=	1;     % initialize subplot count
                figureTitle=['probe alone:  ' num2str(primaryToneFrequency) ' Hz;   ' num2str(primaryLevelDB) ' dB SPL'];
                
            case 2
                primaryLevelDB=-100;
                suppressorLevelDB=suppressorDB;
                plotGuide.subPlotNo=	2;     % initialize subplot count
                figureTitle=['suppressor alone:  ' num2str(suppressorFrequency) ' Hz;   ' num2str(suppressorLevelDB) ' dB SPL'];
            case 3
                primaryLevelDB=primaryDB;
                suppressorLevelDB=suppressorDB;
                plotGuide.subPlotNo=	3;     % initialize subplot count
                figureTitle=['probe + suppressor'];
        end
        
        % primary BF tone
        time1=dt: dt: duration;
        amp=10^(primaryLevelDB/20)*28e-6;
        inputSignal=amp*sin(2*pi*primaryToneFrequency*time1);
        rampDuration=.005; rampTime=dt:dt:rampDuration;
        ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(time1)-length(rampTime))];
        inputSignal=inputSignal.*ramp;
        inputSignal=inputSignal.*fliplr(ramp);
        
        % suppressor
        tone2Duration=duration/2; % s
        time2= dt: dt: tone2Duration;
        % B: tone on
        amp=10^(suppressorLevelDB/20)*28e-6;
        inputSignal2=amp*sin(2*pi*suppressorFrequency*time2);
        rampDuration=.005; rampTime=dt:dt:rampDuration;
        ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ones(1,length(time2)-length(rampTime))];
        inputSignal2=inputSignal2.*ramp;
        % A: initial silence (delay to suppressor)
        silence=zeros(1,length(time2));
        inputSignal2=[silence inputSignal2];
        
        % add tone and suppressor components
        inputSignal=inputSignal+inputSignal2;
        
        % specify model parameters
        method=MAPparamsDEMO(BFlist, sampleRate);
        % parameter change (must be global to take effect)
        global   AN_IHCsynapseParams
        AN_IHCsynapseParams.mode=	'probability';
%         method.useEfferent=0;
        
        method.plotGraphs=	1;	   % please plot
        
        figure(4), subplot(3,1,conditionNo)
        plot(time1, inputSignal, 'k')% *************
        
        [ANresponse, method, A]=MAPsequenceSeg(inputSignal, method, moduleSequence);
        response{conditionNo}=ANresponse;
% 	F(frameCount) = getframe(gcf);
    end
    
%     peakResponse=max(max([response{1} response{2} response{3}]));
%     figure(3), subplot(4,1,1)
%     mesh((response{1}), [0 peakResponse])
%     title(['primary: ' num2str(primaryDB) ' dB SPL'])
%     zlim([0 5000]), view([-28 36])
%     
%     figure(3), subplot(4,1,2)
%     mesh((response{2}), [0 peakResponse])
%     title(['suppressor: ' num2str(suppressorDB) ' dB SPL'])
%     zlim([0 5000]),  view([-28 36])
%     
%     figure(3), subplot(4,1,3)
%     mesh((response{3}), [0 peakResponse])
%     title([' primary + suppressor '])
%     zlim([0 5000]), view([-28 36])
%     
%     figure(3), subplot(4,1,4)
%     sum= response{3}-response{1}-response{2};
%     mesh(sum)
%     zlim([-2000 2000])
%     title(' difference matrix (combined - primary - suppressor)')
%     view([-28 20])
%     disp( [num2str([ suppressorDB primaryDB max(max(sum))])] )
%     
%     colormap(jet)
    pause(1)
end

figure(4)
xlabel('time (s)')
figure(3)
xlabel('time (s)')

%  figure(99), movie(F, 1, 1)
% UTIL_showAllMAPStructures