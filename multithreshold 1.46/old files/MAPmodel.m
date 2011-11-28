function [modelResponse, MacGregorResponse]=MAPmodel( MAPplot, method)

global experiment stimulusParameters audio withinRuns
global outerMiddleEarParams DRNLParams AN_IHCsynapseParams

savePath=path;
addpath(['..' filesep 'MAP'], ['..' filesep 'utilities'])
modelResponse=[];
MacGregorResponse=[];

% mono only (column vector)
audio=audio(:,1)';

% if stop button pressed earlier
if experiment.stop, return, end

% -------------------------------------------------------------- run Model
MAPparamsName=experiment.name;
showPlotsAndDetails=experiment.MAPplot;
AN_spikesOrProbability='spikes';

% [response, method]=MAPsequenceSeg(audio, method, 1:8);
global ICoutput dtSpikes
    MAP1_14(audio, 1/method.dt, method.nonlinCF,...
        MAPparamsName, AN_spikesOrProbability);
    
if showPlotsAndDetails
    options.printModelParameters=0;
    options.showModelOutput=1;
    options.printFiringRates=1;
    options.showACF=0;
    options.showEfferent=0;
    UTIL_showMAP(options)
end

% No response,  probably caused by hitting 'stop' button
if isempty(ICoutput), return, end

% MacGregor response is the sum total of all final stage spiking
MacGregorResponse= sum(ICoutput,1);                 % use IC

% ---------------------------------------------------------- end model run

dt=dtSpikes;
time=dt:dt:dt*length(MacGregorResponse);

% group delay on unit response
MacGonsetDelay= 0.004;
MacGoffsetDelay= 0.022;

% now find the response of the MacGregor model during the target presentation + group delay
switch experiment.threshEstMethod
    case {'2I2AFC++', '2I2AFC+++'}
        idx= time>stimulusParameters.testTargetBegins+MacGonsetDelay ...
            & time<stimulusParameters.testTargetEnds+MacGoffsetDelay;
        nSpikesTrueWindow=sum(MacGregorResponse(:,idx));
        idx=find(time>stimulusParameters.testNonTargetBegins+MacGonsetDelay ...
            & time<stimulusParameters.testNonTargetEnds+MacGoffsetDelay);
        nSpikesFalseWindow=sum(MacGregorResponse(:,idx));
        % nSpikesDuringTarget is +ve when more spikes are found
        %   in the target window
        difference= nSpikesTrueWindow-nSpikesFalseWindow;

        if difference>0
            % hit
            nSpikesDuringTarget=experiment.MacGThreshold+1;
        elseif    difference<0
            % miss (wrong choice)
            nSpikesDuringTarget=experiment.MacGThreshold-1;
        else
            if rand>0.5
                % hit (random choice)
                nSpikesDuringTarget=experiment.MacGThreshold+1;
            else
                % miss (random choice)
                nSpikesDuringTarget=experiment.MacGThreshold-1;
            end
        end
        disp(['level target dummy decision: ' ...
            num2str([withinRuns.variableValue nSpikesTrueWindow ...
            nSpikesFalseWindow  nSpikesDuringTarget], '%4.0f') ] )

    otherwise
        % idx=find(time>stimulusParameters.testTargetBegins+MacGonsetDelay ...
        %         & time<stimulusParameters.testTargetEnds+MacGoffsetDelay);
        % no delay at onset
        idx=find(time>stimulusParameters.testTargetBegins +MacGonsetDelay...
            & time<stimulusParameters.testTargetEnds+MacGoffsetDelay);
        nSpikesDuringTarget=sum(MacGregorResponse(:,idx));
        
        % find(MacGregorResponse)*dt-stimulusParameters.stimulusDelay
        timeX=time(idx);
end

% now find the response of the MacGregor model at the end of the masker
idx2=find(time>stimulusParameters.testTargetBegins-0.02 ...
    & time<stimulusParameters.testTargetBegins);
if ~isempty(idx2)
    maskerRate=mean(mean(MacGregorResponse(idx2)));
else
    %e.g. no masker
    maskerRate=0;
end

if experiment.MAPplot
    % add vertical lines to indicate target region
    figure(99), subplot(6,1,6)
    hold on
    yL=get(gca,'YLim');
    plot([stimulusParameters.testTargetBegins + MacGonsetDelay ...
        stimulusParameters.testTargetBegins   + MacGonsetDelay],yL,'r')
    plot([stimulusParameters.testTargetEnds   + MacGoffsetDelay ...
        stimulusParameters.testTargetEnds     + MacGoffsetDelay],yL,'r')
end

% specify unambiguous response
switch experiment.paradigm
    case 'gapDetection'
        gapResponse=(maskerRate-nSpikesDuringTarget)/maskerRate;
        if gapResponse>0.2
            modelResponse=2;    % gap detected
        else
            modelResponse=1;    % gap not detected
        end
        [nSpikesDuringTarget maskerRate gapResponse modelResponse]
        figure(22), plot(timeX,earObject(idx))
    otherwise
        if nSpikesDuringTarget>experiment.MacGThreshold
            modelResponse=2;    % stimulus detected
        else
            modelResponse=1;    % nothing heard (default)
        end
end


path(savePath)
