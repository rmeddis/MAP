function Exp_Tutorial_2(isMasterNode)

% This tutorial recycles a HMM

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the basic experiment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expName = 'Tutorial';
dataFolderPrefix = 'recycle_featR';
if isunix
    expFolderPrefix = '/scratch/nrclark/exps/';
else
    expFolderPrefix = 'D:\Exps';
end

% expFolderPrefix = pwd;
expFolder = fullfile(expFolderPrefix,expName);
hmmFolder = fullfile(expFolder,'hmm');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort out the training (LEARNING) condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
learnFolder = fullfile(expFolder,'featL');

xL = cJob('L', learnFolder);

xL.participant = 'Normal';
xL.MAPparamChanges= {'DRNLParams.rateToAttenuationFactorProb=0;', 'OMEParams.rateToAttenuationFactorProb=0;' };

xL.noiseLevToUse   =  -200;
xL.speechLevToUse  =  60;

xL.MAPopHSR = 1;
xL.MAPopMSR = 0;
xL.MAPopLSR = 0;


xL.numCoeff = 14;
xL.removeEnergyStatic = 0;

%%%%% Group of params that will influence simulation run time %%%%%%%
xL.numWavs = 12; %MAX=8440
testWavs = 6; %MAX = 358
nzLevel = [-200 40:10:70];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xL.noisePreDur = 1;
xL.noisePostDur = 0.1;
xL.truncateDur  = xL.noisePreDur-0.1; 
xL.noiseName = 'pink_demo';


% if isMasterNode && ~isdir(xL.opFolder)
%     mkdir(xL.opFolder);
%     xL = xL.assignFiles;
%     xL.storeSelf;
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort out the testing (RECOGNITION) conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recConditions = numel(nzLevel);

tmpIdx=0;
for nn = 0*recConditions+1:1*recConditions    
    tmpIdx=tmpIdx+1;
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,[dataFolderPrefix num2str(nn)]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = testWavs; %MAX = 358
    xR{nn}.noiseLevToUse = nzLevel(tmpIdx); 
    xR{nn}.MAPparamChanges= {'DRNLParams.a=400;'};

    
    %Now just to wrap it up ready for processing
    if isMasterNode && ~isdir(xR{nn}.opFolder)
        mkdir(xR{nn}.opFolder);
        xR{nn} = xR{nn}.assignWavPaths('R');
        xR{nn} = xR{nn}.assignFiles;
        xR{nn}.storeSelf;
    end
end

tmpIdx=0;
for nn = 1*recConditions+1:2*recConditions    
    tmpIdx=tmpIdx+1;
    xR{nn} = xL; %simply copy the "Learn" object and change it a bit below
    recFolder = fullfile(expFolder,[dataFolderPrefix num2str(nn)]);
    xR{nn}.opFolder = recFolder;    
    
    %These are the interesting differences between training and testing
    xR{nn}.numWavs = testWavs; %MAX = 358
    xR{nn}.noiseLevToUse = nzLevel(tmpIdx); 
    xR{nn}.MAPparamChanges= {'DRNLParams.a=400;'};
    
    xR{nn}.mainGain = [27.2013;   26.0797;   26.0939;   26.7997;   26.0520];  
    xR{nn}.TCdBO    = [37;   37;   37;   37;   37];      %Compression thresholds (in dB OUTPUT from 2nd filt)
    xR{nn}.TMdBO    = [20;   20;   20;   20;   20];      %MOC thresholds (in dB OUTPUT from 2nd filt)
    xR{nn}.ARthresholddB = 85;       % dB SPL (input signal level) =>200 to disable
    xR{nn}.MOCtau = 1;
    xR{nn}.useAid = 1;
    
    %Now just to wrap it up ready for processing
    if isMasterNode && ~isdir(xR{nn}.opFolder)
        mkdir(xR{nn}.opFolder);
        xR{nn} = xR{nn}.assignWavPaths('R');
        xR{nn} = xR{nn}.assignFiles;
        xR{nn}.storeSelf;
    end
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Generate features **
% This is the time consuming, processing intensive portion of the program.
% Nodes that are not the master node are only interested in the opFolder
% member of the jobjects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% worker(xL.opFolder);
maxConds = nn;
if ~isMasterNode %dont bother wasting master node effort on generating testing features (for now)
    for nn = 1:maxConds
        worker(xR{nn}.opFolder);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Train and test the recogniser - a job for the master node only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isMasterNode    
%     while(~all(xL.todoStatus==2))        
%         disp('Waiting on straggler nodes to complete their jobs before HMM is trained . . .')
%         pause(30); %Wait for 30 seconds before looking again
%         xL.lockJobList;
%         xL = xL.loadSelf; %Reload incase changed
%         xL.unlockJobList;
%     end
    y = cHMM(hmmFolder);    
    y.numCoeff = 14*3;
%     y.createSCP(xL.opFolder)
%     y.createMLF(xL.opFolder)
%     y.train(xL.opFolder) %This node can be busy training, even if other jobs are being processed for testing
    
    % ALLOW MASTER NODE TO MUCK IN WITH GENERATING TESTING FEATURES ONCE
    % HMM HAS BEEN TRAINED
    for nn = 1:maxConds
        worker(xR{nn}.opFolder);
    end    
    
    xR{end}.lockJobList;
    xR{end} = xR{end}.loadSelf; %Reload changes
    xR{end}.unlockJobList;
    while(~all(xR{end}.todoStatus==2))        
        disp('Waiting on straggler nodes to complete their jobs before HMM is tested . . .')
        pause(30); %Wait for 30 seconds before looking again
        xR{end}.lockJobList;
        xR{end} = xR{end}.loadSelf; %Reload incase changed
        xR{end}.unlockJobList;
    end
      
    for nn = 1:maxConds
        y.createSCP(xR{nn}.opFolder);
        y.test(xR{nn}.opFolder);
    end
    
    %Show all of the scores in the command window at the end
    for nn = 1:maxConds
        y.score(xR{nn}.opFolder);
    end
end
