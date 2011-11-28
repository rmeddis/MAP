function testPhysiologyProb(BF,paramsName, paramChanges)
% e.g.
% testPhysiologyProb(1000,'Normal', [])

restorePath=path;
addpath (['..' filesep 'MAP'])

if nargin<3
    paramChanges=[];
end

if nargin==0
    error('testPhysiologyProb must be called from the command line')
end


disp('testPhysiologyProb...........computing')

disp('testOME...........computing')
testOME(paramsName, paramChanges)

disp('testBM...........computing')
relativeFrequencies=[0.25    .5   .75  1  1.25 1.5    2];
testBM (BF, paramsName,relativeFrequencies,'probability', paramChanges)

disp('testRP...........computing')
testRP2(paramsName,paramChanges)

disp('testSynapse...........computing')
testSynapse(BF,paramsName, 'probability', paramChanges)

disp('testFM...........computing')
testFM(BF,paramsName,'probability', paramChanges)

disp('testANprob...........computing')
testANprob(BF,BF, -10:10:80,paramsName, paramChanges);

figure(4)

path(restorePath)
