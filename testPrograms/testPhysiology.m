function testPhysiology(BF,paramsName, paramChanges)
% e.g.
% testPhysiology(1000,'Normal', [])

restorePath=path;
addpath (['..' filesep 'MAP'])

if nargin<3
    paramChanges=[];
end

testOME(paramsName,paramChanges)
relativeFrequencies=[0.25    .5   .75  1  1.25 1.5    2];
testBM (BF, paramsName,relativeFrequencies,'spikes', paramChanges)
testRP(BF,paramsName,paramChanges)
testSynapse(BF,paramsName,paramChanges)
testFM(BF,paramsName,'spikes', paramChanges)
testPhaseLocking(paramsName,paramChanges)
testAN(BF,BF, -10:10:80,paramsName,paramChanges);
figure(14)
figure(15)

path(restorePath)
