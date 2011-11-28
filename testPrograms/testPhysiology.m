function testPhysiology(BF,paramsName, paramChanges)
% e.g.
% testPhysiology(1000,'Normal', {})

restorePath=path;
addpath (['..' filesep 'MAP'])

if nargin<3,  paramChanges=[]; end

if nargin<2, paramsName='Normal'; end

if nargin<1, BF=1000; end

disp('testPhysiology...........computing')
disp('testOME...........computing')
testOME(paramsName,paramChanges)

relativeFrequencies=[0.25    .5   .75  1  1.25 1.5    2];
disp('testBM...........computing')
testBM (BF, paramsName,relativeFrequencies,'spikes', paramChanges)

disp('testRP...........computing')
testRP2(paramsName,paramChanges)

disp('testSynapse...........computing')
testSynapse(BF,paramsName, 'spikes', paramChanges)

disp('testFM...........computing')
testFM(BF,paramsName,'spikes', paramChanges)

disp('testPhaseLocking...........computing')
testPhaseLocking(paramsName,paramChanges)

disp('testAN...........computing')
testAN(BF,BF, -10:10:80,paramsName,paramChanges);

% put these figures on top
figure(14)
figure(15)

path(restorePath)
