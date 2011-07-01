function UTIL_testPhysiology(BF)
testOME('Normal')
relativeFrequencies=[0.25    .5   .75  1  1.25 1.5    2];
testBM (BF, 'Normal',relativeFrequencies)
testRP(BF,'Normal')
testSynapse(BF,'Normal')
testFM(BF,'Normal',1)
testPhaseLocking
testAN(BF,BF, 10:10:80);
figure(14)
figure(15)
