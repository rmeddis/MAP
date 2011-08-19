foreground = MTprofile17_14hr18_Aug_2011

fileName=['MTprofile' Util_timeStamp];
longTone=foreground.LongTone;
shortTone=foreground.ShortTone;
gaps=foreground.Gaps;
BFs=foreground.BFs;
TMC=foreground.TMC;
offBFs=foreground.MaskerRatio;
IFMCs=foreground.IFMCs;
profile2mFile(longTone, shortTone, gaps', BFs, TMC', offBFs', IFMCs',...
    fileName)
pause(1)
plotProfile(fileName, 'profile_CMA_L')
