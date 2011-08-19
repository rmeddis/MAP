function timeStamp=UTIL_timeStamp

% returns a string showing time now and current date
% the string should be suitable for creating fileNames
% i.e. no ':' and no'-' characters

today=date;
idx=findstr('-',today);
today(idx)='_';

timeNow=clock;
timeNow=[num2str(timeNow(4)) ':' num2str(timeNow(5))];
idx=findstr(':',timeNow);
timeNow(idx)='_';

timeStamp=[timeNow 'hr' today ];
