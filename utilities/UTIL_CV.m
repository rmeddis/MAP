function [cv, cvTimes, allTimeStamps, allISIs]=  UTIL_CV(earObject, dt, timeStep)
% UTIL_CV computes coefficient of variation for multiple spike trains
% earObject must be logical 0/1.  Each row is a separate run
% CV is computed for successive time regions specified by timeStep

if ~islogical(earObject),error('UTIL_CV: input is not logical/ spikes'), end

[rows cols]=size(earObject);
totalDuration=cols*dt;

if nargin<3, timeStep=totalDuration/5; end

% identify all intervals
allISIs=[]; allTimeStamps=[];
for i=1:rows
    temp=find(earObject(i,:))*dt;    % find spikes
    isi{i}=diff(temp);                      % find ISIs
    timeStamps{i}=temp(2:end); % time of isi is time of second spike
    allISIs=[allISIs isi{i}];
    allTimeStamps=[allTimeStamps timeStamps{i}];
end

count=0;
cvTimes=0: timeStep:totalDuration-timeStep; % bin edges
for t= cvTimes
    % sort ISIs according to when they happened
    idx=find(allTimeStamps>t & allTimeStamps<=t+timeStep);
    count=count+1;
    if ~isempty(allISIs(idx))
        cv(count)=std(allISIs(idx))/mean(allISIs(idx));
    else
        cv(count)=0;
    end
end

