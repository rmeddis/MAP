function [PH, binTimes]=UTIL_periodHistogram(A, dt, frequency)
% UTIL_makePeriodHistogram converts a time sequence into a period histogram
%*********************
% The period is 1/frequency.
%
% usage:
%	PH=UTIL_periodHistogram(A, dt, frequency)
%
% Input:
% A is a channel x time matrix of spikes (or other stuff)
% frequency determines the period of the histogram
%
% Output
%  PH is a channel by periodhistogram matrix
% bintimes is useful for plotting the output
%

periodInSeconds=1/frequency;
[numChannels signalNpoints]=size(A);

% retrict data array to a multiple of the period.
pointsPerPeriod= round(periodInSeconds/dt);
NcompletePeriods=floor(signalNpoints/pointsPerPeriod);
totalPointsUsed=NcompletePeriods*pointsPerPeriod;

% check that the period is a whole number of epochs
aliasing=NcompletePeriods*(periodInSeconds/dt-pointsPerPeriod);

if aliasing>.1
    error('UTIL_periodHistogram: irregular period length')
end

if NcompletePeriods<1
    error('UTIL_periodHistogram: too few datapoints')
end

% transpose data so that time is down a column
A=A(:,1:totalPointsUsed)';

% knock it into shape
A=reshape(A,pointsPerPeriod, NcompletePeriods, numChannels);
% each period is a separate column
% imagesc(squeeze(A)) % should have horizontal stipe

% channels are now the third dimension.
PH=squeeze(sum(A,2))';
% PH=PH/NcompletePeriods;

binTimes=dt:dt:pointsPerPeriod*dt;
