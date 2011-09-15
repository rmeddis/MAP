function [PSTH dt]=UTIL_PSTHmaker(inputData, dt, PSTHbinWidth)
% UTIL_PSTHmaker accumulates values into bins.
% No corrections are applied
% usage:
%	PSTH=UTIL_PSTHmaker(inputData, method)
%
% arguments
%	inputData is a channel x time matrix
%   dt is the sampling interval of the input signal (not clear why this is
%   returned?!)
%	PSTH is the reduced matrix, the sum of all elements in the bin
%
% e.g.
% [PSTH dt]=UTIL_PSTHmaker(inputData, dt, PSTHbinWidth)

[numChannels numDataPoints]= size(inputData);

% Multiple fibers is the same as repeat trials
% Consolidate data into a histogram 
dataPointsPerBin=round(PSTHbinWidth/dt);
if dataPointsPerBin<=1;
% 	Too few data points
	PSTH=inputData;
	return
end

numBins=floor(numDataPoints/dataPointsPerBin);
PSTH=zeros(numChannels,numBins);

% take care that signal length is an integer multiple of bin size
%  by dropping the last unuseable values
useableDataLength=numBins*dataPointsPerBin;
inputData=inputData(:,1:useableDataLength);

for ch=1:numChannels
	% Convert each channel into a matrix where each column represents 
	% the content of a single PSTH bin
	PSTH2D=reshape (inputData(ch,:), dataPointsPerBin, numBins );
	% and sum within each bin (across the rows
	PSTH(ch,:)=sum (PSTH2D,1);
end

