% ------------------------------------------------------- psychometricFunction
function [psy, levelsBinVector, binFrequencies, nNo, nYes]=...
    psychometricFunction(levels, responses, psyBinWidth)
% computes a psychometric function from rwo vectors (levels and responses)
% responses is a binary vectory (1=yes)
% psyBinWidth is the bin width (dB) of the output psychometric function
% this function is called by both fitPsychometricFunction and subjGUI
% for this reason it should notbe bundled with fitPsychometricFunction

binFrequencies=[]; nNo=[]; nYes=[];
if min(levels)+abs(psyBinWidth) < max(levels)
    % create a set of bins for the psychometric function
    levelsBinVector= min(levels): abs(psyBinWidth): max(levels);
else
    psy=[]; levelsBinVector=[];
    return
end

idx0=find(responses==0);	% isolate 'no'
z=levels(idx0);
nNo=hist(z, levelsBinVector);

idx1=find(responses>0);	% isolate 'yes'
y=levels(idx1);
nYes=hist(y, levelsBinVector);

if sum(nNo)==0 | sum(nYes)==0
        psy=[]; levelsBinVector=[];
    return             % yesses and nos required for a function
end

binFrequencies=nNo+nYes;

warning off MATLAB:divideByZero
psy=nYes./binFrequencies;        % psy is the proportion of 'yes' responses
warning on MATLAB:divideByZero
lastwarn('');

idx=~isnan(psy);	%remove empty bins
idx=(nYes>0) |(nNo>0);	%remove empty bins
psy=psy(idx);
levelsBinVector=levelsBinVector(idx);
binFrequencies=binFrequencies(idx);
nNo=nNo(idx);
nYes=nYes(idx);

% [nNo' nYes']
% [levelsBinVector' psy']
% plot(levelsBinVector,psy,['o'])