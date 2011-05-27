function VS=UTIL_vectorStrength(PH)
% UTIL_vectorStrength converts a period histogram to a vector strength measure.
% **************
% using the Johnson(1980) formula.
% usage:
% 	VS=UTIL_vectorStrength(PH);
%
% Input arguments:
% 	PH is channels x  period histogram matrix
%
% Output argumkents:
% VS is a column of vector strengths, one per channel
[numChannels K]=size(PH);

% a=cos(2*pi*(1/K:1/K:1));
% b=sin(2*pi*(1/K:1/K:1));
% VS=(sum(PH.*repmat(a,numChannels,1),2).^2 + sum(PH.*repmat(b,numChannels,1),2).^2).^0.5;
% warning off MATLAB:divideByZero
% VS=VS./sum(PH,2);
% warning on MATLAB:divideByZeroa=cos(2*pi*(1/K:1/K:1));


a=cos(2*pi*(1/K:1/K:1));
b=sin(2*pi*(1/K:1/K:1));
N=sum(PH,2);
warning off MATLAB:divideByZero
VS=((sum(PH.*repmat(a,numChannels,1),2)./N).^2 ...
    + (sum(PH.*repmat(b,numChannels,1),2)./N).^2).^0.5;
warning on MATLAB:divideByZero
