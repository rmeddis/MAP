function UTIL_cascadePlot(toPlot, colValues)
% % useful code
[nChannels nLags]=size(toPlot);

% cunning code to represent channels as parallel lines
[nRows nCols]=size(toPlot);
if nChannels<2
    error('UTIL_cascadePlot: only one row found')
end

% a is the height to be added to each channel
a=max(max(toPlot))*(0:nRows-1)';

% peakGain emphasises the peak height
% peaks can be higher than the space between channels
peakGain=10;
x=peakGain*toPlot+repmat(a,1,nCols);
x=nRows*x/max(max(x));

for row=1:nRows-1
    x(row,:)=min(x(row,:), x(row+1,:));
end

plot(colValues,   x','k')
ylim([0 nRows])

