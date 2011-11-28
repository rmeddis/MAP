function [pooledIPIdata,ctr] = poolIPI_across_channels(IPIhisttime,IPIhistweight)

%function that pools IPIdata across all channels in order to plot the data
%as in Secker-Walker JASA 1990 Fig. 6

%plots IPI-Histograms in a matrix display as in Secker-Walker JASA 1990
%Fig. 4
%
%   Tim Jürgens, February 2011
%
% input:  IPIhisttime:     matrix containing the interval times found in the
%                          analysis
%                          first dimension: frequency channel
%                          second dimension: time step (hop position)
%                          third dimension: interval dimension (max 3)
%         IPIhistweight:   matrix containing the interval weights found in the
%                          analysis, same dimensions as IPIhisttime
% output: pooledIPIdata:   matrix containing the pooled and weighted IPI
%                          histograms as a function of time (hop position)
%                          first dimension: time step (hop position)                           
%                          second dimension: IPI intervals
%         ctr:             center values of classes of IPI intervals
% ATTENTION: depending on the amplitude of the signals the value
% 'verticalshift' might have to be adjusted in order to see structures in
% the plot

verticalshift =   3000; %0.008;%vertical shift of single time series in z-coordinate units
ctr = [0:0.02:5].*1e-3; %classes (center-values) of 20 microsec width

%preallocation of variables
pooledIPIdata = zeros(size(IPIhisttime,2),length(ctr));
smoothed_pooledIPI = zeros(size(IPIhisttime,2),length(ctr)-5);

for iCounter = 1:size(IPIhisttime,2) %one for time spacing
    %cannot use matlabs hist function because weighting must be applied
   
    tmpIPIhisttime = squeeze(IPIhisttime(:,iCounter,:));
    tmpIPIhistweight = squeeze(IPIhistweight(:,iCounter,:));
    tmpIPIhisttime = tmpIPIhisttime(:);
    tmpIPIhistweight = tmpIPIhistweight(:);
    for jCounter = 1:length(tmpIPIhisttime)
        %look which class
        [tmp1,classindex] = min(abs(tmpIPIhisttime(jCounter)-ctr));
        pooledIPIdata(iCounter,classindex) = pooledIPIdata(iCounter,classindex)+tmpIPIhistweight(jCounter);
    end
    
end

%smooth data using a 5-point hamming window
hamm_window = hamming(5);
for iCounter = 1:size(pooledIPIdata,1)
    for jCounter = 3:length(ctr)-2 %start with 3 and end with 2 samples
        %less the length of ctr in order not to get in conflict with the length of
        %the hamm_window
        smoothed_pooledIPI(iCounter,jCounter-2) = ...
            pooledIPIdata(iCounter,(jCounter-2):(jCounter+2))*hamm_window./sum(hamm_window);
    end
end
smoothed_ctr = ctr(3:end-2);


figure;

Tickvector = [];
TickLabels = [];
for iCounter = 1:size(IPIhisttime,2)
    hold on;
    verticalposition = verticalshift*(iCounter-1);
    %multiply by 1000 to set abscissa to ms units
    plot(1000.*smoothed_ctr, ...
        smoothed_pooledIPI(iCounter,:)+verticalposition, ...
        'k','LineWidth',2);  
    if mod(iCounter,5) == 1 %set best frequency as a label every 10 channels
        Tickvector = [Tickvector verticalposition];
        TickLabels = [TickLabels; (iCounter-1)*3]; %time spacing is every 3ms
    end
end
set(gca,'yTick',Tickvector,'yTickLabel',num2str(TickLabels,'%4.0f'));
ylabel('Stimulus Time (ms)');
xlabel('Interval (ms)');
xlim([min(1000*ctr) max(1000*ctr)]);
ylim([-verticalshift verticalposition+3*verticalshift]);
box on;
