function [mappediih,centerfreqs] = map_iih_onto_log(iih,numchannels,sfreq)

%function to map an IPIH from the interval axis onto the frequency axis
%with no overlap and mean firing rate.
%
% Tim Juergens, September 2011
%
% input: iih: IPIH with dimensions interval (1) and time step (2)
%             the first dimension translates to time using the actual
%             sampling frequency
%        numchannels: number of desired channels
%        sfreq: sampling frequency
% output: mappediih: IPIH with dimensions frequency channel (1) and time
%                    step (2)
%         centerfreqs: center frequencies of the channels of mappediih

ctr_intervals = [1/sfreq:1/sfreq:size(iih,1)/sfreq];
lowestBF=1/ctr_intervals(end);
highestBF=10000;

borderfreqs = logspace(log10(lowestBF),log10(highestBF),numchannels+1);

for iCounter = 1:numchannels
    centerfreqs(iCounter)=(borderfreqs(iCounter)+borderfreqs(iCounter+1))/2;
end

for iCounter = 1:length(borderfreqs) %find the indices that correspond to the borderfrequencies of the BF filters
    [tmp,channelbordersIPIindex(iCounter)]=min(abs(borderfreqs(iCounter)-1./ctr_intervals));
end

for iCounter = 1:length(centerfreqs)
    %mapping with one interval sample overlap
    mappediih(iCounter,:) = mean(iih(channelbordersIPIindex(end-iCounter+1):channelbordersIPIindex(end-iCounter),:));
end



% OPTIONAL PLOTTING
            figure
            YTickIdx = 1:floor(numel(centerfreqs)/6):numel(centerfreqs);
            YTickIdxRev = numel(centerfreqs)+1-YTickIdx;
            if ~isempty(gca)
                axes(gca);  %#ok<MAXES>
                imagesc(mappediih);         
                set(gca, 'YTick', YTickIdx);                
                set(gca, 'YTickLabel', num2str(   centerfreqs(YTickIdxRev)', '%0.0f' ));
                ylabel('cf in Hz')                
            end
