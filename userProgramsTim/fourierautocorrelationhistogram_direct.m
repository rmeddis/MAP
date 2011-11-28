function fach=fourierautocorrelationhistogram_direct(ANpattern,sfreq,plothandle)


time_axis = 0:1/sfreq:(size(ANpattern,2)-1)/sfreq;

%find how many samples of AN_pattern are 10ms and 3ms
%one_sample_is_a_time_of = time_axis(2);
[tmp, start_time_index] = min(abs(0-time_axis));
[tmp, stop20_time_index] = min(abs(0.020-time_axis));
number_of_samples20ms = stop20_time_index - start_time_index;

[tmp, stop3_time_index] = min(abs(0.003-time_axis));
number_of_samples3ms = stop3_time_index - start_time_index;
every_3ms = 1:number_of_samples3ms:size(ANpattern,2)-number_of_samples20ms;

hamm_window = hamming(11);
halfHamming = (length(hamm_window)-1)/2;

% window normalization

norm = conv(ones(1,floor(number_of_samples20ms/2)),hamm_window);
norm = norm(5+1:end-5)';
win_size = number_of_samples20ms;
half_win_size = floor(win_size/2);
hop_size = number_of_samples3ms;

%preallocation due to speed
fach = zeros(half_win_size,size(every_3ms,2)+1);

for iCounter = 1:size(ANpattern,1) %each channel
    fprintf('Channel No. %i\n',iCounter);
    %time_counter = 1;
    %for jCounter = every_3ms %every 3ms time segment
    
    
    
    %% Guy's code
    % enframe this signal
    
    frames = enframe(ANpattern(iCounter,:),win_size,hop_size);
    
    % compute the autocorrelation
    
    %acf = real(ifft(abs(fft(frames,[],2)).^2,[],2));
    %acf(acf<0)=0;
    %acf = sqrt(acf(:,1:half_win_size));
    
    % smooth with hamming window and take the root
    
    for frame=1:size(frames,1)
        
               
        smoothed_frame = conv(frames(frame,:),hamm_window);
        smoothed_frame = smoothed_frame(halfHamming+1:end-halfHamming);
        fsra = 20*log10(abs(fft(smoothed_frame-mean(smoothed_frame))));
        fsra = fsra(1:floor(length(fsra)/2));
        
        t = [0:1/sfreq:length(smoothed_frame)/sfreq-1/sfreq];
        frequency = [0:1/t(end):1/(2*(t(2)-t(1)))];
        %identify peaks in the fft
        df = [0 ; diff(fsra')];
        idx = find((df(1:end-1)>=0)&(df(2:end)<0));
%         % interpolate
%         a=df(idx);
%         b=df(idx+1);
%         idx = (idx-1+a./(a-b));
        [sorted,sortedindex]=sort(fsra(idx),'descend');
        % just take the three highest values of the fourier-transform
         valid_peak_index = sortedindex(1:min([length(sortedindex) 1]));
         amp = sorted(1:min([length(sortedindex) 1]));
         
         %store valid peaks according to amplitude in a histogram
         if (~isempty(valid_peak_index))
            for k=1:length(valid_peak_index),
                fach(idx(valid_peak_index(k)),frame) = fach(idx(valid_peak_index(k)),frame)+amp(k);            
            end
        end
         %transform index into frequencies
         
    end
end


%plot the result
if ~exist('plothandle'), plothandle=figure; end
maxfrequency = 4000;
[tmp,number_of_channels_to_display] = min(abs(frequency-maxfrequency));
frequency = frequency(1:number_of_channels_to_display);

set(gcf,'Currentaxes',plothandle);

YTickIdx = 1:floor(numel(frequency)/6):numel(frequency);
XTickIdx = 1:floor(numel(every_3ms)/6):numel(every_3ms);
YTickIdxRev = numel(frequency)+1-YTickIdx;
if ~isempty(gca)
    axes(gca);  %#ok<MAXES>
    imagesc(fach(1:number_of_channels_to_display,:));
    axis xy
    set(gca, 'YTick', YTickIdx);
    set(gca, 'YTickLabel', num2str(   frequency(YTickIdx)', '%0.0f' ));
    ylabel('frequency (Hz)')
    set(gca, 'XTick', XTickIdx);
    set(gca, 'XTickLabel', XTickIdx.*3);
    xlabel('Time (ms)')
end




