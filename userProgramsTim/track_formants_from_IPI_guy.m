function [iih,IPIhisttime,IPIhistweight]=track_formants_from_IPI_guy(IFRAN_pattern, sfreq)
% 
% tracks the formants according to an analysis proposed in Secker-Walker
% JASA 1990, section V.A
% Tim Jürgens, February 2011, code from Guy Brown included
% 
% input:  IFRAN_pattern: pattern of the auditory model (dependend on the number of modules used)
%                        first dimension: frequency channel, 
%                        second dimension: time (samples)
%         sfreq:         sampling frequency
% output: iih:           interpeak-interval histogram, matrix very similar
%                        the plot 5 in the Secker-Walker paper
%
%
%


time_axis = 0:1/sfreq:(size(IFRAN_pattern,2)-1)/sfreq;

%find how many samples of AN_pattern are 10ms and 3ms
%one_sample_is_a_time_of = time_axis(2);
[tmp, start_time_index] = min(abs(0-time_axis));
[tmp, stop10_time_index] = min(abs(0.01-time_axis));
number_of_samples10ms = stop10_time_index - start_time_index;

[tmp, stop3_time_index] = min(abs(0.003-time_axis));
number_of_samples3ms = stop3_time_index - start_time_index;
every_3ms = 1:number_of_samples3ms:size(IFRAN_pattern,2)-number_of_samples10ms;

hamm_window = hamming(11);
halfHamming = (length(hamm_window)-1)/2;

% window normalization

norm = conv(ones(1,floor(number_of_samples10ms/2)),hamm_window);
norm = norm(5+1:end-5)';
win_size = number_of_samples10ms;
half_win_size = floor(win_size/2);
hop_size = number_of_samples3ms;


%pre-allocation due to speed
%Acorr = zeros(size(IFRAN_pattern,1),size(every_3ms,2),number_of_samples10ms*2+1);
%RAcorr = zeros(size(IFRAN_pattern,1),size(every_3ms,2),number_of_samples10ms*2+1);
%SRAcorr = zeros(size(IFRAN_pattern,1),size(every_3ms,2),number_of_samples10ms*2+1-10);
IPIhisttime = zeros(size(IFRAN_pattern,1),size(every_3ms,2),3);
IPIhistweight = zeros(size(IFRAN_pattern,1),size(every_3ms,2),3);  %maximum 3 peaks from the SRA
iih = zeros(half_win_size,size(every_3ms,2)+1);




for iCounter = 1:size(IFRAN_pattern,1) %each channel
    fprintf('Channel No. %i\n',iCounter);
    %time_counter = 1;
    %for jCounter = every_3ms %every 3ms time segment
    
    
    
    %% Guy's code
    % enframe this signal
    
    frames = enframe(IFRAN_pattern(iCounter,:),win_size,hop_size);
    
    % compute the autocorrelation
    
    acf = real(ifft(abs(fft(frames,[],2)).^2,[],2));
    acf(acf<0)=0;
    acf = sqrt(acf(:,1:half_win_size));
    
    % smooth with hamming window and take the root
    
    for frame=1:size(acf,1)
        
        %%debug
        %if iCounter == 130
        %    disp('here');
        %end
        
        
        sra = conv(acf(frame,:),hamm_window);
        sra = sra(halfHamming+1:end-halfHamming)./norm';
        df = [0 ; diff(sra')];
        idx = find((df(1:end-1)>=0)&(df(2:end)<0));
        % interpolate
        a=df(idx);
        b=df(idx+1);
        idx = (idx-1+a./(a-b));
        % get rid of a zero peak, if it exists
        idx = idx(idx>1);
        % peak values corresponding to these intervals
        amp = interp1(1:length(sra),sra,idx,'linear');
        % if required, remove peaks that lie below the mean sra
        % note that we disregard the value at zero delay
        %if (params.removePeaksBelowMean)
        valid = find(amp>mean(sra(2:end)));
        idx = idx(valid);
        amp = amp(valid);
        %end
        % only use the first four peaks (three intervals)
        idx = idx(1:min(4,length(idx)));
        % find the intervals
        interval = diff(idx);
        % now histogram the intervals
        if (~isempty(interval))
            for k=1:length(interval),
                iih(round(interval(k)),frame) = iih(round(interval(k)),frame)+amp(k);
                IPIhisttime(iCounter,frame,k) = interval(k)/sfreq;
                IPIhistweight(iCounter,frame,k) = amp(k);
            end
        end
        
    end
    
    
    
    
    %% end Guy's code
    
    
    %         %take the autocorrelation (ACF) of a 10ms-segment of each channel
    %         Acorr(iCounter,time_counter,:) = xcorr(IFRAN_pattern(iCounter,jCounter:number_of_samples10ms+jCounter),'biased'); %biased scales the ACF by the reciprocal of the length of the segment
    %         %root calculation
    %         RAcorr(iCounter,time_counter,:) = sqrt(abs(Acorr(iCounter,time_counter,:)));
    %
    %         %smoothing using the 11-point hamming window
    %         for kCounter = 6:size(RAcorr(iCounter,time_counter,:),3)-5 %start with 6 and end with 5 samples
    %             %less the length of time_axis not to get in conflict with the length of
    %             %the hamm_window
    %             SRAcorr(iCounter,time_counter,kCounter-5) = ...
    %                 squeeze(RAcorr(iCounter,time_counter,(kCounter-5):(kCounter+5)))'*hamm_window./sum(hamm_window);
    %         end
    %
    %         %mean value of actual SRA
    %         SRA_mean = mean(SRAcorr(iCounter,time_counter,:));
    %
    %         %find signed zero-crossings of the first derivative (=difference)
    %         z_crossings_indices = find(diff(sign(diff(squeeze(SRAcorr(iCounter,time_counter,:))))) < 0)+1; %+1 is necessary, because diff shortens vector by 1
    %         middle_index = ceil(size(SRAcorr(iCounter,time_counter,:),3)/2);
    %
    %         validCounter = 1;
    %         valid_z_crossings_indices = [];
    %         %find valid zero-crossings (peak higher than meanvalue and within first 5 ms of SRA)
    %         for lCounter = 1:length(z_crossings_indices)
    %             if (SRAcorr(iCounter,time_counter,z_crossings_indices(lCounter)) > SRA_mean) && ...
    %                     (abs(z_crossings_indices(lCounter)-middle_index) < round(number_of_samples10ms/2));
    %                 valid_z_crossings_indices(validCounter) = z_crossings_indices(lCounter);
    %                 validCounter = validCounter+1;
    %             end
    %         end
    %
    %         %find main peak in the ACF
    %         [tmp,index_of_z_crossings_main_index] = min(abs(middle_index-valid_z_crossings_indices));
    %         if ~tmp == 0
    %             disp('middle peak not appropriately found');
    %         end
    %
    %         %%% for debugging
    % %         if iCounter == 130
    % %             disp('here');
    % %             figure, plot(squeeze(SRAcorr(iCounter,time_counter,:)));
    % %             hold on, plot([1 length(squeeze(SRAcorr(iCounter,time_counter,:)))],[SRA_mean SRA_mean],'r-');
    % %         end
    %         %%%
    %
    %         %generate IPI-histogram: take the first 3 intervals of SRAcorr
    %         %(positive delay) in the first 5 ms
    %         histcounter = 1;
    %         for lCounter = index_of_z_crossings_main_index+1:min([length(valid_z_crossings_indices(index_of_z_crossings_main_index+1:end)) 3])+index_of_z_crossings_main_index
    %             sampledifference = abs(valid_z_crossings_indices(lCounter)-valid_z_crossings_indices(lCounter-1));
    %             %the difference between two adjacent peaks in the SRA is taken
    %             %as IPI estimate
    %             IPIhisttime(iCounter,time_counter,histcounter) = sampledifference*one_sample_is_a_time_of;
    %             %the amplitude of the SRA at the start of the SRA interval is
    %             %taken as the IPIweight
    %             IPIhistweight(iCounter,time_counter,histcounter) = SRAcorr(iCounter,time_counter,valid_z_crossings_indices(lCounter-1));
    %             histcounter = histcounter + 1;
    %         end
    
    %time_counter = time_counter+1;
end

