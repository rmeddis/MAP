function outsignal = cutsignal(insignal,samplingfrequency,vocabularyset)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This function cuts the silence before and after the audiosignal
%
%  (c) Tim Jürgens, Medizinische Physik, Feb.2006
%
%  usage: outsignal = cutsignal(insignal,samplingfrequency,vocabularyset)
%
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial_energythreshold = 0.11;
% final_energythreshold = 0.05;
% 
% % choose initial and final delay due to different vocabularysets
% if (strcmp(vocabularyset,'f_f') > 0)
%     initial_delay = 0.17; % the speech started delay seconds before trespassing threshold (default 0.01, for f_f: 0.17)
%     final_delay = 0.2; %default 0.1, for f_f: 0.2
% else
%     initial_delay = 0.01; % the speech started delay seconds before trespassing threshold (default 0.01, for f_f: 0.2)
%     final_delay = 0.1; %default 0.1, for f_f: 0.2
% end
% 
% [initialsample_of_frame, energy] = compute_energy(insignal, samplingfrequency);
% for i = 1:length(initialsample_of_frame)
%     if (energy(i) > initial_energythreshold)
%         initialsample = initialsample_of_frame(i)-samplingfrequency*initial_delay;
%         break;
%     end
% end
% 
% time_inverted_energy = energy(end:-1:1); %turn signal around
% time_inverted_initialsample = initialsample_of_frame(end:-1:1);
% 
% for i = 1:length(time_inverted_initialsample)
%     if (time_inverted_energy(i) > final_energythreshold)
%         finalsample = time_inverted_initialsample(i)+samplingfrequency*final_delay;
%         break;
%     end
% end
% 
% if (finalsample > length(insignal))
%     finalsample = length(insignal);
% end
% outsignal = insignal(initialsample:finalsample);

if nargin < 3
    vocabularyset = 'a_a'; %default vocabularyset
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initial_threshold = 0.15; %amplitude threshold for detection of speech
final_threshold = 0.1;

% choose initial and final delay due to different vocabularysets
if (strcmp(vocabularyset,'f_f') > 0)
    initial_delay = 0.2; % the speech started delay seconds before trespassing threshold (default 0.01, for f_f: 0.2)
    final_delay = 0.2; %default 0.1, for f_f: 0.2
else
    initial_delay = 0.01; % the speech started delay seconds before trespassing threshold (default 0.01, for f_f: 0.2)
    final_delay = 0.1; %default 0.1, for f_f: 0.2
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%cut signal with taking initial delay into account
for(i = 1:length(insignal))
    if (abs(insignal(i)) > initial_threshold)
        initialsample = i - samplingfrequency*initial_delay;
        break;
    end
end


time_inverted_insignal = insignal(end:-1:1); %turn signal around
% cut it with taking a final delay into account
for(i = 1:length(time_inverted_insignal))
    if (abs(time_inverted_insignal(i)) > final_threshold)
        finalsample = i - samplingfrequency*final_delay;
        break;
    end
end

if (finalsample < 0)
    finalsample = 0;
end
if (initialsample < 1)
    initialsample = 1;
end

%% output %%%
outsignal = insignal(initialsample:end-finalsample);