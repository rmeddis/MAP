%skript to test several stimuli with MAP1_14
% Tim Juergens, September 2011

%addpath('C:\Documents and Settings\tjurgens\My Documents\Dropbox\matlab');
addpath('C:\MAP\userPrograms');
parameterfile='Normal';

%Read the wavfiles
[da,sfreq] = wavread('..\wavFileStore\new-da-44khz.wav'); %artificial stimulus 
ba = wavread('..\wavFileStore\new-ba-44khz.wav'); %artificial stimulus

%set level of speech stimuli (RMS 1 equals 94 dB)
da_69dB = da./sqrt(mean(da.^2)).*10^(-(94-69)/20);
da_49dB = da./sqrt(mean(da.^2)).*10^(-(94-49)/20);
ba_69dB = ba./sqrt(mean(ba.^2)).*10^(-(94-69)/20);
%20*log10(sqrt(mean(ba_69dB.^2))/20e-6) %reference pressure: 20uPa


cd ..
cd MAP

%MAP and store AN output of HSR fibers in variables
MAP1_14(da_69dB,sfreq,-1,parameterfile,'probability');
global ANprobRateOutput
AN_da_69dB = ANprobRateOutput(1:size(ANprobRateOutput,1)/2,:);

 MAP1_14(da_49dB,sfreq,-1,parameterfile,'probability');
 global ANprobRateOutput
 AN_da_49dB = ANprobRateOutput(1:size(ANprobRateOutput,1)/2,:);
% 
 MAP1_14(ba_69dB,sfreq,-1,parameterfile,'probability');
 global ANprobRateOutput
 global savedBFlist
 AN_ba_69dB = ANprobRateOutput(1:size(ANprobRateOutput,1)/2,:);

%Do the IPIH analysis
 [iih_da_69dB,IPIhisttime_da_69dB,IPIhistweight_da_69dB] = track_formants_from_IPI_guy(AN_da_69dB,sfreq);
 %[iih_da_69dB,IPIhisttime_da_69dB,IPIhistweight_da_69dB] = IPIHextract(AN_da_69dB,sfreq);
 poolIPI_across_channels(IPIhisttime_da_69dB,IPIhistweight_da_69dB);
 map_iih_onto_log(iih_da_69dB,30,sfreq);
  [iih_da_49dB,IPIhisttime_da_49dB,IPIhistweight_da_49dB] = track_formants_from_IPI_guy(AN_da_49dB,sfreq);
%  poolIPI_across_channels(IPIhisttime_da_49dB,IPIhistweight_da_49dB);
%  
  [iih_ba_69dB,IPIhisttime_ba_69dB,IPIhistweight_ba_69dB] = track_formants_from_IPI_guy(AN_ba_69dB,sfreq);
%  [tmp,ctr]=poolIPI_across_channels(IPIhisttime_ba_69dB,IPIhistweight_ba_69dB);
%  
 
 
 
 %%%%%% OLLO stimuli %%%%%%%%%%%%%%%%%
%  OLLOwavfiles = {'dahd','bahb','dehd','behb','died','bieb','dohd','bohb','duhd','buhb','atta','ascha','assa'};
%  OLLO_level = 69; %level of OLLO files in dB SPL
%  paramChanges = [];%{'DRNLParams.rateToAttenuationFactorProb = 0;'};
%  
% for iCounter = 1:length(OLLOwavfiles)
%     
%     %read the stimuli
%     eval(['[' OLLOwavfiles{iCounter} ',sfreq_OLLO] = wavread([''..\wavFileStore\S02M_L' ...
%         sprintf('%3.3i',sloga2iloga(OLLOwavfiles{iCounter})) '_V6_M1_N2_CS0.wav'']);'])
%     
%     %delete preceding and subsequent silence and resample to 44100 Hz
%     %sampling frequency
%     eval([OLLOwavfiles{iCounter} ' = cutsignal(' OLLOwavfiles{iCounter} ',sfreq_OLLO,''d_d'');']);
%     eval([OLLOwavfiles{iCounter} ' = resample(' OLLOwavfiles{iCounter} ',sfreq,sfreq_OLLO);']);
%     
%     %set level
%     partfilename = [OLLOwavfiles{iCounter} num2str(OLLO_level) 'dB'];
%     eval([partfilename ' = ' OLLOwavfiles{iCounter} ...
%         './sqrt(mean(' OLLOwavfiles{iCounter} '.^2)).*10^(-(94-' num2str(OLLO_level) ')/20);']);
%     
%     %use MAP
%     eval(['MAP1_14(' partfilename ',sfreq,-1,parameterfile,''probability'', paramChanges );']);
%     global ANprobRateOutput
%     eval(['AN_' partfilename ' = ANprobRateOutput(1:size(ANprobRateOutput,1)/2,:);']);
%     
%     %Do the IPIH analysis
%     eval(['[iih_' partfilename ',IPIhisttime_' partfilename ',IPIhistweight_' ...
%         partfilename '] = track_formants_from_IPI_guy(AN_' partfilename ',sfreq);']);
%     eval(['poolIPI_across_channels(IPIhisttime_' partfilename ',IPIhistweight_' partfilename ');']);
%     title(partfilename);
%     %set(gca,'Title',partfilename);
%     xlabel('Interval (ms)')
%     ylabel('Stimulus time (ms)');
%     
%     eval(['map_iih_onto_log(iih_' partfilename ',30,sfreq);']);
%     title(partfilename);
% end


%%%%%%%%%%%%% OLLO stimuli from different speakers
%  OLLOwavfiles = {'S01F_L111_V6_M1_N2_CS0.wav','S02M_L111_V6_M1_N2_CS0.wav'};
%  OLLO_level = 69; %level of OLLO files in dB SPL
% for iCounter = 1:length(OLLOwavfiles)
%     
%     partfilename = [OLLOwavfiles{iCounter}(1:end-4) num2str(OLLO_level) 'dB'];
%     %read the stimuli
%     eval(['[' partfilename ',sfreq_OLLO] = wavread([''..\wavFileStore\' OLLOwavfiles{iCounter} ''']);'])
%     
%     %delete preceding and subsequent silence and resample to 44100 Hz
%     %sampling frequency
%     eval([partfilename ' = cutsignal(' partfilename ',sfreq_OLLO,''d_d'');']);
%     eval([partfilename ' = resample(' partfilename ',sfreq,sfreq_OLLO);']);
%     
%     %set level
%     
%     eval([partfilename ' = ' partfilename ...
%         './sqrt(mean(' partfilename '.^2)).*10^(-(94-' num2str(OLLO_level) ')/20);']);
%     
%     %use MAP
%     eval(['MAP1_14(' partfilename ',sfreq,-1,parameterfile,''probability'');']);
%     global ANprobRateOutput
%     eval(['AN_' partfilename ' = ANprobRateOutput(1:size(ANprobRateOutput,1)/2,:);']);
%     
%     %Do the IPIH analysis
%     eval(['[iih_' partfilename ',IPIhisttime_' partfilename ',IPIhistweight_' ...
%         partfilename '] = track_formants_from_IPI_guy(AN_' partfilename ',sfreq);']);
%     eval(['poolIPI_across_channels(IPIhisttime_' partfilename ',IPIhistweight_' partfilename ');']);
%     title(partfilename);
%     %set(gca,'Title',partfilename);
%     xlabel('Interval (ms)')
%     ylabel('Stimulus time (ms)');
%     
%     eval(['map_iih_onto_log(iih_' partfilename ',30,sfreq);']);
%     title(partfilename);
% end



%%%%%% da stimuli with different pitches %%%%%%%%%%%%%%%%%
 
  dawavfiles = { '200ms_da_080Hz.wav','200ms_da_100Hz.wav','200ms_da_120Hz.wav','200ms_da_140Hz.wav', ... 
      '200ms_da_160Hz.wav','200ms_da_180Hz.wav','200ms_da_200Hz.wav','200ms_da_220Hz.wav', ...
      '200ms_da_240Hz.wav', ...
      'noise.wav', ...
      'da_whispered.wav', ...
      };
      
 da_level = 69; %level of OLLO files in dB SPL
for iCounter = 1:length(dawavfiles)
    
    %read the stimuli
    partfilename = ['da' dawavfiles{iCounter}(1:end-4)];
    eval([ partfilename ' = wavread(''..\wavFileStore\' dawavfiles{iCounter} ''');'])
    
    %set level
    eval([partfilename ' = ' partfilename ...
        './sqrt(mean(' partfilename '.^2)).*10^(-(94-' num2str(da_level) ')/20);']);
    
    %use MAP
    eval(['MAP1_14(' partfilename ',sfreq,-1,parameterfile,''probability'');']);
    global ANprobRateOutput
    eval(['AN_' partfilename ' = ANprobRateOutput(1:size(ANprobRateOutput,1)/2,:);']);
    
    %Do the IPIH analysis
    eval(['[iih_' partfilename ',IPIhisttime_' partfilename ',IPIhistweight_' ...
            partfilename '] = IPIHextract(AN_' partfilename ',sfreq);']);%partfilename '] = track_formants_from_IPI_guy(AN_' partfilename ',sfreq);']);%
    eval(['poolIPI_across_channels(IPIhisttime_' partfilename ',IPIhistweight_' partfilename ');']);
    title(partfilename);
    %set(gca,'Title',partfilename);
    xlabel('Interval (ms)')
    ylabel('Stimulus time (ms)');
    
    eval(['map_iih_onto_log(iih_' partfilename ',30,sfreq);']);
    title(partfilename);
end
% 
% 

