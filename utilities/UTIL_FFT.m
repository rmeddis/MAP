function [fft_powerdB, fft_phase, frequencies, fft_ampdB]=...
    UTIL_FFT(getFFTsignal, dt, referenceAmplitude)
% UTIL_FFT
% [fft_powerdB, fft_phase, frequencies, fft_ampdB]= UTIL_FFT(getFFTsignal, dt, referenceAmplitude)

x=length(getFFTsignal);         % adjust the length of the signal to a power of 2 for FFT
n=2;
while n<=x
    n=n*2;
end
n=round(n/2);
getFFTsignal=getFFTsignal(1:n);
frequencies = (1/dt)*(1:n/2)/n;

fft_result = fft(getFFTsignal, n);				    % Compute FFT of the input signal.
fft_power = fft_result .* conj(fft_result); % / n;	% Compute power spectrum.  
% Dividing by 'n' we get the power spectral density.
%     fft_power = fft_result .* conj(fft_result) / n;	% Compute power spectralDensity.  

% Prepare the FFT for display
fft_power=fft_power(1:length(fft_power)/2);         % remove mirror frequencies
%     fft_powerdB = UTIL_amp2dB (fft_power, max(fft_power)); % convert to dB
fft_powerdB = 10*log10(fft_power);                  % convert to dB

% amplitude spectrum
if nargin<3
    referenceAmplitude=28e-6; % for SPL
end
% refAmpdB= referenceAmplitude^0.5;
fft_ampdB = 10*log10(fft_power/referenceAmplitude); % convert to dB
%     peak_fft_powerdBSPL=max(fft_ampdB)


fft_phase = angle(fft_result);			            % Compute the phase spectrum.
fft_phase=fft_phase(1:length(fft_phase)/2);         % remove mirror phases
jags=find(diff(fft_phase)>0);                       % unwrap phase
for i=jags, 
    fft_phase(i+1:end)=fft_phase(i+1:end)-2*pi; 
end
