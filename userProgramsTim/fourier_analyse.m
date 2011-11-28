function [frequency,out] = fourier_analyse(in, sfreq,color_plot)
%
%plots absolute values of fft of input signal on a semilogarithmic scale
%long term spectrum
%
%use:      [frequency,out] = fourier_analyse(in, sfreq)
%          in: input signal
%          sfreq: samplingfrequency
%          frequency: frequency vector
%          out: fourier-spectrum (complex)
if ~exist('color_plot')
    color_plot = 'b';
end
out = fft(in)/(length(in));  %%Normierung auf Wurzel der Länge, da Matlab intern ohne Normierung der Hin-Fouriertransformation arbeitet
t = [0:1/sfreq:length(in)/sfreq-1/sfreq];
frequency = [0:1/t(end):1/(2*(t(2)-t(1)))];
%spektrale leistungsdichte wird geplottet, wobei eine amplitude von 1 100
%dB entspricht
plot(frequency,20*log10(sqrt(2)*abs(out(1:round(length(in)/2)))),color_plot);
%sqrt(2) weil Gesamtenergie auch in Spiegelfrequenzen enthalten
xlabel('frequency / Hz');
ylabel('fourier amplitude / dB');