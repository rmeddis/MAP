function out = mellin_trafo(inx,iny)

%This function computes the Mellin Transformation of a one-dimensional
%Signal in analytical terms given as
%
% S(p)=Integral from 0 to infty of s(t) t^(p-1)dt
% Which equals 
% S(c) = Integral from -Infty to Infty s(t)*exp(j*c*log(t))d(log(t))
%
% taken from Irino and Patterson, Speech Communication 2002 Eq. (1) and (3)
% Tim Juergens, September 2011
%

%Resample the Signal onto a log(t) axis
minimalx=min(inx);
maximalx=max(inx);
logarithmicx= exp([log(minimalx):(log(maximalx)-log(minimalx))/length(inx):log(maximalx)]);
logarithmicy= interp1(inx,iny,logarithmicx,'linear','extrap');
%figure, semilogx(logarithmicx,logarithmicy);

%Absolute of the inverse Fourier-Transform of the resampled signal
out = abs(ifft(logarithmicy));
%figure, plot(out(1:40));


%%just to show that the mellin transform results in invariable patterns if
%%the formants are a constant ratio
% frequencies_short = [1:8:8000];
% frequencies_long = [1:10:10000];
% Intervals_long = 1./frequencies_long;
% Intervals_short = 1./frequencies_short;
% contoursin = sin(2*pi*0.00015.*frequencies_long).^2
% figure, plot(frequencies_long,contoursin), hold on, plot(frequencies_short,contoursin,'r')
% xlabel('Frequency (Hz)')
% ylabel('Rate (arb. units)')
% figure, semilogx(Intervals_long,contoursin), hold on, plot(Intervals_short,contoursin,'r')
% xlabel('Interval (s)')
% ylabel('Rate (arb. units)')
% m1 = mellin_trafo(Intervals_long,contoursin);
% m2 = mellin_trafo(Intervals_short,contoursin);
% figure, plot(m1(1:40)), hold on, plot(m2(1:40),'r');
% xlabel('Mellin variable c');
% ylabel('Magnitude');


