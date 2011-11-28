function x=UTIL_Butterworth (x, dt, fl, fu, order)
% UTIL_Butterworth (x, dt, fu, fl, order)
% Taken from Yuel and beauchamp page 261
% NB error in their table for K (see their text)
% x is original signal
% x is the filtered output (approx 3 dB down at cutoff for first order filter)
% fu, fl upper and lower cutoff
% order is the number of times the filter is applied
% (approx 6 dB attenuation is corrected)

sampleRate=1/dt;
if 4*fu>sampleRate  
    error(['UTIL_Butterworth: sample rate ' num2str(sampleRate) ' is too low.  Sampling rate should be >' num2str(4*fu)])
end


q=(pi*dt*(fu-fl));
J=1/(1+ cot(q));
K= (2*cos(pi*dt*(fu+fl)))/((1+tan(q))*cos(q));
L= (tan(q)-1)/(tan(q)+1);
b=[J 0 -J];
a=[1 -K  -L];

[numChannels nDataPoints]=size(x);
for channelNumber=1:numChannels
	for j=1:order
		% x2 compensate for 6 dB attenuation
		x(channelNumber,:)=2*filter(b, a, x(channelNumber,:));
	end
end

