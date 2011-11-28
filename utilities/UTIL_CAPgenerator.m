function wholeNerveCAP  = UTIL_CAPgenerator...
    (ANresponse, dt, BFlist, numANfibers, plotCAP)
% 
% Generates a compound action potential by convolving an impulse repsonse,
% as defined by Mark Chertoff (2004, JASA), with the response of the 
% auditory nerve.
%
%
% -e(-k*time)*SIN(2*PI()*f*time)
%
% mu(t) = e^-kt * sin(omega*t)
% omega = 2 * pi * f
%
%
% Robert T. Ferry
% 01st May 2008
%
%

nChannels=length(BFlist);
[r nSpikeEpochs]=size(ANresponse);

wholeNerveCAP       = [];
channelCAP          = [];
e                   = exp(1);
k                   = 1000;
impulseDuration     = 0.01;
impulseFrequency    = 750;
impulseTime         = dt:dt:impulseDuration;
impulseResponse     = -e.^(-k*impulseTime).*sin(2*pi*impulseFrequency*impulseTime);
impulseResponse=impulseResponse-mean(impulseResponse);

% WholeNerveCAP
ANoutput = sum(ANresponse, 1);
convolvedWholeNerveCAP = conv(ANoutput, impulseResponse(1,:));
% truncate
convolvedWholeNerveCAP=convolvedWholeNerveCAP(1:nSpikeEpochs);

% apply measurement time constant
sampleRate=1/dt;
upperFreq=sampleRate/4;
lowPassCutOff=40;
wholeNerveCAP = UTIL_Butterworth(convolvedWholeNerveCAP, dt, lowPassCutOff, upperFreq, 1);
% or do not do this
% wholeNerveCAP = convolvedWholeNerveCAP;

% Plot output?

CAPtime = dt:dt:dt*length(wholeNerveCAP);

if (plotCAP == 1)
    figure(9)

    subplot(3,1,1), plot(impulseTime, impulseResponse)
    title('Impulse response')
    xlabel('Time (s)'), ylabel('Amplitude (Pa)')
    xlim([0 max(impulseTime)]), ylim([-inf inf])
	
    subplot(3,1,2), plot(CAPtime, wholeNerveCAP)
    title(['AN CAP (whole-nerve)  '   num2str(length(BFlist)) ' channels' num2str(numANfibers) ' fibers/ch'])
    xlabel('Time (s)'), ylabel('Amplitude (Pa)')
    xlim([0 max(CAPtime)]), ylim([-inf inf])   
	
end
