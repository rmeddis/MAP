frequencies=[1000 1250];
amplitudesdB=[20 23];
nFrequencies=length(frequencies);

dt=0.0001;

toneDuration=.010;
time=dt:dt:toneDuration;

% fixed ramp, silenceDuration, toneDuration
rampDuration=0.005;
rampTime=dt:dt:rampDuration;
ramp=[0.5*(1+cos(2*pi*rampTime/(2*rampDuration)+pi)) ...
    ones(1,length(time)-length(rampTime))];
ramp=ramp.*fliplr(ramp);

silenceDuration=0.010;
silenceDurationLength=round(silenceDuration/dt);
initialSilence=zeros(1,silenceDurationLength);

silenceToneDuration=toneDuration + silenceDuration;
silenceToneDurationLength=round(silenceToneDuration/dt);

totalDuration=silenceToneDuration*nFrequencies;
totalDurationLength=round(totalDuration/dt);
stimulus=zeros(1,totalDurationLength);
toneBeginPTR=1;

for i=1:nFrequencies
    frequency=frequencies(i);
    dBSPL=amplitudesdB(i);
    amplitude=28e-6* 10.^(dBSPL/20);
    tone=amplitude*sin(2*pi*frequency*time);
    tone=tone.*ramp;
    stimulus(toneBeginPTR:toneBeginPTR+silenceToneDurationLength-1)=...
        [initialSilence tone];    
    toneBeginPTR=toneBeginPTR+silenceToneDurationLength;
end
figure(2), plot( stimulus')
