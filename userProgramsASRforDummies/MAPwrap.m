% This function wraps up whatever version of MAP I want to call. It is
% implemented partly because I want to avoid messing with jobject too much
% and mostly because I dont want to declare globals in my class.

function [myANprobRateOutput, mydt, myBF] = MAPwrap(stimulus, sampleRate, BFlist, participant, AN_spikesOrProbability, paramChanges)


global ANprobRateOutput  dt BFlist
% disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))
MAP1_14(stimulus, sampleRate, BFlist, participant, AN_spikesOrProbability, paramChanges);
% disp(20*log10(sqrt(mean(stimulus.^2))/20e-6))
myANprobRateOutput   = ANprobRateOutput;
mydt = dt;
myBF = BFlist;
