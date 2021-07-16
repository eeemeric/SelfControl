function [raster,histogram,pt,SDF]= getSDF(spike,trials,event,peritime,KernelType)
%spikeDensityFunction, returns a perievent spikeDensityFunction.
%Input Arguments:
% spike, spike times (in ms) relative to the start of trials for the entire
%   recording (double matrix. ntrials x max spikes in a trial).
%   Will also accept spike times in seconds relative to the start of the
%   recording (double vector. nSpikes x 1)
% trials, index of trials used to produce the spikeDensityFunction (vector
%   of integers).
% event, event times (in ms) relative to the start of trials for the entire
%   recording. Will also accept event times in seconds relative
%   to the start of the recording (double vector. nEvents x 1)
% peritime, the interval of time (in ms) relative to the event from which
%   the spikeDensityFunction will be produced. If spike is a double vector,
%   the interval of time needs to be in seconds. 
% KernelType, the kernel type used to convolve the spike histogram. If 1, a
%   the spikeDensityFunction will be convolved with a PSP filter (see
%   below for specifics). If 2, a half gaussian filter will be used to convolve
%   the spikeDensityFunction.
%
%   PSP filter:
%   Spike density functions were constructed by convolving spike
%   trains with a combination of growth and decay exponential functions that
%   resembled a postsynaptic potential given by the equation
%   R(t)=(1?exp(?t/?g))?(exp(?t/?d))
%   where rate as a function of time [R(t)] varies according to ?g, the time
%   constant for the growth phase, and ?d, the time constant for the decay
%   phase. Physiological data from excitatory synapses indicate that 1 and
%   20 ms are good values for ?g and ?d, respectively (Kim and Connors 1993;
%   Mason et al. 1991; Sayer et al. 1990; Thomson et al. 1993). The rationale
%   for this approach has been described previously (Hanes and Schall 1996;
%   Thompson et al. 1996); its motivation was to derive physiologically
%   plausible spike density functions.
%
%Output Arguments:
% raster, spike times on each trial relative to the event (nTrials x max
%   spikes across trials. NaN padded)
% histogram, histogram (spike counts in millisecond time bins) of the times
%   at which the neuron fires relative to the event.
% pt, a vector of the corresponding time bins of the histogram and
%   spikeDensityFunction (in ms).
% SDF,  smooth and continuous function of spike rate as a function of time
%   relative to the event.
%
% erik.emeric@gmail.com
if size(spike,2) > 1
    % spike times and events in miliseconds relative to the start of each
    % trial
    [raster,histogram,pt,SDF]=spikeDensityFunctionTrialBased(spike,trials,event,peritime,KernelType);
else
    % spike times and events in seconds relative to the start of recording
    [raster,histogram,pt,SDF]=spikeDensityFunction(spike,trials,event,peritime,KernelType);
end
end