# Derive spike density functions

**getSDF**, returns perievent spike density function, raster, and histogram.

e.g., [raster,histogram,pt,SDF]= getSDF(spike,trials,event,peritime,KernelType)

## Input Arguments:
- *spike*, spike times (in ms) relative to the start of trials for the entire recording (double matrix. ntrials x max spikes in a trial). Will also accept spike times in seconds relative to the start of the recording (double vector. nSpikes x 1) trials, index of trials used to produce the spikeDensityFunction (vector of integers).
- *event*, event times (in ms) relative to the start of trials for the entire recording. Will also accept event times in seconds relative to the start of the recording (double vector. nEvents x 1) peritime, the interval of time (in ms) relative to the event from which the spikeDensityFunction will be produced. If spike is a double vector, the interval of time needs to be in seconds. 
- *KernelType*, the kernel type used to convolve the spike histogram. If 1, a the spikeDensityFunction will be convolved with a PSP filter (see
  below for specifics). If 2, a half gaussian filter will be used to convolve the spike density function.

## Output Arguments:
- *raster*, spike times on each trial relative to the event (nTrials x max spikes across trials. NaN padded)
- *histogram*, histogram (spike counts in millisecond time bins) of the times at which the neuron fires relative to the event.
- *pt*, a vector of the corresponding time bins of the histogram and spike density function (in ms).
- *SDF*,  smooth and continuous function of spike rate (spikes/sec) as a function of time relative to the event.
