function [raster,histogram,pt,SDF]=spikeDensityFunction(spike,trials,event,peritime,KernelType);
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

plotflag = 0;
raster=[];histogram=[];SDF=[];pt = [];

% if ~nargin
    error('spikeDensityFunction needs to be updated')
% end

if nargin==4
    KernelType=1; % PSP kernel
end
binsize=1; % ms 
%% for every trial ...
nspks = [];
for ii = 1:length(trials)
    %% find all spikes times... subtract the event of interest from spike times
    allTrialSpikes = nonzeros(spike(trials(ii),:))- event(trials(ii),1);
    perieventTrialSpikes = allTrialSpikes(allTrialSpikes>=peritime(1) & allTrialSpikes<=peritime(2));
    %%
    raster(ii,1:length(perieventTrialSpikes))=perieventTrialSpikes;
    nspks(ii,1)=length(perieventTrialSpikes);
end
for ii = 1:length(trials)
    if nspks(ii)<size(raster,2)
        raster(ii,nspks(ii)+1:end)=nan;
    end
end
if isempty(raster)
    raster =nan(size(nspks));
end
pt=peritime(1):binsize:peritime(2);
[histogram,bin]=histc(raster(:),pt);

%%
if KernelType ==1
    Growth=1; Decay=20;
    Half_BW=round(Decay*8);
    BinSize=(Half_BW*2)+1;
    Kernel=[0:Half_BW];
    Half_Kernel=(1-(exp(-(Kernel./Growth)))).*(exp(-(Kernel./Decay)));
    Half_Kernel=Half_Kernel./sum(Half_Kernel);
    Kernel(1:Half_BW)=0;
    Kernel(Half_BW+1:BinSize)=Half_Kernel;
    Kernel=Kernel;
    Kernel=Kernel';
elseif KernelType ==2
    sigma=30;
    X=-50:50;
    Kernel = 1/(sqrt(2*pi)*sigma)*exp(-0.5*X.^2/(sigma^2));
    Kernel=Kernel/sum(Kernel);
    Kernel=Kernel;
    Kernel=Kernel';
    npts =length(Kernel);%npts = 50;
end


%hold on
SDF = (convn(histogram(~isnan(histogram)), Kernel))/length(trials)*1000;
SDF(1:ceil(length(Kernel)/2))=[];
SDF=SDF(1:length(histogram));
SDF= SDF(:);
if plotflag
    subplot(1,2,1)
    plot(pt,SDF,'r')
    hold on
    plot([0 0], ylim)
    
    subplot(1,2,2)
    plot(pt,histogram,'b')
    hold on
    plot([0 0], ylim)
    keyboard
end
end
