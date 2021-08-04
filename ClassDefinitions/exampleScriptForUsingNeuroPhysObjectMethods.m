clc
clear all
% close all
% on my laptop
% load('D:\SelfControl\dotMat\Aragorn-1.mat')
% on my desktop
load('D:\SelfControl\data\dotMat\Aragorn-1.mat')
% the datafile contains a NeuroPhys object 
% Properties for class NeuroPhysObject:
%     HeaderInformation
%     Digital, spikes and events
%     Analog, LFPs
% Methods for class NeuroPhysObject:
%     NeuroPhysObject, constructor. returns an object of class NeuroPhysObject
%     getSpikes, returns list of spike names and the spike timestamps         
%     getAnalog, returns list of LFP names and data
%     getWFs, returns list of spike names and the spike waveforms           
%     getEvents, returns all of the behavioral events and timestamps        
%     getSDFRaster, returns rasters and histograms aligned on an event of interest. 
%
% **** for all of these methods, you can use the help function to get the documentation (F1 line) for the object class. 
% You can do this 2 different ways from the matlab command line
% 1) You can type "help classDefinition.method" (e.g., help objNeuroPhys.getSpikes)
% 2) Once the object is loaded into the workspace, you can type "help
% objectName.method" (e.g., help objNeuroPhys.getSpikes)


%%
trials = nonzeros(DerivedData.NoTemptationTrials(:)); % trials of interest  
event = DerivedData.fix_cue_on; % alignment event
peritime = [-1000 5000]; % time of alignment
KernelType =1; % the filter used convolve the histogram
whichSpike=2; % index of the spike data used

%% get list of spike names and spike timestamps
[spikeList,spikes] = objNeuroPhys.getSpikes;
disp('Spike data in NeuroPhysObject:')
disp(spikeList)

%% get waveforms of spike data from object
close all
[spikeList,WaveForms] = objNeuroPhys.getWFs;
figure
plot(WaveForms{2}')
title('Spike Waveforms')
%% get analog data from object
[adList,AnalogData] = objNeuroPhys.getAnalog;
close all
figure
plot(AnalogData{2}')
%% get events from object
% TODO method getEvents not working 
%[eventValues,eventTimes] = objNeuroPhys.getEvents;
%% get raster, histogram, and SDF from object
[raster,histogram,pt,SDF] = objNeuroPhys.getSDFRaster(whichSpike,trials,event,peritime,KernelType);

% plot results
subplot(2,1,1)
cla
plot(raster{whichSpike},1:size(raster{whichSpike},1),'.r')
hold on
plot([0 0], ylim,'k')
xlim(peritime)
title('Spike Raster')

subplot(2,1,2)
cla
plot(pt{whichSpike},SDF{whichSpike})
hold on
plot([0 0], ylim,'k')
title('Spike density function')

