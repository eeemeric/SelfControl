clear all
close all
load('D:\SelfControl\dotMat\Aragorn-1.mat')
%%
trials = nonzeros(DerivedData.NoTemptationTrials(:));  
event = DerivedData.fix_cue_on;
peritime = [-1000 1000];
KernelType = 1;
whichSpike=1;
[raster,histogram,pt,SDF] = objNeuroPhys.getSDFRaster(whichSpike,trials,event,peritime,KernelType);

close all
% plot results
subplot(2,1,1)
plot(raster{1},1:size(raster{1},1),'.r')
xlim(peritime)
subplot(2,1,2)
plot(pt{1},SDF{1})
