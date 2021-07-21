%% read_orion
% read accdb file created by the selfcontrol experiment run on the orion
% platform into matlab
clc
% cla
ACCDBPATH = 'D:\SelfControl\data\SelfControl\AG130820_1.accdb';
[TRIAL,EYE]=read_orion(ACCDBPATH,'trial','eye');

%% read TDT data into MatLab
tankPath = 'D:\SelfControl\data\SelfControl\Aragorn-1';
% get Spikes 
data = TDTbin2mat(tankPath);
%% parse TDT and Orion data
% [NeuroPhys,ornTrialdata]=parseTDTdata_Orion_SC(TDT_data,orionEvent,ornTrialdata);
[NeuroPhys,ornTrialdata]=parseTDTdata_Orion_SC(data,EYE,TRIAL);