%%
clc
clear all

%% read_orion
% read accdb file created by the selfcontrol experiment run on the orion
% platform into matlab
ACCDBPATH = 'D:\SelfControl\data\SelfControl\AG130820_1.accdb';
[TRIAL,EYE]=read_orion(ACCDBPATH,'trial','eye');
% TRIAL is a struct that contains the trial parameters
%   TRIAL.len   single value of class double. not sure what this represents
%   TRIAL.ver   single value of class double that indicates which version of the Orion experiment control code was used to collect the data.
%   TRIAL.data  nTrials x up to 28 columns. each row contains the trial parameters for each trial. each column represents a different trial parameter.
%                   ## column index of trial information contained in TRIAL.data ##
%                   block               = 1;     
%                   trial               = 2;
%                   error_type          = 3;     
%                   cond                = 4;     
%                   direction           = 5;
%                   choice1             = 6;	    
%                   choice2             = 7;     
%                   reward              = 8;
%                   mind_change         = 9;     
%                   ss_delay            =10;    
%                   ss_reward           =11;
%                   ss_color            =12;    
%                   ll_delay            =13;    
%                   ll_reward           =14;
%                   ll_color            =15;    
%                   max_delay           =16;    
%                   tick                =17;
%                   fix_min             =18;    
%                   fix_total           =19;    
%                   ss_angle            =20;
%                   ss_direction        =21;    
%                   ll_angle            =22;    
%                   ll_direction        =23;
%                   pattern             =24;    
%                   forced_min          =25;    
%                   forced_total        =26;
%                   repetitionCtr       =27;   
%                   temperature         =28;    
%% read TDT data into MatLab
tankPath = 'D:\SelfControl\data\SelfControl\Aragorn-1';
% read the binary tank data
data = TDTbin2mat(tankPath);
% Data is a struct that contains the following fields
%   data.epocs      contains all epoc store data (onsets, offsets, values)
%   data.snips      contains all snippet store data (timestamps, channels, and raw data)
%   data.streams    contains all continuous data (sampling rate and raw data)
%   data.scalars    contains all scalar data (samples and timestamps)
%   data.info       contains additional information about the block
%% parse TDT and Orion data
% [NeuroPhys,ornTrialdata]=parseTDTdata_Orion_SC(TDT_data,orionEvent,ornTrialdata);
[NeuroPhys,ornTrialdata] = parseTDTdata_Orion_SC(data,EYE,TRIAL);
% NeuroPhys contains the TDT data.
%   struct with fields:
%       NeuroPhys.HeaderInformation: [1×1 struct]
%           struct with fields:
%       	NeuroPhys.HeaderInformation.tankpath        string, directory where the TDT data is stored   
%           NeuroPhys.HeaderInformation.blockname       string, TDT file name
%           NeuroPhys.HeaderInformation.date            string, date data was acquired
%           NeuroPhys.HeaderInformation.utcStartTime    string, start time of the recording
%           NeuroPhys.HeaderInformation.utcStopTime     string, end time of the recording
%           NeuroPhys.HeaderInformation.duration        string, duration of the recording
%           NeuroPhys.HeaderInformation.streamChannel   double, unclear what this means
%           NeuroPhys.HeaderInformation.snipChannel     double, unclear what this means
%           NeuroPhys.HeaderInformation.LFPfs           double, LFP sample rate
%       NeuroPhys.Digital: [1×1 struct]
%           struct with fields:
%           NeuroPhys.Digital.CodeNumbers         1 x nTrials cell array. Contains the event codes sent from Orion to the TDT
%           NeuroPhys.Digital.CodeTimes           1 x nTrials cell array. Contains the times event codes were received by the TDT  
%           NeuroPhys.Digital.SpikeTimes          1 x nTrials cell array. Contains the times spikes were acquired by the TDT
%           NeuroPhys.Digital.SpikeWaveForms      1 x nTrials cell array. Contains the waveforms of the spikes that were acquired by the TDT
%       NeuroPhys.Analog: [1×1 struct]
%           struct with fields:
%           NeuroPhys.LFP   struct with fields ts (timestamps of the samples), AD01,AD02,... LFP samples   