function [ DerivedData ] = SortSelfControlTrials( datafile, Behavior, SaccadeData, DerivedData,objNeuroPhys,AnatomicalInfo,SpikeInfo)
%SortSelfControlTrials Summary of this function goes here
%   Detailed explanation goes here
%close all
% datafile='S:\Data\OptimusData\SelfControl\OP140214_4.accdb';
% 
if nargin
    [PATHSTR,NAME,EXT] = fileparts(datafile);
    dotMatFile=fullfile(PATHSTR,[NAME,'.mat']);
    
else
    datafile ='C:\Data\Isidur\IS150313_1.accdb';    
    dotMatFile  ='C:\Data\Isidur\IS150313_1.plx';
%     datafile ='K:\IS141109_1.accdb';
%     plxfile ='K:\IS141109_1.plx';
end
% if exist(dotMatFile,'file') == 2
%     [PATHSTR,NAME,EXT] = fileparts(datafile);
%     [inpath,infile,ext] = fileparts(dotMatFile);
%     % read data from matlab file
%     Plx2matSC(inpath,[infile ext]);   
% else
%     [PATHSTR,NAME,EXT] = fileparts(datafile);
%     [inpath,infile,ext] = fileparts(datafile);
% end
outfile = [NAME,'.mat'];
outpath=PATHSTR;

% read behavioral data from the Orion ACCDB (*.mdb;*.accdb;)
if nargin == 1
    load(datafile,'Behavior', 'SaccadeData', 'DerivedData','objNeuroPhys','AnatomicalInfo','SpikeInfo')
end

%%
lc ='gr';
lcd ='rkbc';
ms ='o*';
% close all
% datafile = 'C:\Data\Optimus\SelfControl\OP140130_2.mat';
% [pathstr, name, ext] = fileparts(datafile)
% load(datafile)
%% column index of trial information
block               =1;     trial               =2;
error_type          =3;     cond                =4;     direction           =5;
choice1             =6;		choice2             =7;     reward              =8;
mind_change         =9;     ss_delay            =10;    ss_reward           =11;
ss_color            =12;    ll_delay            =13;    ll_reward           =14;
ll_color            =15;    max_delay           =16;    tick                =17;
fix_min             =18;    fix_total           =19;    ss_angle            =20;
ss_direction        =21;    ll_angle            =22;    ll_direction        =23;
pattern             =24;    forced_min          =25;    forced_total        =26;
repetitionCtr       = 27;   temperature          =28;
%% get choice fx on normal intertemporal choice trials
% blob->cond = 1
TRIAL_.data = Behavior.TrialInfo.data;
ll_delays = unique(TRIAL_.data(:,ll_delay));
ll_delays(ll_delays==0)=[];
ll_delays(ll_delays> 10000)=[];
xval=(ll_delays(1):10:ll_delays(end))/1000;
directions = unique(TRIAL_.data(:,direction));
rewards = unique([TRIAL_.data(:,ss_reward),TRIAL_.data(:,ll_reward)],'rows');
%%

[rs idx]=sort(diff(rewards,1,2));
rewards = rewards(idx,:);
% lc = jet(128);
% idx =round(linspace(1,size(lc,1),size(rewards,1)));
% lc=lc(idx,:)
lc = 'grb';
IT_trials_=[];
for gg =1:size(rewards,1)
    ChFxNormal_ITC_=[]; ChFxNormal_ITC_withErrors_=[];
    ChFxNormal_ITC_(1:length(ll_delays),1)=ll_delays;
    ChFxNormal_ITC_withErrors_(1:length(ll_delays),1)=ll_delays;
    for ii = 1:length(ll_delays)
        curSS_trls = find(TRIAL_.data(:,cond) == 1 ...      % condition: normal intertemporal choice
            & TRIAL_.data(:,reward) > 0 ...                 % correct trial thus rewarded
            & TRIAL_.data(:,ll_delay) == ll_delays(ii) ...  % current ll_delay
            & TRIAL_.data(:,choice1) == 0 ...               % ss target chosen
            & TRIAL_.data(:,repetitionCtr)==0 ...           % not a repeated trial
            & TRIAL_.data(:,ss_reward)==rewards(gg,1) ...
            & TRIAL_.data(:,ll_reward)==rewards(gg,2) ...
            );
        IT_trials_(1:length(curSS_trls),ii,1)=curSS_trls;
        
        curSS_err_trls = find(TRIAL_.data(:,cond) == 1 ...  % condition: normal intertemporal choice
            & TRIAL_.data(:,reward) == 0 ...                % error trial thus not rewarded
            & TRIAL_.data(:,ll_delay) == ll_delays(ii) ...  % current ll_delay
            & TRIAL_.data(:,choice1) == 0 ...               % ss target chosen
            & TRIAL_.data(:,repetitionCtr)==0 ...           % not a repeated trial
            & TRIAL_.data(:,ss_reward)==rewards(gg,1) ...
            & TRIAL_.data(:,ll_reward)==rewards(gg,2) ...
            );
        IT_error_trials_(1:length(curSS_err_trls),ii,1)=curSS_err_trls;
        
        curLL_trls = find(TRIAL_.data(:,cond) == 1 ...      % condition: normal intertemporal choice
            & TRIAL_.data(:,reward) > 0 ...                 % correct trial thus rewarded
            & TRIAL_.data(:,ll_delay) == ll_delays(ii) ...  % current ll_delay
            & TRIAL_.data(:,choice1) == 1 ...               % ll target chosen
            & TRIAL_.data(:,repetitionCtr)==0 ...           % not a repeated trial
            & TRIAL_.data(:,ss_reward)==rewards(gg,1) ...
            & TRIAL_.data(:,ll_reward)==rewards(gg,2) ...
            );
        IT_trials_(1:length(curLL_trls),ii,2)=curLL_trls;
        
        curLL_err_trls = find(TRIAL_.data(:,cond) == 1 ...  % condition: normal intertemporal choice
            & TRIAL_.data(:,reward) == 0 ...                % error trial thus not rewarded
            & TRIAL_.data(:,ll_delay) == ll_delays(ii) ...  % current ll_delay
            & TRIAL_.data(:,choice1) == 1 ...               % ll target chosen
            & TRIAL_.data(:,repetitionCtr)==0 ...           % not a repeated trial
            & TRIAL_.data(:,ss_reward)==rewards(gg,1) ...
            & TRIAL_.data(:,ll_reward)==rewards(gg,2) ...
            );
        
        IT_error_trials_(1:length(curLL_err_trls),ii,2)=curLL_err_trls;
        
        ChFxNormal_ITC_(ii,2:3) =[length(curSS_trls)/(length(curSS_trls)+length(curLL_trls)) length(curSS_trls)+length(curLL_trls)];
        ChFxNormal_ITC_withErrors_(ii,2:3) =[(length(curSS_trls)+ length(curSS_err_trls)) ...
            /(length(curSS_trls)+length(curLL_trls) + length(curSS_err_trls)+length(curLL_err_trls))...
            length(curSS_trls)+length(curLL_trls) + length(curSS_err_trls)+length(curLL_err_trls)];
    end
        
    %% get choice fx on maintenance trials
    % blob->cond = 3
    
    ChFx_Maintenance_=[];
    for hh=1:2
        switch hh
            case 1
                curChoice = choice1;
            case 2
                curChoice = choice2;
        end
        for ii = 1:length(ll_delays)
            Numerator = sum(TRIAL_.data(:,cond) == 3 ... % condition: normal intertemporal choice
                & TRIAL_.data(:,ll_delay) == ll_delays(ii) ... % current ll_delay
                & TRIAL_.data(:,reward) > 0 ... % correct trial thus rewarded
                & TRIAL_.data(:,curChoice) == 0 ... % ss target chosen
                & TRIAL_.data(:,ss_reward)==rewards(gg,1) ...
                & TRIAL_.data(:,ll_reward)==rewards(gg,2) ...
                );
            
            Denominator = sum(TRIAL_.data(:,cond) == 3 ... % condition: normal intertemporal choice
                & TRIAL_.data(:,reward) > 0 ... % correct trial thus rewarded
                & TRIAL_.data(:,ll_delay) == ll_delays(ii) ... % current ll_delay
                & TRIAL_.data(:,ss_reward)==rewards(gg,1) ...
                & TRIAL_.data(:,ll_reward)==rewards(gg,2) ...
                );
            
            ChFx_Maintenance_(ii,:,hh) =[ ll_delays(ii) Numerator/Denominator Denominator];
        end
    end
    %% maintenance trial types
    % blob->cond = 3
    
    % L-S choice1 =1  choice2=0
    % L-L choice1 =1  choice2=1    
    % S-L choice1 =0  choice2=1
    % S-S choice1 =0  choice2=0
    
    MaintenanceTrls_=[];
    MaintenanceTrlIdx_=[];
    choices = [1 0        
        1 1
        0 1
        0 0];
    legendCell ={'L-S';'L-L';'S-L';'S-S'};
    C=fliplr(['r','b','y','g']); % make a colors list
    for hh=1:size(choices,1)
        for ii = 1:length(ll_delays)
            curTrls =find(TRIAL_.data(:,cond) == 3 ... % condition: normal intertemporal choice
                & TRIAL_.data(:,ll_delay) == ll_delays(ii) ... % current ll_delay
                & TRIAL_.data(:,reward) > 0 ... % correct trial thus rewarded
                & TRIAL_.data(:,choice1) == choices(hh,1) ... % which target chosen initially
                & TRIAL_.data(:,choice2) == choices(hh,2) ... % which target chosen finally
                & TRIAL_.data(:,ss_reward)==rewards(gg,1) ...
                & TRIAL_.data(:,ll_reward)==rewards(gg,2) ...
                );
%             if ~isempty(curTrls)
%                Sacc TargetOn_(:,1)
%             end
            nTrls = sum(TRIAL_.data(:,cond) == 3 ... % condition: normal intertemporal choice
                & TRIAL_.data(:,ll_delay) == ll_delays(ii) ... % current ll_delay
                & TRIAL_.data(:,reward) > 0 ... % correct trial thus rewarded
                & TRIAL_.data(:,choice1) == choices(hh,1) ... % which target chosen initially
                & TRIAL_.data(:,choice2) == choices(hh,2) ... % which target chosen finally
                & TRIAL_.data(:,ss_reward)==rewards(gg,1) ...
                & TRIAL_.data(:,ll_reward)==rewards(gg,2) ...
                );
            % trial types 1st dim
            % ll delays   2nd dim
            MaintenanceTrls_(ii,hh)=nTrls;
            MaintenanceTrlIdx_(1:length(curTrls),ii,hh)=curTrls;
        end
    end
    %% control trial types
    curFSw_trls = find(TRIAL_.data(:,cond) == 4 ... % condition: forced switch
        & TRIAL_.data(:,reward) > 0 ... % correct trial thus rewarded
        & TRIAL_.data(:,ss_delay) > 0 ... % current ll_delay
        );
    ForcedSwitchTrls_(1:length(curFSw_trls),1)=curFSw_trls;
    
    for ii = 1:length(ll_delays)        
        curFSt_trls = find(TRIAL_.data(:,cond) == 5 ... % condition: forced switch
            & TRIAL_.data(:,reward) > 0 ... % correct trial thus rewarded
            & TRIAL_.data(:,ll_delay) == ll_delays(ii) ... % current ll_delay
            );
        ForcedStayTrls_(1:length(curFSt_trls),ii)=curFSt_trls;
    end

end
% add to output variable
DerivedData.ChoiceFunctionNoTemptation          = ChFxNormal_ITC_;  
DerivedData.ChoiceFunctionTemptation            = ChFx_Maintenance_;
DerivedData.ChoiceFunctionNoTemptationErrors    = ChFxNormal_ITC_withErrors_;

DerivedData.NoTemptationTrials                  = IT_trials_;
DerivedData.TemptationTrials                    = MaintenanceTrlIdx_;
DerivedData.ErrorTemptationTrials               = IT_error_trials_;
DerivedData.ForcedSwitchTrials                  = ForcedSwitchTrls_;
DerivedData.ForcedStayTrials                    = ForcedStayTrls_;
DerivedData.MaintenanceTrialIdx                 = MaintenanceTrls_;                                               


% save to .mat file
% save (fullfile(outpath,outfile),'*_','-append')
if exist(fullfile(outpath,outfile),'file')==2
    save(fullfile(outpath,outfile),'DerivedData','-append')
else
    save(fullfile(outpath,outfile),'DerivedData')
end

% SaccadeDetectionPLX(plxfile);

SelfControlSummaryPlot3(datafile, Behavior, SaccadeData, DerivedData,objNeuroPhys,AnatomicalInfo,SpikeInfo)
% SelfControlComparisons(fullfile(outpath,outfile))

%% find saccade of interest
% SaccBegin > TargetOn

%%
end
