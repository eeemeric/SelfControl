function []=SelfControlSummaryPlot3(datafile, Behavior, SaccadeData, DerivedData,objNeuroPhys,AnatomicalInfo,SpikeInfo)
% SELFCONTROLSUMMARYPLOT creates a summary plot of the behavior and
% neurophysiology acquired from the Temptation Task
%
% erik.emeric@gmail.com
if ~nargin    
    datafile = 'D:\SelfControl\data\dotMat\Isildur-17.mat';
    load(datafile,'Behavior', 'SaccadeData', 'DerivedData','objNeuroPhys','AnatomicalInfo','SpikeInfo')
end
if nargin ==1        
    load(datafile,'Behavior', 'SaccadeData', 'DerivedData','objNeuroPhys','AnatomicalInfo','SpikeInfo')
end
figurepath = 'D:\SelfControl\figures\BasicActivityPlot';
analysispath = 'D:\SelfControl\analysisResults\BasicSDFData';
% unpack variables
ChFxNormal_ITC_             = DerivedData.ChoiceFunctionNoTemptation;  
ChFx_Maintenance_           = DerivedData.ChoiceFunctionTemptation;
ChFxNormal_ITC_withErrors_  = DerivedData.ChoiceFunctionNoTemptationErrors;

IT_trials_          = DerivedData.NoTemptationTrials;
MaintenanceTrlIdx_  = DerivedData.TemptationTrials;
IT_error_trials_    = DerivedData.ErrorTemptationTrials;
ForcedSwitchTrls_   = DerivedData.ForcedSwitchTrials;
ForcedStayTrls_     = DerivedData.ForcedStayTrials;
MaintenanceTrls_    = DerivedData.MaintenanceTrialIdx;

if ~isempty(DerivedData.targets_on)
    PlxEvents_.TargetOn_        = DerivedData.targets_on;
    TargetOn_                   = DerivedData.targets_on;
else
    PlxEvents_.TargetOn_       = DerivedData.go_on;
    TargetOn_                  = DerivedData.go_on;
end

SaccOfInterest      = [SaccadeData.SaccadeChoice1 SaccadeData.SaccadeChoice2];
PlxEvents_.Fixate_  = DerivedData.fixation_acquired;
PlxEvents_.Reward_  = DerivedData.reward;
TRIAL_.data         = Behavior.TrialInfo;
% 'fixation_acquired'                   }
%     {'fixation_broken'                     }
%     {'reward'                              }
%     {'fix_cue_on'                          }
%     {'targets_on'                          }
%     {'target_chosen'                       }
%     {'target2_on'                          }
%     {'ChoiceFunctionNoTemptation'          }
%     {'ChoiceFunctionNoTemptationWithErrors'}
%     {'NoTemptationTrials'                  }
%     {'NoTemptationErrorTrials'             }
%     {'ChoiceFunctionTemptation'            }
%     {'TemptationTrials' 
load('C:\Users\Emeric\Documents\MATLAB\DelayedReward\SelfControlDatabaseNeurons.mat')
nMinTrls=0;
% plot colors for different trial types
TrlTypeColors= ...
    [34 139   34;       % ss  -  Forest Green
    30  144  255;       % ll  -  Dodger blue
    255   0    0;       % L-S -  Red
    0     0  255;       % L-L -  Blue
    255 215    0;       % S-L -  Gold
    85  107   47;       % S-S -  Dark Olive Green
    205  92   92;       % Forced Switch - Indian Red	
    25   25  112;       % Forced Stay   - Midnight Blue	
    34 139   34;        % ss  -  Forest Green
    30  144  255] ...   % ll  -  Dodger blue
    ./255;

% if ~nargin
%     datafile = 'C:\Data\Isidur\IS150212_1.mat';
% end
[PATHSTR,NAME,EXT] = fileparts(datafile);

close all
%% plotting params
lc ='gr';
lcd ='rkbc';
ms ='o*';
legendCell ={'L{\rightarrow}S';'L{\rightarrow}L';'S{\rightarrow}L';'S{\rightarrow}S'};

%%
% filename='C:\Data\OptimusData\SelfControl\OP140214_4.mat';
% check if a valid .mat file
[PATHSTR,NAME,EXT] = fileparts(datafile);
if isequal(EXT,'.mat')
    disp(datafile)
else
    error('File is not a .mat file')
end
% load(datafile)
if exist('TankList')
    fileIdx = strcmp(NAME,TankList);
    chUnt{fileIdx};
    trials{fileIdx};
end
% find units
% spkIDs=who('sig*');
spkIDs = fieldnames(objNeuroPhys.Digital.SpikeTimes);
% create list of spikse from SpikeInfo.SingleUnits
unitlabel = 'abcdefghjklmnop';
for jj =1:size(SpikeInfo.SingleUnits,1)
   allUnits{jj,1}= ['Sig',num2str(SpikeInfo.SingleUnits(jj,1)),unitlabel(SpikeInfo.SingleUnits(jj,2))];
end
% if ~isempty(spkIDs)
%     if exist(fullfile(PATHSTR,[NAME,'_SpikeInfo.mat']))
%         load(fullfile(PATHSTR,[NAME,'_SpikeInfo.mat']))
%     else
% %         [allUnits,strtStpUnits,strtStpUnitsRewarded]=findSpikesSC(filename);
%         if exist(fullfile(PATHSTR,[NAME,'_SpikeInfo.mat']))==2
%             save(fullfile(PATHSTR,[NAME,'_SpikeInfo.mat']),'allUnits','strtStpUnits','strtStpUnitsRewarded','-append')
%         else
%             save(fullfile(PATHSTR,[NAME,'_SpikeInfo.mat']),'allUnits','strtStpUnits','strtStpUnitsRewarded')
%         end
%     end
%     if ~exist('allUnits')
%         [allUnits,strtStpUnits,strtStpUnitsRewarded]=findSpikesSC(filename);
%         if exist(fullfile(PATHSTR,[NAME,'_SpikeInfo.mat']))==2
%             save(fullfile(PATHSTR,[NAME,'_SpikeInfo.mat']),'allUnits','strtStpUnits','strtStpUnitsRewarded','-append')
%         else
%             save(fullfile(PATHSTR,[NAME,'_SpikeInfo.mat']),'allUnits','strtStpUnits','strtStpUnitsRewarded')
%         end
%     end
% end
%% task parameters
ll_delay            =13; % column index of ll_delay info for each trial
ss_reward           =11; % column index of ss_reward info for each trial
ll_reward           =14; % ""
pattern             =24; % target configuration
direction           =5;  % direction
choice1             =6;  % 0, if ss target chosen. 1, if the ll target chosen
TRIAL_.data=Behavior.TrialInfo.data;
if size(TRIAL_.data,2)==23
    TRIAL_.data(:,24) = SaccadeData.SaccadeDir;
end
directions  = SaccadeData.SaccadeDir; % direction of choice1
ndirections = unique(directions(~isnan(directions)));

% ndirections = unique(TRIAL_.data(:,[pattern direction]),'rows');
% if size(ndirections,1)>4
%     ndirections = unique([TRIAL_.data(:,pattern) SaccadeData.SaccadeDir],'rows');
% end

ll_delays = unique(TRIAL_.data(:,ll_delay));
ll_delays(ll_delays==0)=[];
ll_delays(ll_delays > 10000)=[];
xval=(ll_delays(1):10:ll_delays(end))/1000;
rewards = unique([TRIAL_.data(:,ss_reward),TRIAL_.data(:,ll_reward)],'rows');

%% plot behavior resultsfigure;
subplotrows = 2;
subplotcols = 2;

mp = get(0, 'MonitorPositions');
%%
figure
set(gcf,'position',[-1919           1        1920        1004] ...%[-1079           9        1068        1828]...
    ,'units','inches'...
    ,'paperposition',[0.25 0.25 8   10.5])
%%
orient portrait
for gg =1:size(rewards,1)
    
    annotation('textbox',[.4 .9 .2 .1],'string',NAME,'FontWeight','bold','interpreter','none','LineStyle','none')
    subplot(subplotrows,subplotcols,1)
    cla
    hold on
    % plot lines outside of range for legend
    plot(-1,-1,'-ok')
    plot(-1,-1,'-og')
    plot(-1,-1,'-*r')
    plot(-1,-1,'--sk')
    
    plot(ChFxNormal_ITC_(:,1)/1000,ChFxNormal_ITC_(:,2),'ok',...
        'markersize',6)
    % plot curvefit
    params = glmfit(ChFxNormal_ITC_(:,1)/1000,ChFxNormal_ITC_(:,2),'binomial', 'logit');
    yhat = glmval(params,xval,'logit');
    plot(xval,yhat,'k')
    
    [m index]=min(abs(yhat-.5));
    y = @(x) glmval(params,x,'logit');% your function of interest
    indifferencePoint = fminsearch(@(x) (y(x)-0.5)^2,xval(index)); % find x which makes y=0.5
    
     xlim([0 xval(end)+ .5]); ylim([-.05 1.05])
    plot(xlim,[.5 .5],'--k')
    
    if ~isempty(nonzeros(ChFx_Maintenance_(:,3,:)))
        % curve fit for maintenance trials
        for hh = 1:size(ChFx_Maintenance_,3)
            plot(ChFx_Maintenance_(:,1,hh)/1000,ChFx_Maintenance_(:,2,hh),[lc(hh),ms(hh)],'markersize',7)
            % plot curvefit
            params = glmfit(ChFx_Maintenance_(:,1,hh)/1000,ChFx_Maintenance_(:,2,hh),'binomial', 'logit');
            yhat = glmval(params,xval,'logit');
            plot(xval,yhat,lc(hh),'linewidth',3)
            [m index]=min(abs(yhat-.5));
            y = @(x) glmval(params,x,'logit');% your function of interest
            indifferencePoint(hh+1) = fminsearch(@(x) (y(x)-0.5)^2,xval(index));
        end
    end
    
    if ~isempty(nonzeros(ChFxNormal_ITC_withErrors_(:,3,:)))
        % include err trials
        plot(ChFxNormal_ITC_withErrors_(:,1)/1000,ChFxNormal_ITC_withErrors_(:,2),'sk',...
            'markersize',6)
        % plot curvefit
        params = glmfit(ChFxNormal_ITC_withErrors_(:,1)/1000,ChFxNormal_ITC_withErrors_(:,2),'binomial', 'logit');
        yhat = glmval(params,xval,'logit');
        plot(xval,yhat,'--k')
        
        [m index]=min(abs(yhat-.5));
        y = @(x) glmval(params,x,'logit');% your function of interest
        indifferencePoint(4) = fminsearch(@(x) (y(x)-0.5)^2,xval(index)); % find x which makes y=0.5
    end
    % legend info
    indifferencePoint=ceil(indifferencePoint*100)/100;
    trlType={'Normal','T_{ch1}','T_{ch2}','Normal_{wErrors}'};
    legendstr={};
    for hh =1:length(indifferencePoint)
        legendstr{hh}=[trlType{hh},': ',num2str(indifferencePoint(hh))];
    end
    legend(legendstr,'location','EastOutside')
    legend('boxoff')
    
   
    title({['Rewards -', num2str(rewards(gg,1)),':',num2str(rewards(gg,2))];
        ['No Temptation, ', num2str(sum(ChFxNormal_ITC_(:,3)))];
        ['Temptation, ', num2str(sum(ChFx_Maintenance_(:,3,1)))]})
    xlabel('LL delay (sec)')
    ylabel('P(SS)','fontsize',8)
    set(gca, 'tickdir','out','ytick',0:.25:1)
    
    if ~isempty(nonzeros(MaintenanceTrls_))
        temp = MaintenanceTrls_;
        idx = sum(MaintenanceTrls_,2)>0;
        temp=temp(idx,:); 
        dly = ll_delays(idx);
    % proportion of trial types by delay
    subplot(subplotrows,subplotcols,3);
    bar(temp,'stacked');
    set(gca,'xticklabel',ll_delays/1000)
    lP= legend (legendCell,'Location','EastOutside');
    legend('boxoff')
    axis('tight')
    % set patch colors
    P=findobj(gca,'type','bar');
    for n=1:length(P)
        set(P(n),'facecolor',TrlTypeColors(7-n,:) ...
            ,'edgecolor',TrlTypeColors(7-n,:));
%         set(lP(n),'facecolor',TrlTypeColors(7-n,:) ...
%             ,'edgecolor',TrlTypeColors(7-n,:));
    end
    
    set(gca, 'tickdir','out')
    % xlabel('LL delay (sec)')
    ylabel('N trials','fontsize',8)
    title(['Maintenance, ', num2str(sum(ChFx_Maintenance_(:,3,1)))])
    oa = get(gca,'position');
    set(gca,'position',[oa(1) (oa(2)+oa(4)/2 +.05)  oa(3) (oa(4)/2-.05)])
    
    h = axes('position',[oa(1) oa(2)  oa(3) (oa(4)/2-.05)]);
    sumMT= repmat(sum(temp,2),1,size(temp,2));
    bar(temp./sumMT,'stacked');
    % set patch colors
    P=findobj(gca,'type','bar');
    for n=1:length(P)
        set(P(n),'facecolor',TrlTypeColors(7-n,:) ...
            ,'edgecolor',TrlTypeColors(7-n,:));
    end
    
    set(gca,'xticklabel',dly/1000)
    set(gca, 'tickdir','out')
    xlabel('LL delay (sec)')
    ylabel('Proportion')
    axis('tight')
    end
end

%% plot RT distributions
subplot(subplotrows,subplotcols,4);
cla
trls = nonzeros(IT_trials_(:,:,1));%choose ss
[rts idx] = sort(SaccOfInterest(trls,1) - TargetOn_(trls,1));
rts = rts/1000;
trls=trls(idx);
if any(rts<0.1)
    idx =find(rts<0.1);
    expSaccs = trls(idx);
    rts(idx)=[];
    for ii=1:length(idx)
        [m n] = find(IT_trials_(:,:,1)==trls(idx(ii)));
        IT_trials_(m,n,1)=0;
    end
end
hold on
plot(rts,linspace(1/length(rts),1,length(rts)),'color',TrlTypeColors(1,:),'linewidth',2)

trls = nonzeros(IT_trials_(:,:,2));%choose ll
[rts idx] = sort(SaccOfInterest(trls,1) - TargetOn_(trls,1));
rts = rts/1000;
trls=trls(idx);
if any(rts<0.1)
    idx =find(rts<0.1);
    expSaccs = trls(idx);
    rts(idx)=[];
    for ii=1:length(idx)
        [m n] = find(IT_trials_(:,:,2)==trls(idx(ii)));
        IT_trials_(m,n,2)=0;
    end
end
plot(rts,linspace(1/length(rts),1,length(rts)),'color',TrlTypeColors(2,:),'linewidth',2)

trls = nonzeros(MaintenanceTrlIdx_(:,:,1));%L-S 1st and 2nd sacc
[rts idx] = sort(SaccOfInterest(trls,1) - TargetOn_(trls,1));
rts = rts/1000;
rts2 = sort(SaccOfInterest(trls,2)/1000 -SaccOfInterest(trls,1)/1000);
trls=trls(idx);
if any(rts<0.1)
    idx =find(rts<0.1);
    expSaccs = trls(idx);
    rts(idx)=[];
    for ii=1:length(idx)
        [m n] = find(MaintenanceTrlIdx_(:,:,1)==trls(idx(ii)));
        MaintenanceTrlIdx_(m,n,2)=0;
    end
end
plot(rts,linspace(1/length(rts),1,length(rts)),'color',TrlTypeColors(3,:),'linewidth',2)
plot(rts2,linspace(1/length(rts2),1,length(rts2)),'--','color',TrlTypeColors(3,:),'linewidth',2)
axis([0 1 -.1 1.1])

trls = nonzeros(MaintenanceTrlIdx_(:,:,2));%L-L
[rts idx] = sort(SaccOfInterest(trls,1) - TargetOn_(trls,1));
rts = rts/1000;
trls=trls(idx);
if any(rts<0.1)
    idx =find(rts<0.1);
    expSaccs = trls(idx);
    rts(idx)=[];
    for ii=1:length(idx)
        [m n] = find(MaintenanceTrlIdx_(:,:,2)==trls(idx(ii)));
        MaintenanceTrlIdx_(m,n,2)=0;
    end
end
plot(rts,linspace(1/length(rts),1,length(rts)),'color',TrlTypeColors(4,:),'linewidth',2)

trls = nonzeros(MaintenanceTrlIdx_(:,:,3));%S-L 1st and 2nd sacc
[rts idx] = sort(SaccOfInterest(trls,1) - TargetOn_(trls,1));
rts = rts/1000;
rts2 = sort(SaccOfInterest(trls,2)/1000 -SaccOfInterest(trls,1)/1000);
trls=trls(idx);
if any(rts<0.1)
    idx =find(rts<0.1);
    expSaccs = trls(idx);
    rts(idx)=[];
    for ii=1:length(idx)
        [m n] = find(MaintenanceTrlIdx_(:,:,3)==trls(idx(ii)));
        MaintenanceTrlIdx_(m,n,3)=0;
    end
end
plot(rts,linspace(1/length(rts),1,length(rts)),'color',TrlTypeColors(5,:),'linewidth',2)
plot(rts2,linspace(1/length(rts2),1,length(rts2)),'--','color',TrlTypeColors(5,:),'linewidth',2)

trls = nonzeros(MaintenanceTrlIdx_(:,:,4));%S-S
[rts idx] = sort(SaccOfInterest(trls,1) - TargetOn_(trls,1));
rts = rts/1000;
trls=trls(idx);
if any(rts<0.1)
    idx =find(rts<0.1);
    expSaccs = trls(idx);
    rts(idx)=[];
    for ii=1:length(idx)
        [m n] = find(MaintenanceTrlIdx_(:,:,4)==trls(idx(ii)));
        MaintenanceTrlIdx_(m,n,4)=0;
    end
end
plot(rts,linspace(1/length(rts),1,length(rts)),'--','color',TrlTypeColors(6,:),'linewidth',2)

% legendCell ={'ss';'ll';'L-S 1st';'L-S 2nd';'L-L';'S-L 1st';'S-L 2nd';'S-S'};
% legend (legendCell,'Location','EastOutside');
% legend('boxoff')

set(gca,'tickdir','out', ...
    'ytick',0:.25:1)
xlim([0 .6])
xlabel('Reaction time(sec)')
ylabel('Probability')

%
subplot(subplotrows,subplotcols,2);
cla
hold on
SS_rts_byDelay=[];
LL_rts_byDelay=[];
LS_rts_byDelay =[];
%     % L-S choice1 =1  choice2=0
%     % L-L choice1 =1  choice2=1
%     % S-L choice1 =0  choice2=1
%     % S-S choice1 =0  choice2=0

for curD = 1:size(IT_trials_,2)
    temp = SaccOfInterest(nonzeros(IT_trials_(:,curD,1)),1)-TargetOn_(nonzeros(IT_trials_(:,curD,1)),1);
    
    SS_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
    
    temp = SaccOfInterest(nonzeros(IT_trials_(:,curD,2)),1)-TargetOn_(nonzeros(IT_trials_(:,curD,2)),1);
    LL_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
end

for curD = 1:size(MaintenanceTrlIdx_,2)
    if ~isempty(nonzeros(MaintenanceTrlIdx_(:,curD,1)))
        temp = SaccOfInterest(nonzeros(MaintenanceTrlIdx_(:,curD,1)),2)-SaccOfInterest(nonzeros(MaintenanceTrlIdx_(:,curD,1)),1);       
        LS_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
        
        temp = SaccOfInterest(nonzeros(MaintenanceTrlIdx_(:,curD,1)),1)-TargetOn_(nonzeros(MaintenanceTrlIdx_(:,curD,1)),1);
        LS_1stCh_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
    else
        LS_rts_byDelay(curD,1:2) = nan;
        LS_1stCh_rts_byDelay(curD,1:2)=nan;
    end
    
    if ~isempty(nonzeros(MaintenanceTrlIdx_(:,curD,3)))% S-L
        temp = SaccOfInterest(nonzeros(MaintenanceTrlIdx_(:,curD,3)),2)-SaccOfInterest(nonzeros(MaintenanceTrlIdx_(:,curD,3)),1);
        SL_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
        
        temp = SaccOfInterest(nonzeros(MaintenanceTrlIdx_(:,curD,3)),1)-TargetOn_(nonzeros(MaintenanceTrlIdx_(:,curD,3)),1);
        SL_1stCh_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
    else
        SL_rts_byDelay(curD,1:2) = nan;
        SL_1stCh_rts_byDelay(curD,1:2)=nan;
    end
    
    if ~isempty(nonzeros(MaintenanceTrlIdx_(:,curD,2)))%L-L
        temp = SaccOfInterest(nonzeros(MaintenanceTrlIdx_(:,curD,2)),1)-TargetOn_(nonzeros(MaintenanceTrlIdx_(:,curD,2)),1);
        LL_1stCh_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
    else
        LL_1stCh_rts_byDelay(curD,1:2)=nan;
    end
    
    if ~isempty(nonzeros(MaintenanceTrlIdx_(:,curD,4)))%S-S
        temp = SaccOfInterest(nonzeros(MaintenanceTrlIdx_(:,curD,4)),1)-TargetOn_(nonzeros(MaintenanceTrlIdx_(:,curD,4)),1);
        SS_1stCh_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
    else
        SS_1stCh_rts_byDelay(curD,1:2)=nan;
    end
end

temp = ChFxNormal_ITC_;
idx = ChFxNormal_ITC_(:,3)>0;
temp=temp(idx,:);
SS_rts_byDelay = SS_rts_byDelay(idx,:);

errorbar(temp(:,1)/1000,SS_rts_byDelay(:,1),SS_rts_byDelay(:,2),'-o',...
    'color',TrlTypeColors(1,:),'linewidth',2)

hold on

LL_rts_byDelay=LL_rts_byDelay(idx,:);

errorbar(temp(:,1)/1000,LL_rts_byDelay(:,1),LL_rts_byDelay(:,2),'-o',...
    'color',TrlTypeColors(2,:),'linewidth',2)

temp = ChFx_Maintenance_;
idx = ChFx_Maintenance_(:,3,:)>0;
idx = idx(:,1);
temp=temp(idx,:,:);
dly = ll_delays(idx);

LS_1stCh_rts_byDelay = LS_1stCh_rts_byDelay(idx,:);
errorbar(temp(:,1)/1000,LS_1stCh_rts_byDelay(:,1),LS_1stCh_rts_byDelay(:,2),...
    '-o','color',TrlTypeColors(3,:),'linewidth',2)

LS_rts_byDelay=LS_rts_byDelay(idx,:);
errorbar(temp(:,1)/1000,LS_rts_byDelay(:,1),LS_rts_byDelay(:,2),...
    '-.s','color',TrlTypeColors(3,:),'linewidth',2)

LL_1stCh_rts_byDelay=LL_1stCh_rts_byDelay(idx,:);
errorbar(temp(:,1)/1000,LL_1stCh_rts_byDelay(:,1),LL_1stCh_rts_byDelay(:,2),'-.o',...
    'color',TrlTypeColors(4,:),'linewidth',2)

SL_1stCh_rts_byDelay=SL_1stCh_rts_byDelay(idx,:);
errorbar(temp(:,1)/1000,SL_1stCh_rts_byDelay(:,1),SL_1stCh_rts_byDelay(:,2),'-o',...
    'color',TrlTypeColors(5,:),'linewidth',2)

SL_rts_byDelay=SL_rts_byDelay(idx,:);
errorbar(temp(:,1)/1000,SL_rts_byDelay(:,1),SL_rts_byDelay(:,2),'-.s',...
    'markeredgecolor',TrlTypeColors(5,:),'color',TrlTypeColors(5,:),'linewidth',2)

SS_1stCh_rts_byDelay=SS_1stCh_rts_byDelay(idx,:);
errorbar(temp(:,1)/1000,SS_1stCh_rts_byDelay(:,1),SS_1stCh_rts_byDelay(:,2),'-.o',...
    'color',TrlTypeColors(6,:),'linewidth',2)



legend({'ss','ll','L{\rightarrow}S 1st','L{\rightarrow}S 2nd','L{\rightarrow}L 1st','S{\rightarrow}L 1st','S{\rightarrow}L 2nd','S{\rightarrow}S 1st',},'Location','EastOutside')
legend('boxoff')
set(gca,'tickdir','out')
xlim([0 max(ll_delays/1000)+0.5])
xlabel('ll delay (sec)')
ylabel('RT (ms)')
%%
annotation('textbox',[.4 .9 .2 .1],'string',NAME,'FontWeight','bold','interpreter','none','LineStyle','none')

print(fullfile(figurepath,[NAME,'_BasicActivityPlot.ps']),'-dpsc2','-r300')

clf

if isempty(spkIDs)
    return
end
evts ={'fixate';'target';'saccade';'reward'};
plotxlims = [-.2 1;
    -.2  1;
    -.2  1;
    -.8 .4]*1000;
subplotrows = length(ll_delays)+1;
subplotcols = 4;
for cSpk = 1:size(allUnits) 
    allTrialTypes   = {};
    allSDF          = {};
    allPeriTime     = {};
    allRaster       = {};
    allHistogram    = {};
    allEventLabels  = {};
    
    [~,curSpk] = getSpikes(objNeuroPhys, allUnits(cSpk));
    curSpk = curSpk{1};
    for curEvt = 1:length(evts)
%         SDF = cell(length(ll_delays)+1,size(ndirections,1)+1,10);
%         Marker  = SDF;
%         raster  = SDF;
%         pt      = SDF;
        switch curEvt
            case 1
                evt = PlxEvents_.Fixate_;
                evtLabel{curEvt,1} = 'fixate';;
            case 2
                evt = PlxEvents_.TargetOn_;
                evtLabel{curEvt,1} = 'targets on';
            case 3
                evt = SaccOfInterest;
                evtLabel{curEvt,1} = 'saccade';
            case 4
                evt = PlxEvents_.Reward_;
                evtLabel{curEvt,1} = 'reward';
        end       
        
        maxSDF=0;
%         strtTrl = strtStpUnits(cSpk,1);
%         endTrl = strtStpUnits(cSpk,2);
        strtTrl = SpikeInfo.TrialsPresent(cSpk,1);
        endTrl = SpikeInfo.TrialsPresent(cSpk,2);
        if endTrl > size(TRIAL_.data,1)
            endTrl = size(TRIAL_.data,1);
        end
        
%         ndirections = unique(TRIAL_.data(strtTrl:endTrl,[pattern direction]),'rows');
%         if any(isnan(ndirections(:)))
%            ndirections(isnan(ndirections(:,1) ),:)=[];
%            ndirections = unique(ndirections,'rows');
%         end
%         if size(ndirections,1)>4
%             ndirections = unique([TRIAL_.data(:,pattern) SaccadeData.SaccadeDir],'rows');
%         end
        % get chfxs for the time the cell was present
        tempITtrls=IT_trials_;
        tempMaintTrls = MaintenanceTrlIdx_;
        tempIT_err_trls=IT_error_trials_;
        
        tempITtrls(tempITtrls<strtTrl)=0;
        tempMaintTrls(tempMaintTrls<strtTrl)=0;
        tempIT_err_trls(tempIT_err_trls<strtTrl)=0;
        tempITtrls(tempITtrls>endTrl)=0;
        tempMaintTrls(tempMaintTrls>endTrl)=0;
        tempIT_err_trls(tempIT_err_trls>endTrl)=0;
        
        tempFcdSwTrls = ForcedSwitchTrls_;
        tempFcdStTrls = ForcedStayTrls_;
        
        tempFcdSwTrls(tempFcdSwTrls<strtTrl)=0;
        tempFcdStTrls(tempFcdStTrls<strtTrl)=0;
        tempFcdSwTrls(tempFcdSwTrls>endTrl)=0;
        tempFcdSwTrls(tempFcdSwTrls>endTrl)=0;
        
        for ii=1:size(tempITtrls,2)
            %% normal IT trials
            totTrls = length(nonzeros(tempITtrls(:,ii,1)))+length(nonzeros(tempITtrls(:,ii,2)));
            pSS = length(nonzeros(tempITtrls(:,ii,1)))/totTrls;
            tempChFxIT(ii,:)=[ChFxNormal_ITC_(ii,1)/1000 pSS totTrls];
        end
        for ii =1:size(tempMaintTrls,2)
%             ChFx_Maintenance_%
            totTrls = length(nonzeros(tempMaintTrls(:,ii,:)));
            pSS1 = length(nonzeros(tempMaintTrls(:,ii,3:4)))/totTrls;
            pSS2 = length(nonzeros(tempMaintTrls(:,ii,[1 4])))/totTrls;
            tempChFxTMPT(ii,:,1)=[ChFx_Maintenance_(ii,1)/1000 pSS1 totTrls];
            tempChFxTMPT(ii,:,2)=[ChFx_Maintenance_(ii,1)/1000 pSS2 totTrls];
        end
        %%
        for ii =1:size(tempMaintTrls,3)
            for hh = 1:size(tempMaintTrls,2)
                nTrls = length(nonzeros(tempMaintTrls(:,hh,ii)));
                tempNMaintTrls_(hh,ii)=nTrls;
            end
        end
        
        %%
        for curDelay=1:length(ll_delays)+1
            for curDir = 1:size(ndirections,1)+1
                for curTrlType=1:10
                    switch curTrlType
                        case 1 % ss
                            if curDelay == length(ll_delays)+1
                                trls = nonzeros(tempITtrls(:,:,1));
                            else
                                trls = nonzeros(tempITtrls(:,curDelay,1));
                            end
                            trlTypeFieldName{curTrlType,1}= 'sooner smaller';
                        case 2 % ll
                            if curDelay == length(ll_delays)+1
                                trls = nonzeros(tempITtrls(:,:,2));
                            else
                                trls = nonzeros(tempITtrls(:,curDelay,2));
                            end
                            trlTypeFieldName{curTrlType,1}= 'larger later';
                        case 3 % L-S
                            if curDelay == length(ll_delays)+1
                                trls = nonzeros(tempMaintTrls(:,:,1));
                            else
                                trls = nonzeros(tempMaintTrls(:,curDelay,1));
                            end
                            trlTypeFieldName{curTrlType,1}= 'switch L-S';
                        case 4 % L-L
                            if curDelay == length(ll_delays)+1
                                trls = nonzeros(tempMaintTrls(:,:,2));
                            else
                                trls = nonzeros(tempMaintTrls(:,curDelay,2));
                            end
                            trlTypeFieldName{curTrlType,1}= 'stay L-L';
                        case 5 % S-L
                            if curDelay == length(ll_delays)+1
                                trls = nonzeros(tempMaintTrls(:,:,3));
                            else
                                trls = nonzeros(tempMaintTrls(:,curDelay,3));
                            end
                            trlTypeFieldName{curTrlType,1}= 'switch S-L';
                        case 6 % S-S
                            if curDelay == length(ll_delays)+1
                                trls = nonzeros(tempMaintTrls(:,:,4));
                            else
                                trls = nonzeros(tempMaintTrls(:,curDelay,4));
                            end
                            trlTypeFieldName{curTrlType,1}= 'stay S-S';
                        case 7 % Forced Switch
                            if curDelay == length(ll_delays)+1
                                
                            else
                                trls = nonzeros(tempFcdSwTrls(:));
                            end
                            trlTypeFieldName{curTrlType,1}= 'forced switch';
                        case 8 % Forced Stay
                            if curDelay == length(ll_delays)+1
                                trls = nonzeros(tempFcdStTrls(:));
                            else
                                trls = nonzeros(tempFcdStTrls(:,curDelay));
                            end
                            trlTypeFieldName{curTrlType,1}= 'forced stay';
                        case 9 % error ss trials
                            if curDelay == length(ll_delays)+1
                                trls = nonzeros(tempIT_err_trls(:,:,1));
                            else
                                trls = nonzeros(tempIT_err_trls(:,curDelay,1));
                            end
                            trlTypeFieldName{curTrlType,1}= 'error ss';
                        case 10 % error ll trials
                            if curDelay == length(ll_delays)+1
                                trls = nonzeros(tempIT_err_trls(:,:,2));
                            else
                                trls = nonzeros(tempIT_err_trls(:,curDelay,2));
                            end
                            trlTypeFieldName{curTrlType,1}= 'error ll';
                    end
                    
                    if isempty(trls)
                        continue                        
                    end
                    % find trials from the current direction
                    if curDir < size(ndirections,1)+1                                                    
                        if (curTrlType==1 | curTrlType==5 | curTrlType==6 | curTrlType==9) % initial choice is the ss target
                            idx = directions(trls) == ndirections(curDir)...
                                & TRIAL_.data(trls,choice1)==0;
                        elseif (curTrlType==2 | curTrlType==3 | curTrlType==4 | curTrlType==7 | curTrlType==8 | curTrlType==10 )
                            % initial choice is the ll target
                            idx =  directions(trls) == ndirections(curDir)...
                                & TRIAL_.data(trls,choice1)==1;
                        end
                    elseif curDir == size(ndirections,1)+1
                        % collapse across target directions
                        idx=trls>0;
                    end
                    if ~any(idx)
                        continue
                    end
                    trls = trls(idx);
                    trls(trls>size(PlxEvents_.TargetOn_,1))=[];
                    if length(trls)>nMinTrls                        
                        % add markers
                        mrkr=[];
                        srtIdx=[];
                        switch curEvt
                            case 1 % fixation
                                mrkr = PlxEvents_.TargetOn_(trls,1) - evt(trls,1);
                                [Y, srtIdx] = sort(mrkr(:,1));
                                mrkr = mrkr(srtIdx,:);
                                trls = trls(srtIdx);
                            case 2 % target
                                mrkr(:,1) = PlxEvents_.Fixate_(trls,1) - evt(trls,1);
                                mrkr(:,2) = SaccOfInterest(trls,1) - evt(trls,1);
                                [Y, srtIdx] = sort(mrkr(:,2));
                                mrkr = mrkr(srtIdx,:);
                                trls = trls(srtIdx);
                            case 3 % saccade
                                mrkr(:,1) = PlxEvents_.TargetOn_(trls,1) - evt(trls,1);
                                mrkr(:,2) = SaccOfInterest(trls,2) - evt(trls,1);
                                [Y, srtIdx] = sort(mrkr(:,1));
                                mrkr = mrkr(srtIdx,:);
                                trls = trls(srtIdx);
                            case 4 % reward
                                if curTrlType<9
                                    mrkr(:,1) = SaccOfInterest(trls,1) - evt(trls,1);
                                    mrkr(:,2) = SaccOfInterest(trls,2) - evt(trls,1);
                                    [Y, srtIdx] = sort(mrkr(:,1));
                                    mrkr = mrkr(srtIdx,:);
                                    trls = trls(srtIdx);
                                else
                                    % should make this aligned on the
                                    % corrective saccade but need to check
                                    % if it's a fixation break or an
                                    % attempt to switch
                                    continue
                                end
                        end
                        Marker{curDelay,curDir,curTrlType}=mrkr;
                        [raster{curDelay,curDir,curTrlType},... 
                            histogram{curDelay,curDir,curTrlType},...
                            pt{curDelay,curDir,curTrlType},...
                            SDF{curDelay,curDir,curTrlType}]...
                            = spikeDensityFunctionTrialBased(curSpk,trls,evt,[-1500 1500],2);
                        
                        idx = pt{curDelay,curDir,curTrlType} >= plotxlims(curEvt,1)...
                            & pt{curDelay,curDir,curTrlType} <= plotxlims(curEvt,2);
                        curMax = max(SDF{curDelay,curDir,curTrlType}(idx));
                        if curMax > maxSDF & size(raster{curDelay,curDir,curTrlType},1)>10
                            maxSDF = curMax;
                        end
                    end
                end
            end
        end
        
        % plot the physiology
%         check for  empty columns
        plotswdata=zeros(size(SDF,1),size(SDF,2));
        for mm =1:size(SDF,1)
            for nn = 1:size(SDF,2)
                for oo = 1:size(SDF,3)
                    if ~isempty(SDF{mm,nn,oo})
                        plotswdata(mm,nn)=1;
                    end
                end
            end
        end
        SDF(:,~sum(plotswdata),:)=[];
        pt(:,~sum(plotswdata),:)=[];
        raster(:,~sum(plotswdata),:)=[];
        Marker(:,~sum(plotswdata),:)=[];
        
        if size(SDF,2)-1 <= 0
            continue
        end
        
         [ MainFigHandle,ActivityPlotsHandles,BehaviorPlotsHandles ] = ...
            SC_plot_figure( size(SDF,1)-1,size(SDF,2)-1 );

        
        for mm =1:size(SDF,1)
            for nn = 1:size(SDF,2)
%                 subplot(size(SDF,1)+1,size(SDF,2),mm*size(SDF,2)+nn)
                axes(ActivityPlotsHandles(mm,nn))
                hold on
                periT={};
                pltSDF={};
                pltRstr={};
                for oo = 1:size(SDF,3)
                    periT{oo,1}  =pt{mm,nn,oo};
                    pltSDF{oo,1} =SDF{mm,nn,oo};
                    pltRstr{oo,1}=raster{mm,nn,oo};
                    pltMrkr{oo,1}=Marker{mm,nn,oo};
                end
                PlotSDFandRaster(periT,pltSDF,pltRstr,maxSDF,plotxlims(curEvt,:)./1000,TrlTypeColors,pltMrkr)
                
                if mm == 1 && nn < size(SDF,2)
                    title(['Target ', num2str(nn)])
                elseif mm == 1 && nn==size(SDF,2)
                    title('Collapsed')
                end
                if mm < size(SDF,1)
                   set(gca,'xticklabel',[]) 
                end
                
                if nn == 1
                    if mm < size(SDF,1)
                        ylabel(['Delay: ',num2str(ll_delays(mm)),' ms'],'fontsize',6) 
                    else
                         ylabel('All Delays','fontsize',6)
                    end
                else 
                   set(gca,'yticklabel',[]) 
                end
                
%                 legendstr= [];
%                 for oo = 1:size(SDF,3)
%                     if ~isempty(SDF{mm,nn,oo})
%                         plot(pt{mm,nn,oo},SDF{mm,nn,oo},'color',TrlTypeColors(oo,:))
%                         legendstr=[legendstr;size(raster{mm,nn,oo},1)];
%                     end
%                 end
%                 axis([plotxlims(curEvt,:) 0 maxSDF])
                set(gca,'tickdir','out','ticklength',[.025 .025])
%                 legend(num2str(legendstr),'location','eastoutside')
%                 legend('boxoff')
            end
        end
        xlabel(['Time from ',evts{curEvt},' (ms)'],'fontsize',12)

         %%
%         subplot(subplotrows,subplotcols,1)
        axes(BehaviorPlotsHandles(1))        
        cla
        plot(tempChFxIT(:,1),tempChFxIT(:,2),'ko')
        hold on
        params = glmfit(tempChFxIT(:,1),tempChFxIT(:,2),'binomial', 'logit');
        yhat = glmval(params,xval,'logit');
        plot(xval,yhat,'k')
        
        [m index]=min(abs(yhat-.5));
        y = @(x) glmval(params,x,'logit');% your function of interest
        indifferencePoint = fminsearch(@(x) (y(x)-0.5)^2,xval(index)); % find x which makes y=0.5
        title({['Rewards -', num2str(rewards(gg,1)),':',num2str(rewards(gg,2))];
            ['No Temptation, ', num2str(sum(tempChFxIT(:,3)))];
            ['Temptation, ', num2str(sum(tempChFxTMPT(:,3,1)))]})
        % pSS 1st choice
        plot(tempChFxTMPT(:,1,1),tempChFxTMPT(:,2,1),'go')
        params = glmfit(tempChFxTMPT(:,1,1),tempChFxTMPT(:,2,1),'binomial', 'logit');
        yhat = glmval(params,xval,'logit');
        plot(xval,yhat,'g')
        
        % pSS second choice
        plot(tempChFxTMPT(:,1,2),tempChFxTMPT(:,2,2),'ro')
        params = glmfit(tempChFxTMPT(:,1,2),tempChFxTMPT(:,2,2),'binomial', 'logit');
        yhat = glmval(params,xval,'logit');
        plot(xval,yhat,'r')
        xlabel('LL delay (sec)','fontsize',6)
        ylabel('P(SS)','fontsize',6)
        ylim([-.1 1.1])
        pos = get(gca,'position');
        
        
        yticklabels{1} = sprintf('%1.1f',0);
        yticklabels{2} = sprintf('%1.1f',0.5);
        yticklabels{3} = sprintf('%1.1f',1);
        set(gca, 'tickdir','out', ...
            'ytick',0:.5:1 ...
            ,'fontsize',8 ...
            ,'yticklabels',yticklabels)
        % proportion of trial types by delay
        axes(BehaviorPlotsHandles(3))   
        bar(tempNMaintTrls_,'stacked');
        set(gca,'xticklabel',ll_delays/1000,'fontsize',8)
        
        % set patch colors
        P=findobj(gca,'type','bar');
        for n=1:length(P)
            set(P(n),'facecolor',TrlTypeColors(7-n,:) ...
                ,'edgecolor',TrlTypeColors(7-n,:));
        end
        
                     
        set(gca,'xticklabel',ll_delays/1000, 'tickdir','out')
        xlim([0 length(ll_delays)+1])
        xlabel('LL delay (sec)','fontsize',8)
        ylabel('N trials','fontsize',8)
        title(['Maintenance, ', num2str(sum(tempNMaintTrls_(:)))])
        axis('tight')
        
        % plot RTs
        axes(BehaviorPlotsHandles(2)) 
        cla
        hold on
        LS_1stCh_rts_byDelay=[];
        LS_rts_byDelay=[];
        SL_1stCh_rts_byDelay=[];
        SL_rts_byDelay=[];
        LL_1stCh_rts_byDelay=[];
        SS_1stCh_rts_byDelay=[];
        SS_rts_byDelay=[];
        LL_rts_byDelay=[];
        %     % L-S choice1 =1  choice2=0
        %     % L-L choice1 =1  choice2=1
        %     % S-L choice1 =0  choice2=1
        %     % S-S choice1 =0  choice2=0
        
        for curD = 1:length(ll_delays)
            temp = SaccOfInterest(nonzeros(tempITtrls(:,curD,1)),1)-TargetOn_(nonzeros(tempITtrls(:,curD,1)),1);
            
            SS_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
            
            temp = SaccOfInterest(nonzeros(tempITtrls(:,curD,2)),1)-TargetOn_(nonzeros(tempITtrls(:,curD,2)),1);
            LL_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
            
            if ~isempty(nonzeros(tempMaintTrls(:,curD,1)))
                temp = SaccOfInterest(nonzeros(tempMaintTrls(:,curD,1)),2)-SaccOfInterest(nonzeros(tempMaintTrls(:,curD,1)),1);
                %             if ~isempty(temp)
                %                 plot(curD,temp,'r.')
                %             end
                LS_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
                
                temp = SaccOfInterest(nonzeros(tempMaintTrls(:,curD,1)),1)-TargetOn_(nonzeros(tempMaintTrls(:,curD,1)),1);
                LS_1stCh_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
            else
                LS_rts_byDelay(curD,1:2) = nan;
                LS_1stCh_rts_byDelay(curD,1:2)=nan;
            end
            
            if ~isempty(nonzeros(tempMaintTrls(:,curD,3)))% S-L
                temp = SaccOfInterest(nonzeros(tempMaintTrls(:,curD,3)),2)-SaccOfInterest(nonzeros(tempMaintTrls(:,curD,3)),1);
                SL_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
                
                temp = SaccOfInterest(nonzeros(tempMaintTrls(:,curD,3)),1)-TargetOn_(nonzeros(tempMaintTrls(:,curD,3)),1);
                SL_1stCh_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
            else
                SL_rts_byDelay(curD,1:2) = nan;
                SL_1stCh_rts_byDelay(curD,1:2)=nan;
            end
            
            if ~isempty(nonzeros(tempMaintTrls(:,curD,2)))%L-L
                temp = SaccOfInterest(nonzeros(tempMaintTrls(:,curD,2)),1)-TargetOn_(nonzeros(tempMaintTrls(:,curD,2)),1);
                LL_1stCh_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
            else
                LL_1stCh_rts_byDelay(curD,1:2)=nan;
            end
            
            if ~isempty(nonzeros(tempMaintTrls(:,curD,4)))%S-S
                temp = SaccOfInterest(nonzeros(tempMaintTrls(:,curD,4)),1)-TargetOn_(nonzeros(tempMaintTrls(:,curD,4)),1);
                SS_1stCh_rts_byDelay(curD,1:2)=[nanmean(temp) nanstd(temp)/sqrt(sum(~isnan(temp)))];
            else
                SS_1stCh_rts_byDelay(curD,1:2)=nan;
            end
        end
        
        errorbar(ll_delays/1000,SS_rts_byDelay(:,1),SS_rts_byDelay(:,2),'-o',...
            'color',TrlTypeColors(1,:))
        
        hold on
        
        errorbar(ll_delays/1000,LL_rts_byDelay(:,1),LL_rts_byDelay(:,2),'-o',...
            'color',TrlTypeColors(2,:))
        
        errorbar(ll_delays/1000,LS_1stCh_rts_byDelay(:,1),LS_1stCh_rts_byDelay(:,2),...
            '-o','color',TrlTypeColors(3,:))
        errorbar(ll_delays/1000,LS_rts_byDelay(:,1),LS_rts_byDelay(:,2),...
            '-.s','color',TrlTypeColors(3,:))
        
        errorbar(ll_delays/1000,LL_1stCh_rts_byDelay(:,1),LL_1stCh_rts_byDelay(:,2),'-.o',...
            'color',TrlTypeColors(4,:))
        
        errorbar(ll_delays/1000,SL_1stCh_rts_byDelay(:,1),SL_1stCh_rts_byDelay(:,2),'-o',...
            'color',TrlTypeColors(5,:))
        errorbar(ll_delays/1000,SL_rts_byDelay(:,1),SL_rts_byDelay(:,2),'-.s',...
            'markeredgecolor',TrlTypeColors(5,:),'color',TrlTypeColors(5,:))
        
        errorbar(ll_delays/1000,SS_1stCh_rts_byDelay(:,1),SS_1stCh_rts_byDelay(:,2),'-.o',...
            'color',TrlTypeColors(6,:))
        %     axis('tight')
        set(gca,'tickdir','out','ytick',100:100:600,'fontsize',8)
        axis([0 max(ll_delays/1000)+0.5 100 600])
        xlabel('ll delay (sec)','fontsize',8)
        ylabel('RT (ms)','fontsize',8)
        
        axes(BehaviorPlotsHandles(4)) 
        hold on
        for ii =1:size(TrlTypeColors,1)
            plot(-1,-1, 'color',TrlTypeColors(ii,:),'linewidth',3)
        end
        axis([0 1 0 1])
        lh = legend({'ss','ll','L{\rightarrow}S','L{\rightarrow}L','S{\rightarrow}L','S{\rightarrow}S','Forced Switch','Forced Stay'},'Location','westOutside');
%         legend('boxoff')
        pos = get(lh,'position');
%         pos(1)=pos(1)+.08;
%         pos(2)=pos(2)+.07;        
        set(lh,'position',pos...
            ,'fontsize',8)
        axis off
                
        %%
        axes(BehaviorPlotsHandles(5)) 
        cla
        hold on
        wt=-30;
        wellradius = 28.87/2;
        welltheta = wt/180*pi;
        ct=1;
        for wellrot=0:60:360
            xl(ct,1)= cos(wellrot/180*pi)*wellradius;
            yl(ct,1)= sin(wellrot/180*pi)*wellradius;
            ct=ct+1;
        end
        plot(xl,yl,'r')
        % get location and channel info
        loc_idx = str2num(allUnits{cSpk}(4));
        DepthInfo = AnatomicalInfo.PenetrationInfo.Depth(cSpk) ...
            - AnatomicalInfo.PenetrationInfo.TONA(loc_idx);
        %%
%         AnatomicalInfo,SpikeInfo
        elecX = AnatomicalInfo.PenetrationInfo.XY(:,1);
        elecY = AnatomicalInfo.PenetrationInfo.XY(:,2);
        %%
%         [elecX elecY]=getElectrodeLocation (plug_number,plug_rotation,proboscis_rotation,...
%             electrode_location(loc_idx),0);

        plot(elecX,elecY,'r.','markersize',8)
        
        axis([-15, 30, -15 30])
        set(gca,'plotboxaspectratio',[1 1 1],...
            'xtick',-15:5:15,...
            'ytick',-15:5:15, ...
            'fontsize',8)
        % penetration info
        text(-14,27,'Penetration Info','fontweight','bold')
        text(-14,22,['(',num2str(round(elecX*100)/100),', ',...
            num2str(round(elecY*100)/100),')'],'fontweight','bold')
        text(-14,15,['Depth: ',num2str(DepthInfo),' {\mu}m'],'fontweight','bold')
        grid on
        %%
        annotation('textbox',[.4 .9 .2 .1],'string',[NAME ', ' allUnits{cSpk} ': #' num2str(cSpk)],'FontWeight','bold','interpreter','none','LineStyle','none')
        print(fullfile(figurepath,[NAME,'_BasicActivityPlot.ps']),'-dpsc2','-append','-r300')
        clf
        % save SDF,pt, and rasters
        allTrialTypes {cSpk,curEvt} = trlTypeFieldName;
        allSDF{cSpk,curEvt}         = SDF;
        allPeriTime{cSpk,curEvt}    = pt;
        allRaster{cSpk,curEvt}      = raster;
        allHistogram{cSpk,curEvt}   = histogram;
        allEventLabels{cSpk,curEvt} = evtLabel;
        SDF = {};
        Marker  = SDF;
        raster  = SDF;
        pt      = SDF;
    end
    % save to file
    fname = [NAME,'-BasicSDFs-', num2str(cSpk),'.mat'];
    save(fullfile(analysispath,fname),'allTrialTypes','allSDF','allPeriTime','allRaster','allHistogram','allEventLabels')
end