%% The Application of a Stepwise Generalized Linear Regression Model to Identifiy Self-Control Neuronal Signals
clear all
close all
clc
disp(mfilename)
warning off
% column index of trial information
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
fillclr =[0.5 0.5 0.5
    0.75 0.75 0.75
    0.5 0.5 0.5
    0.75 0.75 0.75];
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


lcs = hot;
lcs = lcs(linspace(1,size(lcs,1),4),:);
untIDs = 'abcdefghjklmnopqrstuvwxyz';
rocBinSz = 30;
rocPt = -500:rocBinSz:500-rocBinSz;
peritime = [-500 500];
lc ='brgk';
figurePath ='D:\SelfControl\ModelFittingNorm';
% figurePath ='T:\VStupho1\SelfControlModelFitting';
errct=1;
winSize = 50;
stepSz =10;
ME={};
%%
datapath = 'D:\SelfControl\data\dotMat';
filelist = dir(fullfile(datapath,'*mat'));
xlims = [-200 300
    -200 300 %-150 200
    0 300];
evtNms = {'target';'saccade';'reward'};
RedoFits = false;
for ii = 1:length(filelist)
    load(fullfile(datapath,filelist(ii).name))
    disp(filelist(ii).name)
    if ~exist('SaccadeData')
        SaccadeData = SaccadeDetectionOrionOO(Behavior,0);
        SaccadeData = getSaccadeOfInterestSC(Behavior,SaccadeData,DerivedData,0);
        save(fullfile(datapath,filelist(ii).name),'SaccadeData','-append')
    end
    if isempty(DerivedData.targets_on)
        DD = parseOrionBehavioralEvents(Behavior);
        DerivedData.targets_on = DD.targets_on;
        save(fullfile(datapath,filelist(ii).name),'DerivedData','-append')
    end
    
    %% sort by target
    targets = unique(SaccadeData.SaccadeDir);
    targets = targets(~isnan(targets));
    
    spkXevt = [1,2,size(SpikeInfo.TrialsPresent,1)];
    AllSigBinsCh1   = zeros(spkXevt); AllSigModelsCh1 = cell(spkXevt); AllSigCoeffsCh1 = cell(spkXevt);
    AllSigBinsCh2   = zeros(spkXevt); AllSigModelsCh2 = cell(spkXevt); AllSigCoeffsCh2 = cell(spkXevt);
    AllSigMdls1     = cell(spkXevt);  AllSigMdls2     = cell(spkXevt);
    AllSpkCnt_vs_SV1 = cell(spkXevt); AllSpkCnt_vs_SV2 = cell(spkXevt);
    sigAUC_Ch1 = cell([spkXevt(3) spkXevt(2)]);
    sigAUC_Ch2 = cell([spkXevt(3) spkXevt(2)]);
    DeltaSVperSpkCh1 = cell(spkXevt); DeltaSVperSpkCh2 = cell(spkXevt);
    AllSpkCnt_vs_SV1 = cell(spkXevt); AllDeltaSV1 = cell(spkXevt); 
    AllSpkCnt_vs_pll_1 = cell(spkXevt); AllSpkCnt_vs_pll_2 = cell(spkXevt);
    AllSpkCnt_vs_SV1Norm = cell(spkXevt); AllDeltaSV1Norm = cell(spkXevt);
    AllSpkCnt_vs_SV2 = cell(spkXevt);     AllDeltaSV2 = cell(spkXevt);
    AllSpkCnt_vs_SV2Norm = cell(spkXevt); AllDeltaSV2Norm = cell(spkXevt);
    allROCaCh1 = cell([spkXevt(3) spkXevt(2)]);
    allROCaCh2 = cell([spkXevt(3) spkXevt(2)]);
    
    for jj = 1:size( SpikeInfo.TrialsPresent,1)
        disp([num2str(jj),': Sig',...
            num2str(SpikeInfo.SingleUnits(jj,1)),untIDs(SpikeInfo.SingleUnits(jj,2))])
        curTrls = SpikeInfo.TrialsPresent(jj,:);
            eval(['curSpk = objNeuroPhys.Digital.SpikeTimes.Sig',...
                num2str(SpikeInfo.SingleUnits(jj,1)),untIDs(SpikeInfo.SingleUnits(jj,2)),';'])
            % convert cell to a matrix
            spkMat= [];
            nspks = zeros(length(curSpk),1);
            for ct = 1:length(curSpk)
                if ~isempty(curSpk{ct})
                    spkMat(ct,1:length(curSpk{ct}))=curSpk{ct};
                    nspks(ct)=length(curSpk{ct});
                end
            end
            for ct = 1:length(curSpk)
                if nspks(ct)~=size(spkMat,2)
                    spkMat(ct,nspks(ct)+1:end)=nan;
                end
            end
            ModelChoice1 ={};
            ModelChoice2 ={};
            spikeCtCoeffModel1=[];
            spikeCtCoeffModel2=[];
                        
            % find trials where spike present
            pltIdx = reshape(1:10,2,5)';%[[1 3; 2 4],[5 7;6 8],[9 11;10 12]];
%             pltIdx(:,2:3)=[];
            clf
            % get choice functions
            [chFxNoTempt, chFxTempt, TemptChFit,noTemptChFit,bNoTempt, ...
                devNoTempt, bTempt, devTempt,statsTempt, statsNoTempt, chFxMdl] ...
                = getChFxsSpikePresent(DerivedData, ...
                SpikeInfo.TrialsPresent(jj,1),SpikeInfo.TrialsPresent(jj,2));            
            %%
            % find the trials when the neuron was present
            ch1_trls = logical(zeros(size(Behavior.TrialInfo.data,1),1));
            
            % get no temptation trials
            ch1_trls(SpikeInfo.TrialsPresent(jj,1):SpikeInfo.TrialsPresent(jj,2))=true;
            ch1_trls = ch1_trls ...
                & Behavior.TrialInfo.data(:,reward) ~= 0 ... % use only correct trials
                & (Behavior.TrialInfo.data(:,cond) == 1)|... % is not a temptation trial;
                Behavior.TrialInfo.data(:,cond) == 3 ;    % is a temptation trial;
            % trial parameters
            ch1_sd          = SaccadeData.SaccadeDir(ch1_trls);
            ch1_rwrdMag     = Behavior.TrialInfo.data(ch1_trls,reward);
            ch1_dly         = Behavior.TrialInfo.data(ch1_trls,ll_delay);
            ch1_choice      = Behavior.TrialInfo.data(ch1_trls,choice1); % 0 = ss, 1 == ll
            % have to invert the mapping to 1 = ss & 0 = ll
            temp = false(size(ch1_choice)); temp(ch1_choice==0)=true; ch1_choice=temp;
            
            % get temptation trials
            ch2_Trls = logical(zeros(size(Behavior.TrialInfo.data,1),1));
            % get trials when neuron present
            ch2_Trls(SpikeInfo.TrialsPresent(jj,1):SpikeInfo.TrialsPresent(jj,2))=true;
            ch2_Trls = ch2_Trls ... % neuron present
                & Behavior.TrialInfo.data(:,reward) ~= 0 ... % correct trial
                & Behavior.TrialInfo.data(:,cond)==3 ; % is a temptation trial
            
            tempt_sd        = SaccadeData.SaccadeDir(ch2_Trls);
            tempt_rwrdMag   = Behavior.TrialInfo.data(ch2_Trls,reward);
            tempt_dly       = Behavior.TrialInfo.data(ch2_Trls,ll_delay);
            tempt_ch1       = Behavior.TrialInfo.data(ch2_Trls,choice1); % 0 = ss, 1 == ll
            ch2             = Behavior.TrialInfo.data(ch2_Trls,choice2); % 0 = ss, 1 == ll
            
            % have to invert the mapping to 1 = ss & 0 = ll
            temp  = false(size(tempt_ch1)); temp(tempt_ch1==0)=true; tempt_ch1=temp;
            temp  = false(size(ch2)); temp(ch2==0)=true; ch2=temp;
            % encode as switch and not p(ss)
            chSwitch = tempt_ch1 ~= ch2; 
            
            % find the trials when the neuron was present
            trls = logical(zeros(size(spkMat,1),1));
            trls(SpikeInfo.TrialsPresent(jj,1):SpikeInfo.TrialsPresent(jj,2)) = true;
            %% get ss trials
            trls(nonzeros(DerivedData.NoTemptationTrials(:,:,1)),2)=true;
            ssTrls = find(trls(:,1) & trls(:,2));
            trls(:,2) = false;
            
            %% get ll trials
            trls(nonzeros(DerivedData.NoTemptationTrials(:,:,2)),2)=true;
            llTrls = find(trls(:,1) & trls(:,2));
            trls(:,2) = false;
            
            %% get LS trials
            trls(nonzeros(DerivedData.TemptationTrials(:,:,1)),2)=true;
            LSTrls = find(trls(:,1) & trls(:,2));
            trls(:,2) = false;
            
            %% get LL trials
            trls(nonzeros(DerivedData.TemptationTrials(:,:,2)),2)=true;
            LLTrls = find(trls(:,1) & trls(:,2));
            trls(:,2) = false;
            
            %% get SL trials
            trls(nonzeros(DerivedData.TemptationTrials(:,:,3)),2)=true;
            SLTrls = find(trls(:,1) & trls(:,2));
            trls(:,2) = false;
            
            %% get SS trials
            trls(nonzeros(DerivedData.TemptationTrials(:,:,4)),2)=true;
            SSTrls = find(trls(:,1) & trls(:,2));
            trls(:,2) = false;
            SV = {}; SV2 = {};
            
            for ll = 1:2 
                % set up plot                              
                subplot(5,2,pltIdx(3,ll));
                cla
                ylabel('p(SS)')
                xlabel('LL delay(ms)')
                ylim([-0.1 1.1])
                set(gca,'tickdir','out',...
                    'PlotBoxAspectRatio',[1.61 1 1],...
                    'ytick',0:.5:1)
                
                switch ll
                    case 1 % aligned on target
                        evt = DerivedData.targets_on;
                        intervals = -200:10:300;
                    case 2 % aligned on saccade
                        evt = SaccadeData.SaccadeChoice1;
                        intervals = -200:10:300;%-150:10:200;
                    case 3 % aligned on saccade but post-saccadic
                        evt = SaccadeData.SaccadeChoice1;
                        intervals = 100:10:260;
                end
                psth = spkMat-repmat(evt(1:size(spkMat,1),1),1,size(spkMat,2));
                
                % get activity on no temptation trials
                [raster, histogram, pt, SDF]...
                    = spikeDensityFunctionTrialBased(curSpk,find(ch1_trls),evt,peritime,2);
                
                % get activity on temptation trials
                [raster_ch2, histogram_ch2, pt, SDF_ch2]...
                    = spikeDensityFunctionTrialBased(curSpk,find(ch2_Trls),evt,peritime,2);
                
                % get activity on ss trials
                [~, ssHist, pt, ssSDF]...
                    = spikeDensityFunctionTrialBased(curSpk,ssTrls,evt,peritime,2);
                % get activity on no temptation trials
                [~, llHist, pt, llSDF]...
                    = spikeDensityFunctionTrialBased(curSpk,llTrls,evt,peritime,2);
                % get activity on no temptation trials
                [~, LSHist, pt, LSSDF]...
                    = spikeDensityFunctionTrialBased(curSpk,LSTrls,evt,peritime,2);
                % get activity on no temptation trials
                [~, LLHist, pt, LLSDF]...
                    = spikeDensityFunctionTrialBased(curSpk,LLTrls,evt,peritime,2);
                % get activity on no temptation trials
                [~, SLHist, pt, SLSDF]...
                    = spikeDensityFunctionTrialBased(curSpk,SLTrls,evt,peritime,2);
                % get activity on no temptation trials
                [~, SSHist, pt, SSSDF]...
                    = spikeDensityFunctionTrialBased(curSpk,SSTrls,evt,peritime,2);
                %%
                subplot(5,2,pltIdx(1,ll))
                cla
                hold on
                plot(pt, ssSDF,'color',TrlTypeColors(1,:))
                plot(pt, llSDF,'color',TrlTypeColors(2,:),'linewidth',2)
                plot(pt, LSSDF,'color',TrlTypeColors(3,:),'linewidth',2)
                plot(pt, LLSDF,'color',TrlTypeColors(4,:),'linewidth',2)
                plot(pt, SLSDF,'color',TrlTypeColors(5,:))
                plot(pt, SSSDF,'color',TrlTypeColors(6,:))
                xlim(xlims(ll,:))
                ylabel('Firing rate (Hz)')
                
                %%           
                spkcts = [];
                tempt_spkcts = [];
                sigCt=1;
                sigBins=[];
                sigSpkCt={};
                sigCoeffs={};
                sigMdls1 = {}; 
                
                sigCt2=1;
                sigBins2=[];
                sigSpkCt2={};
                sigCoeffs2={};
                sigMdls2 = {};
                
                % Stepwise GLM Sliding window analysis
                ssSpkCts =[];
                llSpkCts =[];
                LSSpkCts =[];
                LLSpkCts =[];
                SLSpkCts =[];
                SSSpkCts =[];
                tb =[];
                
                for mm = 1:length(intervals)-(winSize/stepSz)-1
                    % spike counts on no tempation trials in each interval
                    spkcts(:,mm) = sum(psth(ch1_trls,:) >= intervals(mm) ...
                        & psth(ch1_trls,:) < intervals(mm) + winSize,2);                   
                    ssSpkCts(:,mm) = sum(psth(ssTrls,:)  >= intervals(mm) ...
                        & psth(ssTrls,:) < intervals(mm) + winSize,2);
                    llSpkCts(:,mm) = sum(psth(llTrls,:)  >= intervals(mm) ...
                        & psth(llTrls,:) < intervals(mm) + winSize,2);
                    LSSpkCts(:,mm) = sum(psth(LSTrls,:)  >= intervals(mm) ...
                        & psth(LSTrls,:) < intervals(mm) + winSize,2);
                    LLSpkCts(:,mm) = sum(psth(LLTrls,:)  >= intervals(mm) ...
                        & psth(LLTrls,:) < intervals(mm) + winSize,2);
                    SLSpkCts(:,mm) = sum(psth(SLTrls,:)  >= intervals(mm) ...
                        & psth(SLTrls,:) < intervals(mm) + winSize,2);
                    SSSpkCts(:,mm) = sum(psth(SSTrls,:)  >= intervals(mm) ...
                        & psth(SSTrls,:) < intervals(mm) + winSize,2);
                    tb(1,mm) = intervals(mm);% + winSize/2;
                    
                    % spike counts on tempation trials in each interval
                    tempt_spkcts(:,mm) = sum(psth(ch2_Trls,:) >= intervals(mm) ...
                        & psth(ch2_Trls,:) < intervals(mm) + winSize,2);                                        
                end
                
                %% normalize 
                spkctsNorm   = (spkcts  - min(spkcts(:)))/range(spkcts(:));
                normParams   = [min(spkcts(:)) range(spkcts(:))];
                tempt_spkctsNorm =(tempt_spkcts-min(tempt_spkcts(:)))/range(tempt_spkcts(:));
                % choice probability for first choice
                [ROCaCh1, BOOTroca, sigROC_Ch1,ROCci_ch1 ] = ROC_bootstrap([ssSpkCts;SLSpkCts;SSSpkCts], ...
                    [llSpkCts;LSSpkCts;LLSpkCts], 1000, 1, 1);
                
                % choice probability for second choice
                [ROCaCh2, BOOTroca, sigROC_Ch2,ROCci_ch2] = ROC_bootstrap(LSSpkCts, ...
                    LLSpkCts, 1000, 1, 1);
                    
                subplot(5,2,2+ll)
                cla
                hold on
                plot(tb,ROCaCh1,'k')
                plot(tb,ROCci_ch1,'k-.')
                plot(tb,ROCaCh2,'k:')
                plot(tb,ROCci_ch2,'k-.')
                plot(tb(sigROC_Ch1),ROCaCh1(sigROC_Ch1),'or','markerfacecolor','r')
                plot(tb(sigROC_Ch2),ROCaCh2(sigROC_Ch2),'og','markerfacecolor','g')
                xlim(xlims(ll,:))
                ylim([0 1])
                plot(xlim,[0.5 0.5],'k')
                ylabel('AUC')
                
                sigAUC_Ch1{jj,ll} = tb(sigROC_Ch1)';
                sigAUC_Ch2{jj,ll} = tb(sigROC_Ch2)';    
                allROCaCh1{jj,ll} = ROCaCh1;
                allROCaCh2{jj,ll} = ROCaCh2;
                tb =[];
                for mm = 1:length(intervals)-(winSize/stepSz)-1                      
                    % full model for initial choice only uses no temptation
                    modelspec = 'pSS ~ direction*ll_delay*spikeCount - direction:ll_delay:spikeCount';
                    vnames = {'direction' 'll_delay' 'spikeCount' 'pSS' };
                    
                    % Xds = mat2dataset([ch1_sd ch1_rwrdMag ch1_dly spkcts(:,mm) ch1_choice],'VarNames', vnames);
                    Xds = mat2dataset([ch1_sd ch1_dly spkcts(:,mm) ch1_choice],'VarNames', vnames);
                    ModelChoice1{ll,mm} = stepwiseglm(Xds,'constant', ...
                        'upper', 'interactions',...
                        'lower', 'pSS ~ ll_delay', ...
                        'Distribution','binomial', ...
                        'ResponseVar','pSS', ...
                        'CategoricalVars',1, ...
                        'Verbose',0); 
                    tb(1,mm) = intervals(mm);
                    
                    XdsNorm = mat2dataset([ch1_sd ch1_dly spkctsNorm(:,mm) ch1_choice],'VarNames', vnames);
                    ModelChoice1Norm{ll,mm} = stepwiseglm(XdsNorm,'constant', ...
                        'upper', 'interactions',...
                        'lower', 'pSS ~ ll_delay', ...
                        'Distribution','binomial', ...
                        'ResponseVar','pSS', ...
                        'CategoricalVars',1, ...
                        'Verbose',0);
                    
                    %%
                    % test if any of the coefficients are not significant
                    if any(contains(ModelChoice1{ll,mm}.CoefficientNames,'spikeCount'))
                        [ModelChoice1{ll,mm}] = sigCoeffsTest(ModelChoice1{ll,mm});
                    end
                    % test if any of the coefficients are not significant
                    if any(contains(ModelChoice1Norm{ll,mm}.CoefficientNames,'spikeCount'))
                        [ModelChoice1Norm{ll,mm}] = sigCoeffsTest(ModelChoice1Norm{ll,mm});
                    end
                    if any(contains(ModelChoice1Norm{ll,mm}.CoefficientNames,'spikeCount'))
                        %%
                        evtTm = [num2str(intervals(mm)) ' ms on ' evtNms{ll}];
                        disp([evtTm ', p(ss1) = ' ModelChoice1{ll,mm}.Formula.LinearPredictor])
                        disp([evtTm ', p(ss1)Norm = ' ModelChoice1Norm{ll,mm}.Formula.LinearPredictor])
                        figure
                        set(gcf,'units','pixels',...
                            'position',[-1215 9 608 988])
                        
                        % ch functions all trials
                        clf
                        subplot(5,2,1)
                        cla
                        xval = min(DerivedData.ChoiceFunctionNoTemptation(:,1)):100:max(DerivedData.ChoiceFunctionNoTemptation(:,1));
                        hold on
                        plot(chFxNoTempt(:,1)/1000,chFxNoTempt(:,2),'ko', ...
                            chFxTempt(:,1,1)/1000,chFxTempt(:,2,1),'bo', ...                            
                            xval/1000, noTemptChFit,'k',...
                            xval/1000, TemptChFit(:,1),'b',...
                            chFxTempt(:,1,2)/1000,chFxTempt(:,2,2),'ro', ...
                            xval/1000, TemptChFit(:,2),'r')
                        xlabel('Delay(ms)')
                        ylabel('p(ss)')
                        title('Choice functions: Initial Choice')
                        set(gca,'tickdir','out')
                        
                        % sdfs
                        subplot(5,2,3)
                        cla
                        hold on
                        plot(pt, ssSDF,'color',TrlTypeColors(1,:))
                        plot(pt, llSDF,'color',TrlTypeColors(2,:),'linewidth',2)
                        plot(pt, LSSDF,'color',TrlTypeColors(3,:),'linewidth',2)
                        plot(pt, LLSDF,'color',TrlTypeColors(4,:),'linewidth',2)
                        plot(pt, SLSDF,'color',TrlTypeColors(5,:))
                        plot(pt, SSSDF,'color',TrlTypeColors(6,:))
                        ys = sort([ylim ylim]);
                        xs = [intervals(mm) intervals(mm) + winSize ...
                            intervals(mm) + winSize intervals(mm)];
                        fill(xs,ys, [1 1 1]*.5)
                        alpha(0.25)
                        xlim(xlims(ll,:))
                        ylabel('Firing rate (Hz)')
                        if ll == 1
                            xlabel('Time from target(ms)')
                        else
                            xlabel('Time from saccade(ms)')
                        end
                        title([num2str(intervals(mm)),' : ', num2str(intervals(mm) + winSize),' ms'])
                        set(gca,'tickdir','out')
                        
                        subplot(5,2,5)
                        cla
                        hold on
                        if ~isempty(SSSpkCts)
                            plot(1, ssSpkCts(:,mm),'.','color',TrlTypeColors(1,:))
                            plot(1,mean(ssSpkCts(:,mm)),'o','color',TrlTypeColors(1,:))
                        end
                        plot(2, llSpkCts(:,mm),'.','color',TrlTypeColors(2,:))
                        plot(2,mean(llSpkCts(:,mm)),'o','color',TrlTypeColors(2,:))
                        if ~isempty(LSSpkCts(:,mm))
                            plot(3, LSSpkCts(:,mm),'.','color',TrlTypeColors(3,:))
                            plot(3,mean(LSSpkCts(:,mm)),'o','color',TrlTypeColors(3,:))
                        end
                        if ~isempty(LLSpkCts(:,mm))
                            plot(4, LLSpkCts(:,mm),'.','color',TrlTypeColors(4,:))
                            plot(4,mean(LLSpkCts(:,mm)),'o','color',TrlTypeColors(4,:))
                        end
                        if ~isempty(SLSpkCts)
                            plot(5, SLSpkCts(:,mm),'.','color',TrlTypeColors(5,:))
                            plot(5,mean(SLSpkCts(:,mm)),'o','color',TrlTypeColors(5,:))
                        end
                        if ~isempty(SSSpkCts)
                            plot(6, SSSpkCts(:,mm),'.','color',TrlTypeColors(6,:))
                            plot(6,mean(SSSpkCts(:,mm)),'o','color',TrlTypeColors(6,:))
                        end
                        xlim([0 7])
                        ylabel('Spike count')
                        xlabel('Trial type')
                        title('Raw Spike Counts')
                        set(gca,'xtick',1:6,...
                            'xticklabel',{'ss','ll','LS','LL','SL','SS'},...
                            'tickdir','out')
                        
                        subplot(5,2,6)
                        cla
                        hold on
                        if ~isempty(SSSpkCts)
                            plot(1, (ssSpkCts(:,mm)-normParams(1))/normParams(2),'.','color',TrlTypeColors(1,:))
                            plot(1,(mean(ssSpkCts(:,mm))-normParams(1))/normParams(2),'o','color',TrlTypeColors(1,:))
                        end
                        plot(2, (llSpkCts(:,mm)-normParams(1))/normParams(2),'.','color',TrlTypeColors(2,:))
                        plot(2,(mean(llSpkCts(:,mm))-normParams(1))/normParams(2),'o','color',TrlTypeColors(2,:))
                        if ~isempty(LSSpkCts)
                            plot(3, (LSSpkCts(:,mm)-normParams(1))/normParams(2),'.','color',TrlTypeColors(3,:))
                            plot(3,(mean(LSSpkCts(:,mm))-normParams(1))/normParams(2),'o','color',TrlTypeColors(3,:))
                        end
                        if ~isempty(LLSpkCts(:,mm))
                            plot(4, (LLSpkCts(:,mm)-normParams(1))/normParams(2),'.','color',TrlTypeColors(4,:))
                            plot(4,(mean(LLSpkCts(:,mm))-normParams(1))/normParams(2),'o','color',TrlTypeColors(4,:))
                        end
                        if ~isempty(SLSpkCts)
                            plot(5, (SLSpkCts(:,mm)-normParams(1))/normParams(2),'.','color',TrlTypeColors(5,:))
                            plot(5,(mean(SLSpkCts(:,mm))-normParams(1))/normParams(2),'o','color',TrlTypeColors(5,:))
                        end
                        if ~isempty(SSSpkCts)
                            plot(6, (SSSpkCts(:,mm)-normParams(1))/normParams(2),'.','color',TrlTypeColors(6,:))
                            plot(6,(mean(SSSpkCts(:,mm))-normParams(1))/normParams(2),'o','color',TrlTypeColors(6,:))
                        end
                        
                        xlim([0 7])
                        ylabel('Spike count')
                        xlabel('Trial type')
                        title('Normalized Spike Counts')
                        set(gca,'xtick',1:6,...
                            'xticklabel',{'ss','ll','LS','LL','SL','SS'},...
                            'tickdir','out')
                        
                        annotation('textbox',...
                            'string',[filelist(ii).name(1:end-4) ', Sig',num2str(SpikeInfo.SingleUnits(jj,1)),untIDs(SpikeInfo.SingleUnits(jj,2)),': Initial Choice'] ,...
                            'position',[0.1 0.9 0.8 0.1],...
                            'linestyle','none',...
                            'fontsize',14)
                        annotation('textbox',...
                            'string',{evtTm; ['p(ss1) = ' ModelChoice1{ll,mm}.Formula.LinearPredictor]},...
                            'position',[0.05 0.3 0.8 0.1],...
                            'linestyle','none',...
                            'fontsize',12,...
                            'interpreter','none')
                        
                        % get the table with the coefficients
                        T = ModelChoice1Norm{ll,mm}.Coefficients;
                        % Get the table in string form.
                        TString = evalc('disp(T)');
                        % Use TeX Markup for bold formatting and underscores.
                        TString = strrep(TString,'<strong>','\bf');
                        TString = strrep(TString,'</strong>','\rm');
                        TString = strrep(TString,'_','\_');
                        % Get a fixed-width font.
                        FixedWidth = get(0,'FixedWidthFontName');
                        % Output the table using the annotation command.
                        annotation(gcf,'Textbox',...
                            'String',TString,...
                            'Interpreter','Tex',...
                            'FontName',FixedWidth,...
                            'linestyle','none',...
                            'Units','Normalized',...
                            'Position',[0.05,0.05,0.9,0.30],...
                            'FontSize',10);
                        
                        % get delta sv
                        [fits] = getDeltaSV_Norm(ModelChoice1Norm{ll,mm});
                        AllSpkCnt_vs_pll_1{sigCt,ll,jj} = fits.AreaUnderChFx;
                        
                        subplot(5,2,2)
                        cla
                        hold on
                        pll_Area = [];
                        
                        plot(xval/1000, 1-TemptChFit(:,1),'b','linewidth',2)
                        plot(xval/1000, 1-TemptChFit(:,2),'r','linewidth',2)
                        for cf = 1:size(fits.yFits,3)
                            plot(fits.delays/1000,1-fits.yFits(:,1,cf))
                            plot(fits.delays/1000,1-fits.yFits(:,2,cf))
                        end                                                
                        %%
                        lgndNums=[];
                        [m n o] = size(fits.predictorValues);
                        tot = cumprod([m n o]); tot = tot(end);
                        cols = ModelChoice1Norm{ll, mm}.NumCoefficients-1;
%                         clc
                        ct = 1; ct1 = 0;
                        for cf = 1:tot
%                             fits.predictorValues(cf)                            
                            if ct1 < cols
                                ct1 = ct1+1;
                            else
                                ct1=1;
                                ct = ct+1;
                            end
%                             if ct1 == cols
%                                 ct = ct+1;
%                             end
                            lgndNums(ct, ct1) = fits.predictorValues(cf);
                        end
                        %%
                        lgnd = ['I';'F'];
                        tmp  = num2str(lgndNums);
                        lgnd(3:3+size(tmp,1)-1,1:size(tmp,2)) = tmp;
                        legend(lgnd,'location','best')
                        legend('boxoff')
                        xlabel('Delay(ms)')
                        ylabel('p(ll)')
                        title({'Choice functions'; 'vs norm Spike count'})
                        set(gca,'tickdir','out')
                        
                        subplot(5,2,4)
                        cla
                        hold on
                        bar(fits.AreaUnderChFx(:))
                        ylabel('Area under p(ll)')
                        xlabel('Normalized Spike Count')
                        title('Delta p(ll)')
                        text(1:length(fits.AreaUnderChFx(:)),fits.AreaUnderChFx(:),num2str(round(fits.AreaUnderChFx(:),2)), ...
                            'horizontalalignment','center',...
                            'verticalalignment','bottom')
                        set(gca,'tickdir','out',...
                            'xtick',1:length(fits.AreaUnderChFx(:)),...
                            'xticklabels', num2str(lgndNums))
                        xtickangle(-30)
                        
                        % set up to print to file
                        set(gcf,'units','inches'...
                            ,'paperposition',[0.25 0.25 8   10.5])
                        docname = [filelist(ii).name(1:end-4) '-Sig',num2str(SpikeInfo.SingleUnits(jj,1)),untIDs(SpikeInfo.SingleUnits(jj,2)),'_Initial_Choice_', evtTm,'.pdf'];
                        % print to file
                        print(fullfile(figurePath,docname),'-dpdf')
                        % delete figure
                        delete(gcf)
                        %%                                                
                        spkIdx = find(contains(ModelChoice1{ll,mm}.CoefficientNames,'spikeCount'));
                        sigBins(sigCt,1)   = intervals(mm);
                        sigSpkCt{sigCt,1}  = ModelChoice1{ll,mm};
                        
                        if any(contains(ModelChoice1{ll,mm}.CoefficientNames,'spikeCount'))
                            sigCoeffs{sigCt,1} = ModelChoice1{ll,mm}.CoefficientNames{contains(ModelChoice1{ll,mm}.CoefficientNames,'spikeCount')};
                        else
                            sigCoeffs{sigCt,1} = '';
                        end
                        sigMdls1{sigCt,1}  = ModelChoice1{ll,mm}.Formula.LinearPredictor;
                        
                        sigCoeffsNorm{sigCt,1} = ModelChoice1Norm{ll,mm}.CoefficientNames{contains(ModelChoice1Norm{ll,mm}.CoefficientNames,'spikeCount')};
                        sigMdls1Norm{sigCt,1} = ModelChoice1Norm{ll,mm}.Formula.LinearPredictor;
                        
                        
                        AllSigBinsCh1(sigCt,ll,jj)    = sigBins(sigCt,1);
                        AllSigModelsCh1{sigCt,ll,jj}  = sigSpkCt{sigCt,1};
                        AllSigCoeffsCh1{sigCt,ll,jj}  = sigCoeffs{sigCt,1};
                        AllSigMdls1{sigCt,ll,jj}      = sigMdls1{sigCt,1};
                        [AllSpkCnt_vs_SV1{sigCt,ll,jj}, AllDeltaSV1{sigCt,ll,jj},Fits1] ...
                            = getDeltaSV(ModelChoice1{ll,mm});
                        
                        fits = getDeltaSV_Norm(ModelChoice1Norm{ll,mm});
                        AllSpkCnt_vs_pll_1{sigCt,ll,jj} = fits.AreaUnderChFx;
                        
                        % plot SV shifts                          
                        for pp = 1:size(fits.yFits,3)
                            subplot(5,2,4+ll)
                            hold on
                            switch pp
                                case 1
                                    cs = parula;
                                case 2
                                    cs = jet;
                                case 3
                                    cs = copper;
                                case 4
                                    cs = hsv;
                                case 5
                                    cs = cool;
                                case 6
                                    cs = spring;
                                case 7
                                    cs = summer;
                                case 8
                                    cs = autumn;
                                case 9
                                    cs = winter;
                            end
                            css = cs(floor(linspace(1,64,size(fits.yFits,2))),:);
                            for oo = 1:size(fits.yFits,2)
                                plot(fits.delays/1000, fits.yFits(:,oo,pp),'color',css(oo,:))
                            end
                            plot(AllSpkCnt_vs_SV1{sigCt,ll,jj}(:,2,pp)/1000,0.5,'o','color',css(1,:))
                            %%
                            subplot(5,2,8+ll)
                            hold on
                            tmp = diff(fits.AreaUnderChFx);
                            [~, idx] = max(abs(tmp));
                            plot(sigBins(sigCt,1),tmp(idx),'or','markerfacecolor','r')
                            xlim(xlims(ll,:))
                        end                                                
                        %%
                        sigCt = sigCt + 1;
                    end
                    % full model for final choice only uses temptation
                    modelspec = 'pSS2 ~ direction*ll_delay*spikeCount - direction:ll_delay:spikeCount';
                    vnames = {'direction' 'll_delay' 'spikeCount' 'pSS2' };
                    
                    Xds = mat2dataset([tempt_sd tempt_dly tempt_spkcts(:,mm) ch2],'VarNames', vnames);
                    ModelChoice2{ll,mm} = stepwiseglm(Xds,'constant', ...
                        'upper', 'interactions',...
                        'lower', 'pSS2 ~ ll_delay', ...
                        'Distribution','binomial', ...
                        'ResponseVar','pSS2', ...
                        'CategoricalVars',1, ...
                        'Verbose',0);
                    
                    Xds = mat2dataset([tempt_sd tempt_dly tempt_spkctsNorm(:,mm) ch2],'VarNames', vnames);
                    ModelChoice2Norm{ll,mm} = stepwiseglm(Xds,'constant', ...
                        'upper', 'interactions',...
                        'lower', 'pSS2 ~ ll_delay', ...
                        'Distribution','binomial', ...
                        'ResponseVar','pSS2', ...
                        'CategoricalVars',1, ...
                        'Verbose',0);
                    
%                     modelspec = 'pSwitch ~ direction*ll_delay*spikeCount - direction:ll_delay:spikeCount';
%                     vnames = {'direction' 'll_delay' 'spikeCount' 'pSwitch' };
%                     
%                     Xds = mat2dataset([tempt_sd tempt_dly tempt_spkcts(:,mm) chSwitch],'VarNames', vnames);
%                     ModelSwitch{ll,mm} = stepwiseglm(Xds,'constant', ...
%                         'upper', 'interactions',...
%                         'lower', 'pSwitch ~ ll_delay', ...
%                         'Distribution','binomial', ...
%                         'ResponseVar','pSwitch', ...
%                         'CategoricalVars',1, ...
%                         'Verbose',0);
                    
                    % test if any of the coefficients are significant
                    if any(contains(ModelChoice2{ll,mm}.CoefficientNames,'spikeCount'))
                        [ModelChoice2{ll,mm}] = sigCoeffsTest(ModelChoice2{ll,mm});
                    end
                    % test if any of the coefficients are significant
                    if any(contains(ModelChoice2Norm{ll,mm}.CoefficientNames,'spikeCount'))
                        [ModelChoice2Norm{ll,mm}] = sigCoeffsTest(ModelChoice2Norm{ll,mm});
                    end
                    
                    if any(contains(ModelChoice2Norm{ll,mm}.CoefficientNames,'spikeCount'))
                        %%
                        evtTm = [num2str(intervals(mm)) ' ms on ' evtNms{ll}];
                        disp([evtTm ', p(ss2) = ' ModelChoice2{ll,mm}.Formula.LinearPredictor])
                        disp([evtTm ', p(ss2)Norm = ' ModelChoice2Norm{ll,mm}.Formula.LinearPredictor])
                        
                        set(gcf,'units','pixels',...
                            'position',[-1215 9 608 988])
                        % ch functions all trials
                        figure
                        set(gcf,'units','pixels',...
                            'position',[-1268 9 662 988])

                        clf
                        subplot(5,2,1)
                        cla
                        xval = min(DerivedData.ChoiceFunctionNoTemptation(:,1)):100:max(DerivedData.ChoiceFunctionNoTemptation(:,1));
                        hold on
                        plot(chFxTempt(:,1,1)/1000,chFxTempt(:,2,1),'bo', ...                            
                            xval/1000, TemptChFit(:,1),'b',...
                            chFxTempt(:,1,2)/1000,chFxTempt(:,2,2),'ro', ...
                            xval/1000, TemptChFit(:,2),'r')
                        xlabel('Delay(ms)')
                        ylabel('p(ss)')
                        title('Choice functions: Final Choice')
                        set(gca,'tickdir','out')
                        
                        % sdfs
                        subplot(5,2,3)
                        cla
                        hold on
                        plot(pt, LSSDF,'color',TrlTypeColors(3,:),'linewidth',2)
                        plot(pt, LLSDF,'color',TrlTypeColors(4,:),'linewidth',2)
                        plot(pt, SLSDF,'color',TrlTypeColors(5,:))
                        plot(pt, SSSDF,'color',TrlTypeColors(6,:))
                        ys = sort([ylim ylim]);
                        xs = [intervals(mm) intervals(mm) + winSize ...
                            intervals(mm) + winSize intervals(mm)];
                        fill(xs,ys, [1 1 1]*.5)
                        alpha(0.25)
                        xlim(xlims(ll,:))
                        ylabel('Firing rate (Hz)')
                        if ll == 1
                            xlabel('Time from target(ms)')
                        else
                            xlabel('Time from saccade(ms)')
                        end
                        title([num2str(intervals(mm)),' : ', num2str(intervals(mm) + winSize),' ms'])
                        set(gca,'tickdir','out')
                        
                        subplot(5,2,5)
                        cla
                        hold on                        
                        plot(3, LSSpkCts(:,mm),'.','color',TrlTypeColors(3,:))
                        plot(3,mean(LSSpkCts(:,mm)),'o','color',TrlTypeColors(3,:))
                        plot(4, LLSpkCts(:,mm),'.','color',TrlTypeColors(4,:))
                        plot(4,mean(LLSpkCts(:,mm)),'o','color',TrlTypeColors(4,:))
                        if ~isempty(SLSpkCts)
                            plot(5, SLSpkCts(:,mm),'.','color',TrlTypeColors(5,:))
                            plot(5,mean(SLSpkCts(:,mm)),'o','color',TrlTypeColors(5,:))
                        end
                        if ~isempty(SSSpkCts)
                            plot(6, SSSpkCts(:,mm),'.','color',TrlTypeColors(6,:))
                            plot(6,mean(SSSpkCts(:,mm)),'o','color',TrlTypeColors(6,:))
                        end
                        xlim([0 7])
                        ylabel('Spike count')
                        xlabel('Trial type')
                        title('Raw Spike Counts')
                        set(gca,'xtick',1:6,...
                            'xticklabel',{'ss','ll','LS','LL','SL','SS'},...
                            'tickdir','out')
                        
                        subplot(5,2,6)
                        cla
                        hold on                        
                        plot(3, (LSSpkCts(:,mm)-normParams(1))/normParams(2),'.','color',TrlTypeColors(3,:))
                        plot(3,(mean(LSSpkCts(:,mm))-normParams(1))/normParams(2),'o','color',TrlTypeColors(3,:))
                        plot(4, (LLSpkCts(:,mm)-normParams(1))/normParams(2),'.','color',TrlTypeColors(4,:))
                        plot(4,(mean(LLSpkCts(:,mm))-normParams(1))/normParams(2),'o','color',TrlTypeColors(4,:))
                        if ~isempty(SLSpkCts)
                            plot(5, (SLSpkCts(:,mm)-normParams(1))/normParams(2),'.','color',TrlTypeColors(5,:))
                            plot(5,(mean(SLSpkCts(:,mm))-normParams(1))/normParams(2),'o','color',TrlTypeColors(5,:))
                        end
                        if ~isempty(SSSpkCts)
                            plot(6, (SSSpkCts(:,mm)-normParams(1))/normParams(2),'.','color',TrlTypeColors(6,:))
                            plot(6,(mean(SSSpkCts(:,mm))-normParams(1))/normParams(2),'o','color',TrlTypeColors(6,:))
                        end
                        xlim([0 7])
                        ylabel('Spike count')
                        xlabel('Trial type')
                        title('Normalized Spike Counts')
                        set(gca,'xtick',1:6,...
                            'xticklabel',{'ss','ll','LS','LL','SL','SS'},...
                            'tickdir','out')
                        
                        annotation('textbox',...
                            'string',[filelist(ii).name(1:end-4) ', Sig',num2str(SpikeInfo.SingleUnits(jj,1)),untIDs(SpikeInfo.SingleUnits(jj,2)),': Final Choice'] ,...
                            'position',[0.1 0.9 0.8 0.1],...
                            'linestyle','none',...
                            'fontsize',14)
                        
                        annotation('textbox',...
                            'string',{evtTm; ['p(ss2) = ' ModelChoice2Norm{ll,mm}.Formula.LinearPredictor]},...
                            'position',[0.05 0.3 0.8 0.1],...
                            'linestyle','none',...
                            'fontsize',12,...
                            'interpreter','none')
                        
                        % get the table with the coefficients
                        T = ModelChoice2Norm{ll,mm}.Coefficients;
                        % Get the table in string form.
                        TString = evalc('disp(T)');
                        % Use TeX Markup for bold formatting and underscores.
                        TString = strrep(TString,'<strong>','\bf');
                        TString = strrep(TString,'</strong>','\rm');
                        TString = strrep(TString,'_','\_');
                        % Get a fixed-width font.
                        FixedWidth = get(0,'FixedWidthFontName');
                        % Output the table using the annotation command.
                        annotation(gcf,'Textbox',...
                            'String',TString,...
                            'Interpreter','Tex',...
                            'FontName',FixedWidth,...
                            'linestyle','none',...
                            'Units','Normalized',...
                            'Position',[0.05,0.05,0.9,0.30],...
                            'FontSize',10);
                        
                        % get delta sv
                        [fits] = getDeltaSV_Norm(ModelChoice2Norm{ll,mm});
                        AllSpkCnt_vs_pll_2{sigCt,ll,jj} = fits.AreaUnderChFx;
                        
                        subplot(5,2,2)
                        cla
                        hold on
                        pll_Area = [];
                        
                        plot(xval/1000, 1-TemptChFit(:,1),'b','linewidth',2)
                        plot(xval/1000, 1-TemptChFit(:,2),'r','linewidth',2)
                        for cf = 1:size(fits.yFits,3)
                            plot(fits.delays/1000,1-fits.yFits(:,1,cf))
                            plot(fits.delays/1000,1-fits.yFits(:,2,cf))
                        end                                                
                        %%
                        lgndNums=[];
                        [m n o] = size(fits.predictorValues);
                        tot = cumprod([m n o]); tot = tot(end);
                        cols = m;

                        ct = 1; ct1 = 0;
                        for cf = 1:tot                           
                            if ct1 < cols
                                ct1 = ct1+1;
                            else
                                ct1=1;
                                ct = ct+1;
                            end
                            lgndNums(ct, ct1) = fits.predictorValues(cf);
                        end
                        %%
                        lgnd = ['I';'F'];
                        tmp  = num2str(lgndNums);
                        lgnd(3:3+size(tmp,1)-1,1:size(tmp,2)) = tmp;
                        legend(lgnd,'location','best')
                        legend('boxoff')
                        xlabel('Delay(ms)')
                        ylabel('p(ll)')
                        title({'Choice functions'; 'vs norm Spike count'})
                        set(gca,'tickdir','out')
                        
                        subplot(5,2,4)
                        cla
                        hold on
                        bar(fits.AreaUnderChFx(:))
                        ylabel('Area under p(ll)')
                        xlabel('Normalized Spike Count')
                        title('Delta p(ll)')
                        text(1:length(fits.AreaUnderChFx(:)),fits.AreaUnderChFx(:),num2str(round(fits.AreaUnderChFx(:),2)), ...
                            'horizontalalignment','center',...
                            'verticalalignment','bottom')
                        set(gca,'tickdir','out',...
                            'xtick',1:length(fits.AreaUnderChFx(:)),...
                            'xticklabels', num2str(lgndNums))
                        xtickangle(-30)
                        
                        % set up to print to file
                        set(gcf,'units','inches'...
                            ,'paperposition',[0.25 0.25 8   10.5])
                        %%
                        docname = [filelist(ii).name(1:end-4) '-Sig',num2str(SpikeInfo.SingleUnits(jj,1)),untIDs(SpikeInfo.SingleUnits(jj,2)),'_Final_Choice_', evtTm,'.pdf'];
                        % print to file
                        print(fullfile(figurePath,docname),'-dpdf')
                        % delete figure
                        delete(gcf)
                        %%                        
                        
                        spkIdx = find(contains(ModelChoice2{ll,mm}.CoefficientNames,'spikeCount'));
                        sigBins2(sigCt2,1)   = intervals(mm);
                        sigSpkCt2{sigCt2,1}  = ModelChoice2{ll,mm};
                        sigCoeffs2{sigCt2,1} = ModelChoice2{ll,mm}.CoefficientNames{contains(ModelChoice2{ll,mm}.CoefficientNames,'spikeCount')};
                        sigMdls2{sigCt2,1}   = ModelChoice2{ll,mm}.Formula.LinearPredictor;
                        
                        sigCoeffs2Norm{sigCt2,1} = ModelChoice2Norm{ll,mm}.CoefficientNames{contains(ModelChoice2Norm{ll,mm}.CoefficientNames,'spikeCount')};
                        sigMdls2Norm{sigCt2,1}   = ModelChoice2Norm{ll,mm}.Formula.LinearPredictor;
                        
                        
                        AllSigBinsCh2(sigCt2,ll,jj)   = sigBins2(sigCt2,1);
                        AllSigModelsCh2{sigCt2,ll,jj} = sigSpkCt2{sigCt2,1};
                        AllSigCoeffsCh2{sigCt2,ll,jj} = sigCoeffs2{sigCt2,1};
                        AllSigMdls2{sigCt2,ll,jj}     = sigMdls2{sigCt2,1};
                        
                        [AllSpkCnt_vs_SV2{sigCt2,ll,jj}, AllDeltaSV2{sigCt2,ll,jj}, Fits] ...
                            = getDeltaSV(ModelChoice2{ll,mm});
                                                
                        %%
                        
                        for pp = 1:size(fits.yFits,3)
                            subplot(5,2,6+ll)
                            hold on
                            switch pp
                                case 1
                                    cs = parula;
                                case 2
                                    cs = jet;
                                case 3
                                    cs = copper;
                                case 4
                                    cs = hsv;
                                case 5
                                    cs = cool;
                                case 6
                                    cs = spring;
                                case 7
                                    cs = summer;
                                case 8
                                    cs = autumn;
                                case 9
                                    cs = winter;
                            end
                            css = cs(floor(linspace(1,64,size(fits.yFits,2))),:);
                            for oo = 1:size(fits.yFits,2)
                                plot(fits.delays/1000, fits.yFits(:,oo,pp),'color',css(oo,:))
                            end
                            plot(AllSpkCnt_vs_SV2{sigCt2,ll,jj}(:,2,pp)/1000,0.5,'o','color',css(1,:))
                            %%
                            subplot(5,2,8+ll)
                            hold on
                            tmp = diff(fits.AreaUnderChFx);
                            [~, idx] = max(abs(tmp));
                            plot(sigBins2(sigCt2,1),tmp(idx),'og','markerfacecolor','g')
                            xlim(xlims(ll,:))
                        end
                        
                        sigCt2 = sigCt2 + 1;                        
                    end
                end                 
            end
            %%
            y=[];
            subplot(5,2,1)
            y(1) = max(ylim);
            subplot(5,2,2)
            y(2) = max(ylim);
            
            subplot(5,2,1)
            ylim([0 max(y)])
            subplot(5,2,2)
            ylim([0 1]*max(y))
            
            subplot(5,2,3)
            xlabel('Time from target (ms)')
            
            subplot(5,2,4)
            xlabel('Time from 1st saccade (ms)')
            
            subplot(5,2,5)
            title('Initial Choice')
            xlabel('ll delay (sec)')
            ylabel('p(ss)')
            
            subplot(5,2,6)
            title('Initial Choice') 
            xlabel('ll delay (sec)')
            ylabel('p(ss)')
            
            subplot(5,2,7)
            title('Final Choice')
            xlabel('ll delay (sec)')
            ylabel('p(ss)')
            
            subplot(5,2,8)
            title('Final Choice')
            xlabel('ll delay (sec)')
            ylabel('p(ss)')
            
            subplot(5,2,9)
            plot(xlim,[0 0],':k')
            xlabel('Time from target (ms)')
            ylabel('{\Delta}p(ll)')
            y = max(ylim);
            if y ~= 0
                ylim([-1 1]*y)
            else
                ylim([-1 1])
            end
            
            subplot(5,2,10)
            xlabel('Time from 1st saccade (ms)')
            ylabel('{\Delta}SV per spike (sec)')
            plot(xlim,[0 0],':k')
            y = max(ylim);
            if y ~= 0
                ylim([-1 1]*y)
            else
                ylim([-1 1])
            end
            %%
            annotation('textbox',...
                'string',[filelist(ii).name(1:end-4) ': Sig',num2str(SpikeInfo.SingleUnits(jj,1)),untIDs(SpikeInfo.SingleUnits(jj,2))] ,...
                'position',[0.1 0.9 0.8 0.1],...
                'linestyle','none',...
                'fontsize',14)
            %         linkaxes(ax)
            set(gcf,'units','inches'...
                ,'paperposition',[0.25 0.25 8   10.5])
            print(fullfile(figurePath,[filelist(ii).name(1:end-4) ,'_', num2str(jj),'.pdf']),'-dpdf')            
    end
    save(fullfile(datapath,filelist(ii).name),'AllSigBinsCh1','AllSigModelsCh1',...
        'AllSigCoeffsCh1','AllSigBinsCh2','AllSigModelsCh2','AllSigCoeffsCh2',...
        'sigAUC_Ch1','sigAUC_Ch2','allROCaCh1',...
        'allROCaCh2','AllSpkCnt_vs_SV1','AllSpkCnt_vs_SV2',...
        'AllSpkCnt_vs_SV1','AllSpkCnt_vs_SV2','AllDeltaSV1','AllDeltaSV2',...
        'AllDeltaSV1Norm','AllDeltaSV2Norm','AllSpkCnt_vs_pll_1','AllSpkCnt_vs_pll_2','-append')
    
    clear AnatomicalInfo Behavior DerivedData objNeuroPhys SpikeInfo SaccadeData AllSigBinsCh1 ...
        AllSigModelsCh1 AllSigCoeffsCh1 AllSigBinsCh2 AllSigModelsCh2 AllSigCoeffsCh2 ...
        AllSpkCnt_vs_SV1 AllSpkCnt_vs_SV2 AllSigBinsCh1 AllSigModelsCh1 ...
        AllSigCoeffsCh1 AllSigBinsCh2 AllSigModelsCh2 AllSigCoeffsCh2 ...
        sigAUC_Ch1 sigAUC_Ch2 DeltaSVperSpkCh1 DeltaSVperSpkCh2 allROCaCh1 ...
        allROCaCh2 AllSpkCnt_vs_SV1 AllSpkCnt_vs_SV2 ...
        AllSpkCnt_vs_SV1 AllSpkCnt_vs_SV2 AllDeltaSV1 AllDeltaSV2 AllDeltaSV1Norm AllDeltaSV2Norm...
        AllSpkCnt_vs_pll_1 AllSpkCnt_vs_pll_2
    
end