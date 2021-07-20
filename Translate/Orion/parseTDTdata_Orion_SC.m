function [NeuroPhys,ornTrialdata]=parseTDTdata_Orion_SC(TDT_data,orionEvent,ornTrialdata)
% to do
% incorporate sessions with multiple associated orion files
disp('Parsing TDT/Orion Data ...')
%% get header info
NeuroPhys.HeaderInformation=TDT_data.info;
idx = regexp(TDT_data.info.blockname,'-');
matPath = 'D:\Data\Gollum\Cooling';
subjName = TDT_data.info.blockname(1:idx(1)-1);
if length(idx)>1
    projName = TDT_data.info.blockname(idx(1)+1:idx(2)-1);
    blockID  = TDT_data.info.blockname(idx(2)+1:end);
    matName  = [subjName,'-',projName,'-',blockID];
else
    projName = 'SelfControl';%data.info.blockname(idx(1)+1:idx(2)-1);
    blockID  = TDT_data.info.blockname;
    matName  = [subjName,'-',projName,'-',blockID];
end
NeuroPhys.HeaderInformation.LFPfs = TDT_data.streams.LFPs.fs;
% how many channels 
unit = unique([TDT_data.snips.eSpk.chan,TDT_data.snips.eSpk.sortcode],'rows');
unitlabels='iabcdefghjklmnopqrstuvwxyz';
% load(fullfile(matPath,[matName,'.mat']))
%% parse digital data 
strobevalues_ = TDT_data.epocs.Evnt.data;  % strobe values
evt_ts_ = TDT_data.epocs.Evnt.onset; % strobe ts
% find task ID
% if ~any(strobevalues_ == 201)
%     msg = ['Error: ''', filename, ''', is not a Self-Control task file.'];
%     error(msg)
% end
% 
% // Initialize the Recorder markers database
% // To modify the event code list, see LibCommon/EventCode.cpp
fixation_acquired =  1;
fixation_broken   =  2;
reward            =  9;
fix_cue_on        = 10;
targets_on        = 20;
target_chosen     = 30;
target2_on        = 40;

% find the beginning and end of each trial and end of the last trial
tsIdx   = find(strobevalues_ >= 16385 & strobevalues_ <= 20000); 
tendIdx = find(strobevalues_ >= 30000 & strobevalues_ <= 32766);
% check for multiple 1st trials 

startTrialOrion = find(strobevalues_ >= 16385, 1); 
% if more than one 1st trials find which trial matches up with the orion
% file
if length(startTrialOrion)>1
    for curStrt=1:length(startTrialOrion)
        if curStrt < length(startTrialOrion)
            howManyTrls = sum(evt_ts_ >= evt_ts_(startTrialOrion(curStrt))...
                & evt_ts_ < evt_ts_(startTrialOrion(curStrt+1))...
                & strobevalues_ >= 16385 ...
                & strobevalues_ >= 20000);
        else
            howManyTrls = sum(evt_ts_ >= evt_ts_(startTrialOrion(curStrt))...
                & evt_ts_ < evt_ts_(end)...
                & strobevalues_ >= 16385 ...
                & strobevalues_ >= 20000);
        end
        
        if abs(howManyTrls - size(ornTrialdata.data,1)) <= 1
            firstTrl = startTrialOrion(curStrt);
        end
    end    
else
    firstTrl = startTrialOrion;
end

firstTrlIdx = find(tsIdx==firstTrl);
if length(tsIdx)== size(ornTrialdata.data,1)
    lastTrlIdx = length(tsIdx);
elseif size(ornTrialdata.data,1) > length(tsIdx)
    lastTrlIdx = length(tsIdx);
    ornTrialdata.data = ornTrialdata.data(1:length(tsIdx),:);
else
    lastTrlIdx = firstTrlIdx + size(ornTrialdata.data,1)-1;
end
ntrials = length(firstTrlIdx:lastTrlIdx);
disp([num2str(ntrials), ' trials'])
disp('Behavioral events timestamps ...')
if ~isfield(orionEvent,'trialevents')
    orionEvent.trialevents = cell2mat(orionEvent.mark);
end
for curTrl=firstTrlIdx:lastTrlIdx
    if curTrl==lastTrlIdx % if last trial
        % find events that happen between trial start and trial end events
        TDTCodes = strobevalues_(evt_ts_ >= evt_ts_(tsIdx(curTrl)) ...
            & evt_ts_ < evt_ts_(tendIdx(curTrl)));
        
        TDTCodeTimes = evt_ts_(evt_ts_ >= evt_ts_(tsIdx(curTrl)) ...
            & evt_ts_ < evt_ts_(tendIdx(curTrl)));
    else
        % find events that happen between consecutive trial start events
        TDTCodes = strobevalues_(evt_ts_ >= evt_ts_(tsIdx(curTrl)) ...
            & evt_ts_ < evt_ts_(tsIdx(curTrl+1)));
        
        TDTCodeTimes = evt_ts_(evt_ts_ >= evt_ts_(tsIdx(curTrl)) ...
            & evt_ts_ < evt_ts_(tsIdx(curTrl+1)));
    end
    cdNumIdx = orionEvent.trialevents(:,1) == strobevalues_(tsIdx(curTrl))-16384;
    if isempty(cdNumIdx)
        keyboard
%         error('Missing Trial Event Data')
    end
    OrnCodeNumbers=orionEvent.trialevents(cdNumIdx,3);
    OrnCodeTimes = orionEvent.trialevents(cdNumIdx,4);
    
    [~,IA,IB] = intersect(OrnCodeNumbers,TDTCodes);
    tdt_trl_Ts = ceil((TDTCodeTimes(IB) - TDTCodeTimes(1)).*1000); 
    % make sure the relative time between events in TDT and Orion is not more than 3 ms 
    if ~any(abs(tdt_trl_Ts-OrnCodeTimes(IA))>3)
        NeuroPhys.Digital.CodeNumbers{1,strobevalues_(tsIdx(curTrl))-16384}= TDTCodes;
        NeuroPhys.Digital.CodeTimes{1,strobevalues_(tsIdx(curTrl))-16384} ...
            = ceil((TDTCodeTimes - TDTCodeTimes(1)).*1000);
    else
        keyboard
        error('TDT-Orion Timing Difference > 3 ms')
    end
    %% get spike timestamps
    
end
%% get spike timestamps
disp('Single unit timestamps ...')
for curUnit=1:length(unit)   
    %% cut up trials from the vector and place them in the struct
    chunitname = ['Sig',num2str(unit(curUnit,1)),unitlabels(unit(curUnit,2)+1)];
    disp(chunitname)
    curSpikes = TDT_data.snips.eSpk.ts(TDT_data.snips.eSpk.chan==unit(curUnit,1) ...
        & TDT_data.snips.eSpk.sortcode == unit(curUnit,2));
    curWFs    = TDT_data.snips.eSpk.data(TDT_data.snips.eSpk.chan==unit(curUnit,1) ...
        & TDT_data.snips.eSpk.sortcode == unit(curUnit,2),:);
    eval(['NeuroPhys.Digital.SpikeTimes.',chunitname,'=cell(1,ntrials);;'])
    eval(['NeuroPhys.Digital.SpikeWaveForms.',chunitname,'=cell(1,ntrials);;'])
        
    for curTrl=firstTrlIdx:lastTrlIdx
        if curTrl==lastTrlIdx
            trlSpksIdx =find(curSpikes >= evt_ts_(tsIdx(curTrl)) ...
                & curSpikes < evt_ts_(tendIdx(curTrl)));
        else
            trlSpksIdx =find(curSpikes >= evt_ts_(tsIdx(curTrl)) ...
                & curSpikes < evt_ts_(tsIdx(curTrl+1)));
        end
        trlSpks = (curSpikes(trlSpksIdx)- evt_ts_(tsIdx(curTrl))).*1000;
        trlWfs  = curWFs(trlSpksIdx,:);
        eval(['NeuroPhys.Digital.SpikeTimes.',chunitname,'{1,strobevalues_(tsIdx(curTrl))-16384}=trlSpks;'])
        eval(['NeuroPhys.Digital.SpikeWaveForms.',chunitname,'{1,strobevalues_(tsIdx(curTrl))-16384}=trlWfs;'])
    end
end

disp('LFPs ...')
% create timestamps for lfps
TS = (1: size(TDT_data.streams.LFPs.data,2)) ...
    / TDT_data.streams.LFPs.fs;
NeuroPhys.Analog.LFP.ts=cell(1,ntrials);
for curLFP=1:size(TDT_data.streams.LFPs.data,1)
    if curLFP<10
        chname = ['AD0',num2str(curLFP)];
    else
        chname = ['AD',num2str(curLFP)];
    end
    disp(chname)
    % preallocate trlLFPs
    eval(['NeuroPhys.Analog.LFP.',chname,' = cell(1,ntrials);'])
    for curTrl=firstTrlIdx:lastTrlIdx        
        if curTrl==lastTrlIdx
            trlLFPIdx =find(TS >= evt_ts_(tsIdx(curTrl)) ...
                & TS < evt_ts_(tendIdx(curTrl)));
        else
            trlLFPIdx =find(TS >= evt_ts_(tsIdx(curTrl)) ...
                & TS < evt_ts_(tsIdx(curTrl+1)));
        end
        
        trlts = (TS(trlLFPIdx)- evt_ts_(tsIdx(curTrl))).*1000;
        trlLFPs  = TDT_data.streams.LFPs.data(curLFP,trlLFPIdx);
        if curLFP==1
            NeuroPhys.Analog.LFP.ts{1,strobevalues_(tsIdx(curTrl))-16384}=trlts;
        end
        eval(['NeuroPhys.Analog.LFP.',chname,'{1,strobevalues_(tsIdx(curTrl))-16384}=trlLFPs;'])
        if curLFP==1
            ts = (1: length(trlLFPs)) / TDT_data.streams.LFPs.fs;
            NeuroPhys.Analog.LFP.ts{1,strobevalues_(tsIdx(curTrl))-16384} = ts;
        end
    end
end

end