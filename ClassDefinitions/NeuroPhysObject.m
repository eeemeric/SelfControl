classdef (ConstructOnLoad = false) NeuroPhysObject %< NeuroPhysExpObject
    %NeuroPhys NeuroPhys Class definition for the Stuphorn lab standard
    % matlab object for storing neurophysiological data
    %
    % The NPstruct returned by Plx2mat_ML_SC has the following properties
    % and is organized in the same fashion as the digital events in the
    % BHV data.
    %       Neurophys.HeaderInformation
    %       Neurophys.Digital
    %           - CodeNumbers: {1xnTrials cell}
    %       	- CodeTimes: {1xnTrials cell}
    %       	- *SpikeTimes: [1x1 struct]
    %       	- *SpikeWaveForms: [1x1 struct]
    %
    %             *SpikeTimes,SpikeWaveForms, and Analog contain fields named using the
    %             plexon data file convention (e.g. sig001a)
    %
    %       obj.Neurophys.Analog contains the analog data which can include eye
    %       position data, joystick data, local field potentials (LFP) etc.
    %           obj.Neurophys.Analog.AD01
    %               adfreq: 1000
    %               nADsamples: 6047781
    %               tsFirstADsample: 0.0101
    %               Voltage: {1x880 cell}
    %
    % erik.emeric@gmail.com
    
    properties (SetAccess = private)
        % Struct with fields for events, analog, and spike data
        % from a Plexon or TDT datafile
        % EVENT DATA
        % CONTINUOUS DATA (ANALOG)
        % Neurophys.Analog.AD01
        
        % Neurophys.HeaderInfo
        HeaderInformation
        % NeuroPhys.Digital.CodeNumbers
        % NeuroPhys.Digital.CodeTimes
        % NeuroPhys.Digital.SpikeTimes
        %   NeuroPhys.Digital.SpikeTimes.sig001i
        %   NeuroPhys.Digital.SpikeTimes.sig001a
        Digital
        Analog
        % NeuroPhys.Analog.AD01
        %         ExperimentNotes
    end
    
    methods
        % Constructor
        function obj = NeuroPhysObject(NeuroPhysStruct)
        %NeuroPhysObject NeuroPhysObject is the constructor function of the
        % class NeuroPhysObject.
        % INPUT ARGUMENTS
        % The NPstruct returned by Plx2matMonkeyLogic.m has the following fields
        % and is organized in similar fashion as the digital events in the
        % BHV data.
        % NeuroPhysStruct.HeaderInformation
        % NeuroPhysStruct.Digital
        %   - CodeNumbers:      {1xnTrials cell}
        % 	- CodeTimes:        {1xnTrials cell}
        % 	- *SpikeTimes:      [1x1 struct]
        % 	- *SpikeWaveForms:  [1x1 struct]
        %
        %   *SpikeTimes,SpikeWaveForms, and Analog contain fields named using the
        %   plexon data file convention (e.g. sig001a)
        %
        % NeuroPhysStruct.Analog contains the analog data which can include eye
        % position data, joystick data, local field potentials (LFP) etc.
        % NeuroPhysStruct.Analog.AD01
        %   - adfreq: 1000
        %   - nADsamples: 6047781
        %   - tsFirstADsample: 0.0101
        %   - Voltage: {1x nTrial cell array}
        %
        % erik.emeric@gmail.com
        %
        
            if isequal(fieldnames(NeuroPhysStruct) ...
                    ,{'HeaderInformation';'Digital';'Analog'})
                obj.HeaderInformation 	= NeuroPhysStruct.HeaderInformation;
                obj.Digital 			= NeuroPhysStruct.Digital;
                obj.Analog 				= NeuroPhysStruct.Analog;
            else
                error('NeuroPhys struct format not recognized')
            end
            % to do
            % import contents of text files and screenshots
            
        end
        %getSpikes(whichSpike)
        function [spikeList,spikes] = getSpikes(obj,whichSpike)
        %getSpikes The getSpikes method of the NeuroPhysObject class
        % returns a cell array of spike timestamp matricies stored in
        % the NeuroPhysObject.
        %
        % INPUT ARGUMENT
        % - whichSpike, can be a cell array of spike names (e.g.,
        % whichSpike = {'sig001a';'sig001b'};) or an integer array which is
        % used as the index of the desired neuron
        % * if no input argument is provided, all of the spikes from all of
        % the neurons are returned.
        %
        % OUTPUT ARGUMENTS
        % - spikeList,list of neuron IDs (e.g., sig001a)
        % - spikes, whichSpike neurons x {nTrials x max number of trial spikes}
        %
        % erik.emeric@gmail.com
            
            if exist('whichSpike') && ~isempty(whichSpike)
                % is whichSpike a cell array of strings?
                if iscellstr(whichSpike)
                    spikeList = whichSpike;
                elseif isnumeric(whichSpike)
                    allSpikes = fieldnames(obj.Digital.SpikeTimes);
                    spikeList = allSpikes(whichSpike);
                end
            else
                spikeList = fieldnames(obj.Digital.SpikeTimes);
            end
            
            disp(spikeList)
            for curSpike =1:length(spikeList)
                eval(['SpikeTimes = obj.Digital.SpikeTimes.',spikeList{curSpike},';'])
                spks=[];
                for curTrl = 1:length(SpikeTimes)
                    spks(curTrl,1:length(SpikeTimes{curTrl}))=SpikeTimes{curTrl};
                end
                % replace 0s with nans
                spks(spks==0)=nan;
                spikes{curSpike,1}=spks;
            end
        end
        %getWFs
        function [spikeList,WaveForms] = getWFs(obj,whichSpike)
            %getWFs The getWFs method of the NeuroPhysObject class
            % returns a cell array of spike waveforms stored in
            % the NeuroPhysObject.
            %
            % INPUT ARGUMENT
            % - whichSpike, can be a cell array of spike names (e.g.,
            % whichSpike = {'sig001a';'sig001b'};) or an integer array which is
            % used as the index of the desired neuron
            % * if no input argument is provided, all of the spikes from all of
            % the neurons are returned.
            %
            % OUTPUT ARGUMENTS
            % - spikeList,list of neuron IDs (e.g., sig001a)
            % - WaveForms, whichSpike neurons x {nWaveForms x n points per waveform}
            %
            % erik.emeric@gmail.com
            if exist('whichSpike') && ~isempty(whichSpike)
                % is whichSpike a cell array of strings?
                if iscellstr(whichSpike)
                    spikeList = whichSpike;
                elseif isnumeric(whichSpike)
                    allSpikes = fieldnames(obj.Digital.SpikeTimes);
                    spikeList = allSpikes(whichSpike);
                end
            else
                spikeList = fieldnames(obj.Digital.SpikeTimes);
            end
            
            disp(spikeList)            
            for curSpike =1:length(spikeList)
                eval(['SpikeWFs = obj.Digital.SpikeWaveForms.',spikeList{curSpike},';'])
                WFs=[];
                for curTrl = 1:length(SpikeWFs)
                    WFs=[WFs;SpikeWFs{curTrl}];%(curTrl,1:length(SpikeWFs{curTrl}))=SpikeWFs{curTrl};
                end
                % store current WF in output variable WaveForms                
                WaveForms{curSpike,1}=WFs;
            end
            
        end
        %getEvents
        function [adList,AnalogData] = getAnalog(obj,whichAD)
        %getAnalog The getAnalog method of the NeuroPhysObject class
        % returns a cell array of spike timestamp matricies stored in
        % the NeuroPhysObject.
        %
        % INPUT ARGUMENT
        % - whichAD, can be a cell array of AD channel names (e.g.,
        % whichAD = {'AD01';'AD02'};) or an integer array which is
        % used as the index of the desired analog channel(s).
        % * if no input argument is provided, all of the AD channels
        % are returned.
        %
        % OUTPUT ARGUMENTS
        % - adList,list of neuron IDs (e.g., sig001a)
        % - AnalogData, whichAD n AD channels x {nTrials x max number of trial AD samples}
        %
        % erik.emeric@gmail.com
            if exist('whichAD') && ~isempty(whichAD)
                % is whichAD a cell array of strings?
                if iscellstr(whichAD)
                    adList = whichAD;
                elseif isnumeric(whichAD)
                    allAD = fieldnames(obj.Analog);
                    adList = allAD(whichAD);
                end
            else
                adList = fieldnames(obj.Analog);
            end
            
            disp(adList)
            for curAD =1:length(adList)
                eval(['adValues = obj.Analog.',adList{curAD},'.Voltage;'])
                AD=[];
                for curTrl = 1:length(adValues)
                    temp=adValues{curTrl};
                    % replace zero values with a very small value
                    temp(temp==0)= 0.0000000000001;
                    AD(curTrl,1:length(temp)) = temp;
                end
                % replace 0s with nans
                AD(AD==0)=nan;
                AnalogData{curAD,1}=AD;
            end
        end
        %getEvents        
        function [eventValues,eventTimes] = getEvents(obj)
            %getEvents The getEvents method of the NeuroPhysObject class
            % returns the code numbers and timestamps stored in the
            % NeuroPhysObject.
            %
            % INPUT ARGUMENT
            % - n/a
            %
            % OUTPUT ARGUMENTS
            % - eventValues,
            % - eventTimes, corresponding timestamps of behavioral events
            %
            % erik.emeric@gmail.com
            
            disp('Retrieving events...')
            eventValues = []; eventTimes = [];
            for curTrl = 1:length(adValues)
                evts = obj.Digital.CodeNumbers{curTrl};
                ts   = obj.Digital.CodeTimes{curTrl};
                eventValues(curTrl,1:length(evts)) = evts;
                eventTimes(curTrl,1:length(ts))    = ts;
            end
            
        end
        
%         function [raster,histogram,pt,SDF] = getSDFRaster(obj,whichSpike,trials,event,peritime,KernelType)
%            [raster,histogram,pt,SDF]=spikeDensityFunction...
%                (spike,trials,event,peritime,KernelType); 
%         end
        
    end

end
