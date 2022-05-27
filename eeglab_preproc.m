% LFP data pre-processing pipeline:
% adapted from Sangtae Ahn (sangtae_ahn@med.unc.edu) lfp_preproc
%  
% 1. add channel info to the EEG structure
% 2. downsample (fd) if needed (lfp is already lp filtered at 300Hz, ds at
% 1000Hz)
% 3. ASR (interpolation method) (A)
% 4. epoching (e)
% 5. removing noisy trials (r)
% 6. get good trialID and save evtTimes
% 7. save validChn info as eeglab_validChn
% 
% edited by Angel Huang on July 26, 2019
% validated by Sangtae Ahn on Feb 6, 2020
% Frohlich Lab.

%%
close all;
clear;
clc;

%% 
skipRec = 0;
doPlot = 1; % to save time, doPlot =0
redo = 1;
doBeforeEpoch    = 1; % step 1-3
doCleanLine      = 0; % before step 3, optional, depend on if signal has line noise
doKeepChn        = 0; % manually filter out channel by using keepChn.m
doPlot_fdA       = 1; % step 3
doEpoching       = 1; % step 4
if redo== 1
    doBeforeEpoch    = 1; % step 1-3
    doCleanLine      = 0; % before step 3, optional, depend on if signal has line noise
    doKeepChn        = 1; % manually filter out channel by using keepChn.m
    doPlot_fdA       = 1; % step 3
    doEpoching       = 1; % step 4
end
doEpochRejection = 1; % step 5 -- manual step
doEventTimes     = 1; % step 6
doValidChn       = 1; % step 7

CodeDir = ['E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/'];
addpath(genpath([CodeDir 'Toolboxes/eeglab2019_0/']));
addpath([CodeDir 'ephys/AH_toolbox'])
regionNames = {'PFC', 'LPl', 'PPC', 'VC'};
%addpath(genpath(CodeDir));
animalCodes = {'0180','0181','0179','0171'};
baseDir = ['Z:/Ferret Data/'];
for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    
PreprocessDir = [baseDir animalCode '/Preprocessed/'];
AnalysisDir   = [baseDir animalCode '/Analyzed/'];
%fileInfo = dir([PreprocessDir animalCode '_Level8*']); %AttentionTask6* detect files to load/convert  '_LateralVideo*'
fileInfo = dir([PreprocessDir animalCode '_baseline*']); %AttentionTask6* detect files to load/convert  '_LateralVideo*'

% parameters
highCut = 300;
desiredFs = 1000;
freqRange = [1 200];
    

for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    level = splitName{2}(6);
    rootPreprocessDir = [PreprocessDir recName '/'];    

    if exist(join([rootPreprocessDir 'lfp/lfp_' num2str(desiredFs) 'fdAer_StimCorD.set']))
        fprintf('Record %s already analyzed \n',recName); 
    %else continue %<<<temp for adding .keep column (2/1/2020)
        if skipRec == 1; continue; end  % to skip already analyzed records
    end
    fprintf('Analyzing record %s \n',recName);
    
if doBeforeEpoch  == 1    
    eeglab;  % call eeglab before assigning fields, create default empty EEG structure   
    pop_editoptions('option_savetwofiles', 1,'option_single', 0); % save both fdt and set file to speed up loading
    
    % load ttl and behavioral data
    [EEG.data, EEG.srate] = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat','lfpFs');
    EEG.chanlocs = readlocs('chnMap_0171.xyz','filetype','xyz'); % default field for EEGLAB
    EEG.urchanlocs = EEG.chanlocs;
    EEG = eeg_checkset(EEG);
    %figure;topoplot([],EEG.chanlocs); % without label
    eeglab redraw

    %% spectra
if doPlot ==1; pop_eegplot(EEG,1,1,1); end % plot raw traces

%     % 2.1 filter and downsample (fd) % lfp is already filtered
%     figure; 
%     pop_spectopo(EEG, 1, [0  EEG.pnts], 'EEG' , 'freqrange',freqRange,'electrodes','off');
%     %EEG = pop_eegfiltnew(EEG, 'hicutoff',300); % lfp is already filtered
    
%     % 2.2 downsample 
%     EEG = pop_resample(EEG, desiredFs);

    % plot spectra
if doPlot == 1
    figure;
    [spec1 freq]=pop_spectopo(EEG, 1, [0  EEG.pnts], 'EEG' , 'freqrange',freqRange,'electrodes','off');
    
    fig = AH_figure_eeg(freq, spec1, regionNames);
    savefig(fig, [rootPreprocessDir 'lfp/spec_' num2str(desiredFs) 'fd.fig'],'compact');
    saveas(fig, [rootPreprocessDir 'lfp/spec_' num2str(desiredFs) 'fd.png']);
    fig = figure;plottopo(spec1(:,1:30)*-1, EEG.chanlocs); % plot 1-30Hz spectra on topoplot % don't show channel 
end
%% Optional: clean line noise
if strcmp(recName,'0180_Level7c_13_20200115')
    doCleanLine = 1;
    iChn = [1:6,9:16];
end
if doCleanLine == 1
    %select channels need to clean line noise on
    %iChn = [1:6,9:16];
    figure; plot(freq,spec1(iChn,:)');
    %noiseFreq = [31 62 94 120 156 187 219 240 281];
    noiseFreq = [60,180,300];
    % method 1: removeLineNoise_SpectrumEstimation
    EEG.data(iChn,:) = removeLineNoise_SpectrumEstimation(EEG.data(iChn,:), EEG.srate, ...
    'LF = 59.6, NH = 4, M = 1024, HW = 3'); % 
%       NH  = number of harmonics. (default: 1)
%       LF  = line frequency. (default: 50 Hz)
%       TOL = error tolerence. (default: 1% tol)
%       HW  = half-width of peak in samples. (default: 2)
%       M   = Size of window. (For fs = 1 kHz, m = 1024). If M is set, TOL
%                              has no effect.
    
    % method 2: cleanline
    %EEG = pop_cleanline(EEG, 'bandwidth',4,'chanlist',iChn ,'computepower',1,'linefreqs',noiseFreq,...
    %    'normSpectrum',0,'p',0.5,'pad',100,'plotfigures',0,'scanforlines',1,'sigtype','Channels',...
    %    'tau',100,'verb',1,'winsize',0.8,'winstep',1); % winsize =1~1.1 
        % bandwidth: width of a spectral peak for a sinusoid at fixed frequency
        % p-value for detection of significant sinusoid,default=0.01 Input Range[0  1]                  
        % winsize: Sliding window length (sec),1s, longer is slower to run and not as strong 
        % winsize=0.5 can't detect noise
    
    %method 3: notch filter -- create a big dip in targeting frequency and around, not good
    %EEG = pop_eegfiltnew(EEG, 'locutoff',59.57-1,'hicutoff',59.57+1,'revfilt',1); 
    
    % inspect result:
    [spec2]=pop_spectopo(EEG, 1, [0  EEG.pnts], 'EEG' , 'freqrange',freqRange,'electrodes','off','plot','off');
    % when only clean 1 channel:
    %hold on; plot(freq,spec2(iChn,:)); legend('before','after'); xlim([20,40]);ylim([-10,0]);
    figure; plot(freq,spec2(iChn,:)');
    AH_figure_eeg(freq,spec2,regionNames);
    
end
%pop_eegplot_w(EEG,1,1,1); 
pop_saveset(EEG,'filepath',[rootPreprocessDir 'lfp/'],'filename',['lfp_' num2str(desiredFs) 'fd.set']);

%% 3. ASR (interpolation method) Artifact Subspace Reconstruction
% https://sccn.ucsd.edu/wiki/Artifact_Subspace_Reconstruction_(ASR)
% clean_rawdata(data, FlatlineCriterion, Highpass, ChannelCriterion,
% LineNoiseCriterion, BurstCriterion, WindowCriterion)
% Option 1: no channel rejection
%EEG = clean_rawdata(EEG, 5, 'off', 'off', 'off', 5, 'off');

% Option 2: channel rejection without channel interpolation, per region for
% better outcome.
pfcChn = [1:16]';
pulChn = [1:16]' + 16;
ppcChn = [1:16]' + 32;
vcChn  = [1:16]' + 48;
EEG.etc.allChn = {pfcChn, pulChn, ppcChn, vcChn};
EEG.etc.regionNames = regionNames;
oldEEG = EEG;
EEG.data = []; % clear out EEG.data (since new EEG will have less channel)
if doKeepChn == 1 % use manually filtered channel list
    [validChn,~] = keepChn(recName);
    for iRegion = 1:numel(regionNames)
        validChnMask{iRegion} = false(16,1); %create logic array
        validChnMask{iRegion}(validChn{iRegion}-(iRegion-1)*16) = 1;
    end
end

% using ASR to clean data region by region
numPrevChn = 0;
for iRegion = 1:numel(regionNames)
    regionName = regionNames{iRegion};
    lfp.(regionName) = oldEEG;
    lfp.(regionName).data = lfp.(regionName).data(EEG.etc.allChn{iRegion},:);
    lfp.(regionName).chanlocs = lfp.(regionName).chanlocs(EEG.etc.allChn{iRegion});
    lfp.(regionName) = eeg_checkset(lfp.(regionName)); % adjust channel number in structure; necessary before clean_rawdata
    %if iRegion == 4 %VC should have less spindle
    %lfp.(regionName) = clean_rawdata(lfp.(regionName), 60*60, 'off', 0.8, 4, 12, 'off'); %flatline, ,correlation,linenoise,BurstCriterion,
    %else
    if doKeepChn == 0
        lfp.(regionName) = clean_rawdata(lfp.(regionName), 'off', 'off', 0.8, 4, 20, 'off'); %flatline, ,correlation,linenoise,BurstCriterion,        
    else
        % turn off channel rejection, and only use validChn
        if ismember(level,{'7','9'}) && (iRegion == 2 || iRegion == 3) % LPl and PPC have opto effect, use more lenient burst criteria           
           lfp.(regionName) = clean_rawdata(lfp.(regionName), 'off', 'off', 'off', 'off', 40, 'off'); %flatline, ,correlation,linenoise,BurstCriterion,        
        else % either PFC or VC normally don't have opto effect
           lfp.(regionName) = clean_rawdata(lfp.(regionName), 'off', 'off', 'off', 'off', 15, 'off'); %flatline, ,correlation,linenoise,BurstCriterion,        
        end
        lfp.(regionName).etc.clean_channel_mask = logical(validChnMask{iRegion}); % has to use logical otherwise double doesn't work as intended!
%         lfp.(regionName).etc.clean_channel_mask(validChnMask{iRegion},1) = 1; % replace good channel mask
%         lfp.(regionName).etc.clean_channel_mask(~validChnMask{iRegion},1) = 0; % replace bad channel mask
        % update data so that vis_artifacts can run properly
        lfp.(regionName).data = lfp.(regionName).data(lfp.(regionName).etc.clean_channel_mask,:);
    end
        %end
    
        
    % To compare the data before and after ASR
%     old = oldEEG; old.data = old.data(EEG.etc.allChn{iRegion},:);old = eeg_checkset(old);
%     vis_artifacts(lfp.(regionName),old)
    % To mannually reject a channel
    %lfp.(regionName).etc.clean_channel_mask = true(16,1);
    %lfp.(regionName).etc.clean_channel_mask(11) = false;
    
    % Concatenate all regions to the original EEG structure    
    if isfield(lfp.(regionName).etc, 'clean_channel_mask')
        chnMask = logical(lfp.(regionName).etc.clean_channel_mask);
    else
        chnMask = true(numel(EEG.etc.allChn{iRegion}),1); %need to be logical otherwise validChn will populate as 1
    end
    % Concat channel mask
    EEG.etc.clean_channel_mask(EEG.etc.allChn{iRegion},1) = (chnMask);    
    EEG.etc.validChn{iRegion} = EEG.etc.allChn{iRegion}(logical(chnMask));
    EEG.etc.reorderedChn{iRegion} = [1:length(EEG.etc.validChn{iRegion})]'+numPrevChn;
    EEG.data(EEG.etc.reorderedChn{iRegion},:) = lfp.(regionName).data; % concatenate data from all regions, save in new EEG struct
    numPrevChn = numPrevChn + length(EEG.etc.validChn{iRegion});
end
EEG.chanlocs = EEG.urchanlocs(logical(EEG.etc.clean_channel_mask));
EEG.etc.removedCh = find(EEG.etc.clean_channel_mask==0);
EEG = eeg_checkset(EEG);
clear lfp % to save memory


if doPlot_fdA == 1
    [spec3 freq] = pop_spectopo(EEG, 1, [0  EEG.pnts], 'EEG' , 'freqrange',freqRange,'electrodes','off');
    fig = figure;
    for iRegion = 1 : numel(regionNames)
        subplot(2,2,iRegion);
        plot(freq,spec3(EEG.etc.reorderedChn{iRegion},:)');
        title(regionNames{iRegion});
        xlabel('Frequency [Hz]');
        ylabel('Power [dB]');    
        xlim([0,150]);
    end
    savefig(fig, [rootPreprocessDir 'lfp/spec_' num2str(desiredFs) 'fdA.fig'],'compact');
    saveas(fig, [rootPreprocessDir 'lfp/spec_' num2str(desiredFs) 'fdA.png']);
end
pop_saveset(EEG,'filepath',[rootPreprocessDir 'lfp/'],'filename',['lfp_' num2str(desiredFs) 'fdA.set']);

% % parameters 
% FlatlineCriterion=5; % Maximum tolerated flatline duration
% Highpass= -1; % disabled
% ChannelCriterion = 0.8; % Minimum channel correlation
% LineNoiseCriterion = 4;
% BurstCriterion = 15;  % default 10-100, Standard deviation cutoff for
% removal of bursts, lower reject more burst
% (interpolation based on neighboring time period)
% WindowCriterion = -1; % disabled
% 
% EEG = clean_rawdata(EEG,FlatlineCriterion,Highpass,ChannelCriterion,LineNoiseCriterion,BurstCriterion,WindowCriterion); % default setting
% pop_eegplot(EEG,1,1,1); % plot raw traces
end

% 4. epoching data
% load event info
if doEpoching == 1
    % big screen visualization of EEG
    if ~exist('EEG'); EEG = pop_loadset([rootPreprocessDir 'lfp/lfp_' num2str(desiredFs) 'fdA.set']);end
    sessionMetaBehav = is_load([rootPreprocessDir 'sessionMetaBehav.mat'], 'sessionMetaBehav');
    EEG.urevent = makeEvent_for_eeglab(sessionMetaBehav);
    EEG.event = EEG.urevent;
    try
    EEG = pop_epoch(EEG, {'StimCorD4' 'StimCorD5' 'StimCorD6'}, [-8 5], 'epochinfo', 'yes');
    pop_saveset(EEG,'filepath',[rootPreprocessDir 'lfp/'],'filename',['lfp_' num2str(desiredFs) 'fdA_StimCorD.set']);
    catch
    fprintf('Record %s no correct trials \n',recName); 
    end
end



% 5. manually trial rejection
if doEpochRejection == 1 
    desiredFs = 1000;
    if ~exist('EEG')
        EEG = pop_loadset([rootPreprocessDir 'lfp/lfp_' num2str(desiredFs) 'fdA_StimCorD.set']);
    	eeglab redraw % required with previous line, otherwise eeglab error
    end
    pop_eegplot_w(EEG,1,1,1); % wider than pop_eegplot(EEG,1,1,1);
    fprintf('Analyzing record %s \n',recName);
    input(''); % to pause until a key is pressed
    % EEG = pop_select( EEG, 'nochannel',{'''PFC11'''}); % to manually reject selected channel

    
    pop_saveset(EEG,'filepath',[rootPreprocessDir 'lfp/'],'filename',['lfp_' num2str(desiredFs) 'fdAer_StimCorD.set']);
    %EEG = pop_rejepoch(EEG, 47,0); % eg. if reject trial 28
    %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
end

if doEventTimes == 1    
    sessionMetaBehav = is_load([rootPreprocessDir 'sessionMetaBehav.mat'], 'sessionMetaBehav');
    if ~exist('EEG')
        EEG = pop_loadset([rootPreprocessDir 'lfp/lfp_' num2str(desiredFs) 'fdAer_StimCorD.set']);
    end
    level = sessionMetaBehav.SessionType(1,6);
    [alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
    alignNames = {'Stim'}; %'Init',
    for ialign = 1:numel(alignNames)
        alignName = alignNames{ialign};
        
        hitMissName = 'Cor';
        if level == '6'
            condNames = delayNames;
            condID = [1,2,3]; % don't include the last condition with all combined
        else % level 7 8 9
            condNames = optoNames;
            condID = [1,2,3,4,5]; % don't include the last condition with all combined
        end
        numCond = numel(condID);

        for iCond = 1:numCond+1
            evtTimes{iCond} = [];
            trialIDs{iCond} = [];
            epochIDs{iCond} = [];
            baseTwins{iCond} = [];
            delays{iCond} = [];        
        end
        if level == '6'
            for i = 1:numel(EEG.event)
                name = EEG.event(i).type;
                trialAlignName = EEG.event(i).type(1:4);
                hitMissName = EEG.event(i).type(5:7);
                condName = EEG.event(i).type(8:9);
                iCond = find(strcmp(condName, condNames));
                delayDuration = str2num(condName(2)); 
                if strcmp(hitMissName, 'Cor') && strcmp(alignName, trialAlignName)
                    if strcmp(alignName, 'Stim'); columnName = 'StimOnset';
                    elseif strcmp(alignName, 'Init'); columnName = 'Init';end
                    
                    ID = EEG.event(i).urevent;
                    sessionMetaBehav.keep(ID) = 1;
                    evtTimes{iCond} = [evtTimes{iCond} sessionMetaBehav.(columnName)(ID)]; % convert to sec
                    trialIDs{iCond} = [trialIDs{iCond} EEG.event(i).urevent]; 
                    epochIDs{iCond} = [epochIDs{iCond} EEG.event(i).epoch];

                    % all conditions combined
                    evtTimes{numCond+1} = [evtTimes{numCond+1} sessionMetaBehav.(columnName)(ID)];
                    trialIDs{numCond+1} = [trialIDs{numCond+1} EEG.event(i).urevent]; 
                    epochIDs{numCond+1} = [epochIDs{numCond+1} EEG.event(i).epoch];
                    
                    if strcmp(alignName, 'Stim')
                        twins{iCond} = [-8,5]; 
                        twins{numCond+1} = [-8,5];
                        baseTwins{iCond} = [-2,-0.5] - delayDuration;
                        baseTwins{numCond+1} = [baseTwins{numCond+1}; [-2,-0.5] - delayDuration];
                        delays{iCond} = delayDuration;
                        delays{numCond+1} = [delays{numCond+1}; delayDuration]; 
                    elseif strcmp(alignName, 'Init')
                        twins{iCond} = [-2,10]; 
                        twins{numCond+1} = [-2,10];                        
                        baseTwins{iCond} = [-2,-0.5];
                        baseTwins{numCond+1} = [baseTwins{numCond+1}; [-2,-0.5]];
                        delays{iCond} = delayDuration;
                        delays{numCond+1} = [delays{numCond+1}; delayDuration]; 
                    end
                end
            end
        save([rootPreprocessDir 'eventTimes_' alignName 'Cor.mat'], 'evtTimes','trialIDs','epochIDs','twins','baseTwins','delays','condNames');
        else % level 7 8 9 
        for i = 1:numel(EEG.event)
            name = EEG.event(i).type;
            trialAlignName = EEG.event(i).type(1:4);
            hitMissName = EEG.event(i).type(5:7);
            ID = EEG.event(i).urevent;
            condName = sessionMetaBehav(ID,:).OptoType;
            iCond = find(strcmp(condName, condNames));
            delayDuration = str2num(name(end)); 
            if strcmp(hitMissName, 'Cor') && ~isempty(iCond) && strcmp(alignName, trialAlignName)% exclude "On" trials without freq               
                if strcmp(alignName, 'Stim'); columnName = 'StimOnset';
                elseif strcmp(alignName, 'Init'); columnName = 'Init';end
                sessionMetaBehav.keep(ID) = 1;
                evtTimes{iCond} = [evtTimes{iCond} sessionMetaBehav.(columnName)(ID)]; % convert to sec
                trialIDs{iCond} = [trialIDs{iCond} EEG.event(i).urevent]; 
                epochIDs{iCond} = [epochIDs{iCond} EEG.event(i).epoch];

                % all conditions combined
                evtTimes{numCond+1} = [evtTimes{numCond+1} sessionMetaBehav.(columnName)(ID)];
                trialIDs{numCond+1} = [trialIDs{numCond+1} EEG.event(i).urevent]; 
                epochIDs{numCond+1} = [epochIDs{numCond+1} EEG.event(i).epoch];

                if strcmp(alignName, 'Stim')
                    twins{iCond} = [-8,5]; 
                    twins{numCond+1} = [-8,5];
                    baseTwins{iCond} = [baseTwins{iCond}; [-2, -0.5] - delayDuration];
                    baseTwins{numCond+1} = [baseTwins{numCond+1}; [-2, -0.5] - delayDuration];
                    delays{iCond} = [delays{iCond}; delayDuration]; 
                    delays{numCond+1} = [delays{numCond+1}; delayDuration];
                elseif strcmp(alignName, 'Init')
                    twins{iCond} = [-2,10]; 
                    twins{numCond+1} = [-2,10];
                    baseTwins{iCond} = [baseTwins{iCond}; [-2,-0.5]];
                    baseTwins{numCond+1} = [baseTwins{numCond+1}; [-2,-0.5]];
                    delays{iCond} = [delays{iCond}; delayDuration]; 
                    delays{numCond+1} = [delays{numCond+1}; delayDuration];
                end
            end
        end
        save([rootPreprocessDir 'optoEventTimes_' alignName 'Cor.mat'], 'evtTimes','trialIDs','epochIDs','twins','baseTwins','delays','condNames');        
    end
    save([rootPreprocessDir 'sessionMetaBehav.mat'], 'sessionMetaBehav'); % added one column of keepMask
    end
end

if doValidChn == 1
    %eeglab;
    %EEG = pop_loadset([rootPreprocessDir 'lfp/lfp_' num2str(desiredFs) 'fdAer_StimCorD.set']);
    lfp = EEG.etc; 
    for iRegion = 1:numel(regionNames)
        if lfp.validChn{iRegion}(1) == lfp.validChn{iRegion}(2)
            lfp.validChn{iRegion} = [1:16]' + (iRegion-1)*16;
            EEG.etc.validChn{iRegion} = lfp.validChn{iRegion};
        end
    end
    save([rootPreprocessDir 'eeglab_validChn.mat'], 'lfp');
    pop_saveset(EEG,'filepath',[rootPreprocessDir 'lfp/'],'filename',['lfp_' num2str(desiredFs) 'fdAer_StimCorD.set']);
end
clear EEG
end % end of record
end % end of animal