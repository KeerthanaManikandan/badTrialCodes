
% function [analogDataFull,powerVsFreq,plotBadTrialsTime,plotBadTrialsPSD,numElectrodes,timeVals,freqVals,checkPeriodIndices]=findBadTrialsWithEEG

clear; clc;

folderIn = 'E:\data\ADGammaProjectRawData';

% subjectName = 'A054R_F0';
% expDate = '041117';

subjectName = '181VM';
expDate = '211217';
gridType = 'EEG';
montageLabels = 'actiCap64';

subjectFolder = fullfile('E:\data\human\EEG\extractedData\ADGammaProject',subjectName,gridType);

protocolName = 'GAV_0003';
if strcmpi(montageLabels,'actiCap64');channelsNotForConsideration = [65 66];end % TODO: [33 34]
if strcmpi(montageLabels,'actiCap31Posterior'); channelsNotForConsideration = [33 34];end
highPassCutOff = 1.6; % Hz
FsEye = 500; % Hz
checkPeriod = [-0.500 0.750]; % s
impedanceTag = 'Start'; % Start or End
ImpedanceCutOff = 25; % KOhm
thresholdTime = 6;%
thresholdFreq  = 6; %  %64; 66
badTrialThreshold = 20; % Percentage % 10, 20, 30
badElecThreshold = 5; % Percentage  % 10, 20, 30
checkTheseElectrodes = [24 26 29 30 31 57 58 61 62 63]; % P3,O1,Oz,O2,P1,P2,PO3,POz,Po4

folderName = fullfile(subjectFolder,expDate,protocolName);
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');

lfpInfo = load(fullfile(folderLFP,'lfpInfo'),'analogChannelsStored','timeVals');
analogChannelsStored = lfpInfo.analogChannelsStored;
timeVals = lfpInfo.timeVals;

if exist('channelsNotForConsideration','var') && ~isempty(channelsNotForConsideration)
    analogChannelsStored(ismember(analogChannelsStored,channelsNotForConsideration))=[];
end
Fs = 1/(timeVals(2) - timeVals(1));

numElectrodes = length(analogChannelsStored);

if exist('highPassCutOff','var') || ~isempty(highPassCutOff)
    d1 = designfilt('highpassiir','FilterOrder',8, ...
        'PassbandFrequency',highPassCutOff,'PassbandRipple',0.2, ...
        'SampleRate',Fs);
end

% 1. get bad trial ist from eye trials
if exist('FsEye','var') && ~isempty(FsEye)
    badEyeTrials = findBadTrialsFromEyeData(subjectFolder,expDate,protocolName,FsEye,checkPeriod)'; % added by MD 10-09-2017
else
    badEyeTrials = [];
end

% Get electrode impedances
filenameAtStart = fullfile(folderIn,[subjectName expDate],[subjectName expDate 'GAV_ImpedanceAt' impedanceTag '.txt']);
[elecNames,elecImpedance] = getImpedanceEEG(filenameAtStart);
EEGelectrodeLabels = load('EEGelectrodeLabels.mat','EEGelectrodeLabels');
EEGelectrodeLabels = EEGelectrodeLabels.EEGelectrodeLabels;
electrodeLabelsList = EEGelectrodeLabels(2:end,(ismember(EEGelectrodeLabels(1,:),'actiCap64')));

clear elecInds
for iML = 1:length(electrodeLabelsList)
    elecInds(iML) = find(strcmp(electrodeLabelsList(iML),elecNames));
end
elecImpedance = elecImpedance(elecInds);
GoodElec_Z = elecImpedance<ImpedanceCutOff;

% Set MT parameters
params.tapers   = [3 5];
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 200];
params.trialave = 0;

% 3. Analysis for each trial and each electrode
allBadTrials = cell(1,numElectrodes);
allBadTrialsTime = cell(1,numElectrodes); allBadTrialsFreq = cell(1,numElectrodes);

hW1 = waitbar(0,'Processing electrodes...');

for iElec=1:numElectrodes
    %     iElec=24;
    
    if ~GoodElec_Z(iElec); allBadTrials{iElec} = NaN; continue; end % Here we want to save computational time by analysing only those electrodes with good impedance
    
    clear analogData
    electrodeNum=analogChannelsStored(iElec);
    waitbar((iElec-1)/numElectrodes,hW1,['Processing electrode: ' num2str(electrodeNum)]);
    
    analogData = load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeNum) '.mat']),'analogData');
    analogData = analogData.analogData;
    
    analogDataAllElecs(iElec,:,:) = analogData;
    analogData(badEyeTrials,:) = [];
    
    if exist('highPassCutOff','var') || ~isempty(highPassCutOff)
        clear analogDataSegment; analogDataSegment = analogData;
        clear analogData; analogData = filtfilt(d1,analogDataSegment')';
        clear analogDataSegment
    end
    
    % determine indices corresponding to the check period
    checkPeriodIndices = timeVals>=checkPeriod(1) & timeVals<=checkPeriod(2);
    analogData = analogData(:,checkPeriodIndices);
    
    % subtract dc
    analogData = analogData - repmat(mean(analogData,2),1,size(analogData,2));
    analogDataFull(iElec,:,:) = analogData;
    
    
    % Check time-domain waveforms
    numTrials = size(analogData,1);
    meanTrialData = nanmean(analogData,1);                     % mean trial trace
    stdTrialData = nanstd(analogData,[],1);                    % std across trials
    
    tDplusTime = (meanTrialData + (thresholdTime)*stdTrialData);    % upper boundary/criterion
    tDminusTime = (meanTrialData - (thresholdTime)*stdTrialData);   % lower boundary/criterion
    
    tBoolTrials = sum((analogData > ones(numTrials,1)*tDplusTime) | (analogData < ones(numTrials,1)*tDminusTime),2);
    
    clear tmpBadTrialsTime;
    tmpBadTrialsTime = find(tBoolTrials>0);
    
    
    % Check PSDs
    clear blPowerVsFreq freqVals
    [powerVsFreq,freqVals] = mtspectrumc(analogData',params);
    powerVsFreq = powerVsFreq';
    
    clear meanTrialData stdTrialData tDplus tDminus
    meanTrialData = nanmean(powerVsFreq,1);                     % mean trial trace
    stdTrialData = nanstd(powerVsFreq,[],1);                    % std across trials
    
    tDplusFreq = (meanTrialData + (thresholdFreq)*stdTrialData);    % upper boundary/criterion
    tDminusFreq = (meanTrialData - (thresholdFreq)*stdTrialData);   % lower boundary/criterion
    
    clear tBoolTrials; tBoolTrials = sum((powerVsFreq > ones(numTrials,1)*tDplusFreq) | (powerVsFreq < ones(numTrials,1)*tDminusFreq),2);
    
    clear tmpBadTrialsPSD;
    tmpBadTrialsPSD = find(tBoolTrials>0);
    tmpBadTrialsPSD = setdiff(tmpBadTrialsPSD,tmpBadTrialsTime); % removing the common bad trials that are already present in tmpBadTrialsTime
    
    tmpBadTrialsAll = unique([tmpBadTrialsTime;tmpBadTrialsPSD]);
    
    % Remap bad trail indices to original indices
    originalTrialInds = 1:size(analogDataAllElecs,2);
    originalTrialInds(badEyeTrials) = [];
    allBadTrials{iElec} = originalTrialInds(tmpBadTrialsAll);
    allBadTrialsTime{iElec}= originalTrialInds(tmpBadTrialsTime);
    allBadTrialsFreq{iElec} = originalTrialInds(tmpBadTrialsPSD);
    
    plotBadTrialsPSD{iElec} = tmpBadTrialsPSD';
    plotBadTrialsTime{iElec} = tmpBadTrialsTime';
    plotBadTrialsAll{iElec}  = tmpBadTrialsAll';
    
    
end
close(hW1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% badTrials=allBadTrials{1};
% badTrialUL = (badTrialThreshold/100)*numTrials;
% for iElec=1:numElectrodes
%     if allBadTrials{iElec}>badTrialUL; allBadTrials{iElec} = NaN; end % Discard electrodes based on %-age of bad trials
%     trialWise{iElec} = allBadTrials{iElec};
%     badTrials=union(badTrials,allBadTrials{iElec}); % in the previous case we took the union
% end
% badTrials(isnan(badTrials))=[];
% 
% % Counting the trials which occurs in more than x% of the electrodes
% for iTrial = 1:length(badTrials)
%     badTrialElecs(iTrial)=0;
%     for iElec = 1:numElectrodes
%         if isnan(allBadTrials{1,iElec}); continue; end % Discarding the electrodes where the bad trials are NaN
%         if find(badTrials(iTrial)==allBadTrials{1,iElec})
%             badTrialElecs(iTrial) = badTrialElecs(iTrial)+1;
%         end
%     end
% end
% badTrials(badTrialElecs<(badElecThreshold/100.*numElectrodes))=[];


% end