% function plotBadTrials(analogDataFull,powerVsFreq)
%%%%%%%% Plotting the time and frequency domain waveforms/PSDs %%%%%%%%%%%

checkTheseElectrodes = [24 26 29 30 31 57 58 61 62 63];k =0;
% Set MT parameters
params.tapers   = [3 5];
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 200];
params.trialave = 0;

% Time domain waveforms
% clear meanTrialData stdTrialData tDplusFreq tDminusFreq
% numTrials = size(analogDataFull,2);
% meanTrialData = nanmean(analogDataFull,2);                     % mean trial trace
% stdTrialData = nanstd(analogDataFull,[],2);                    % std across trials
% 
% tDplusTime = (meanTrialData + (thresholdTime)*stdTrialData);    % upper boundary/criterion
% tDminusTime = (meanTrialData - (thresholdTime)*stdTrialData);
% figure;
% for iElec = 1:numElectrodes
%     iElec = 58;
%     if find(checkTheseElectrodes==iElec)
%         k= k+1;
%         subplot(5,2,k);
%         plot(timeVals(checkPeriodIndices),squeeze(tDplusTime(iElec,:,:)),'lineWidth',3); hold on; plot(timeVals(checkPeriodIndices),squeeze(tDminusTime(iElec,:,:)),'lineWidth',3); hold on;
%         title(['electrode:',EEGelectrodeLabels{iElec+1,2}]); hold on;
%         xlabel('Time in s'); ylabel('Voltage in uV'); axis tight;
%         if isempty(plotBadTrialsTime{iElec}); continue;end
%         plot(timeVals(checkPeriodIndices),squeeze(analogDataFull(iElec,plotBadTrialsAll{iElec},:))); hold on;  %end; hold on;
%     end
% end
% k=0;
% %PSD plots
% figure;
% for iElec=1:numElectrodes
%     
%     clear blPowerVsFreq freqVals
%     [powerVsFreq,freqVals] = mtspectrumc(squeeze(analogDataFull(iElec,:,:))',params);
%     powerVsFreq = powerVsFreq';
%     clear meanTrialData stdTrialData tDplusFreq tDminusFreq
%     meanTrialData = nanmean(powerVsFreq,1);                     % mean trial trace
%     stdTrialData = nanstd(powerVsFreq,[],1);                    % std across trials
%     
%     tDplusFreq = (meanTrialData + (thresholdFreq)*stdTrialData);    % upper boundary/criterion
%     tDminusFreq = (meanTrialData - (thresholdFreq)*stdTrialData);
%     if find(checkTheseElectrodes==iElec)
%         k= k+1;
%         subplot(5,2,k);
%         plot(freqVals,log10(tDplusFreq),'lineWidth',3); hold on; plot(freqVals,log10(tDminusFreq),'lineWidth',3);
%         title(['electrode:',EEGelectrodeLabels{iElec+1,2}]); hold on;% text(80,50,['Threshold used: ',num2str(thresholdFreq),'sigma']);
%         xlabel('Frequency in Hz'); ylabel('Power in dB'); xlim([0 100]); axis tight;
%         if isempty(plotBadTrialsPSD{iElec}); continue;end
%         plot(freqVals,log10(powerVsFreq(plotBadTrialsAll{iElec},:))); %end; hold on;
%     end
% end
badTrials=[];
badTrialUL = (badTrialThreshold/100)*numTrials;
for iElec=1:numElectrodes
    if length(plotBadTrialsAll{iElec})>badTrialUL; plotBadTrialsAll{iElec} = NaN; badElec(iElec)=1; else badElec(iElec)=0;end % Discard electrodes based on %-age of bad trials
    badTrials=union(badTrials,plotBadTrialsAll{iElec}); % in the previous case we took the union
end
badTrials(isnan(badTrials))=[];
totalTrials = size(analogDataAllElecs,2);
% Counting the trials which occurs in more than x% of the electrodes
for iTrial = 1:totalTrials
    badTrialElecs(iTrial)=0;
    for iElec = 1:numElectrodes
        if isnan(allBadTrials{1,iElec}); badTrialElecs(iTrial)= 0;continue; end % Discarding the electrodes where the bad trials are NaN
        if find(badTrials(iTrial)==allBadTrials{1,iElec})
            badTrialElecs(iTrial) = badTrialElecs(iTrial)+1;
        end
    end
end
tempBadTrials = originalTrialInds(badTrials); 
badTrials(badTrialElecs<(badElecThreshold/100.*numElectrodes))=[];
% badTrials = originalTrialInds(badTrials); %Remapping 
% plot of trial vs electrodes

p= zeros(numElectrodes,totalTrials);
for iElec = 1:numElectrodes
    if GoodElec_Z(iElec)==0; p(iElec,:) = 3; end
    if badElec(iElec)==1; p(iElec,:) = 7; end
    if isnan(allBadTrials{iElec}); continue; end
    if ~isempty(allBadTrials{iElec}); p(iElec,allBadTrials{iElec})= 9; end
end
for  iTrial = 1:totalTrials
    if find(badTrials==iTrial)
        p(:,iTrial)=11;
    end
end
figure;
imagesc(1:totalTrials,1:numElectrodes,p);

%{
if exist('highPassCutOff','var') || ~isempty(highPassCutOff)
    d1 = designfilt('highpassiir','FilterOrder',8, ...
        'PassbandFrequency',highPassCutOff,'PassbandRipple',0.2, ...
        'SampleRate',Fs);
end

for iElec = 1:numElectrodes
    
    analogData = load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeNum) '.mat']),'analogData');
    analogData = analogData.analogData;
    
       
    if exist('highPassCutOff','var') || ~isempty(highPassCutOff)
        clear analogDataSegment; analogDataSegment = squeeze(analogData);
        clear analogData; analogData = filtfilt(d1,analogDataSegment')';
        clear analogDataSegment
    end
    
    analogDataAllElecs(iElec,:,:) = analogData;
    analogData(badEyeTrials,:) = [];
    % subtract dc
    
    analogDataAllElecs = analogDataAllElecs - repmat(mean(analogDataAllElecs,2),1,size(analogDataAllElecs,2));
    numTrials = size(analogData,1);
    meanTrialData = nanmean(analogData,1);                     % mean trial trace
    stdTrialData = nanstd(analogData,[],1);                    % std across trials
    
    tDplusTime = (meanTrialData + (thresholdTime)*stdTrialData);    % upper boundary/criterion
    tDminusTime = (meanTrialData - (thresholdTime)*stdTrialData);   % lower boundary/criterion
    
    tBoolTrialsTime = sum((analogData > ones(numTrials,1)*tDplusTime) | (analogData < ones(numTrials,1)*tDminusTime(iElec)),2);
    clear tmpBadTrialsTime;
    tmpBadTrialsTime = find(tBoolTrials>0);
    %PSD
    clear blPowerVsFreq freqVals
%     params.trialav = 0;
    [powerVsFreq,freqVals] = mtspectrumc(squeeze(analogDataAllElecs(iElec,:,checkPeriodIndices))',params);
    powerVsFreq = powerVsFreq';
    
    clear meanTrialData stdTrialData tDplus tDminus
    meanTrialData = nanmean(powerVsFreq,1);                     % mean trial trace
    stdTrialData = nanstd(powerVsFreq,[],1);                    % std across trials
    
    tDplusFreq(iElec) = (meanTrialData + (thresholdFreq)*stdTrialData);    % upper boundary/criterion
    tDminusFreq(iElec) = (meanTrialData - (thresholdFreq)*stdTrialData);   % lower boundary/criterion
    
    clear tBoolTrials; tBoolTrials = sum((powerVsFreq > ones(numTrials,1)*tDplusFreq(iElec)) | (powerVsFreq < ones(numTrials,1)*tDminusFreq(iElec)),2);
    clear tmpBadTrialsPSD;
    tmpBadTrialsPSD(iElec) = find(tBoolTrials>0);
    tmpBadTrialsPSD(iElec) = setdiff(tmpBadTrialsPSD(iElec),tmpBadTrialsTime(iElec));

end
figure; k = 0;
for i= 1:numElectrodes
    if find(checkTheseElectrodes==iElec)
        k= k+1;
        subplot(5,2,k);
        plot(timeVals(checkPeriodIndices),tDplusTime{iElec},'lineWidth',3); hold on; plot(timeVals(checkPeriodIndices),tDminusTime{iElec},'lineWidth',3); hold on;
        title(['electrode:',EEGelectrodeLabels{iElec+1,2}]); hold on;
        xlabel('Time in s'); ylabel('Voltage in uV'); axis tight;
        if isempty(allBadTrialsTime{iElec}); continue;end
        plot(timeVals(checkPeriodIndices),squeeze(analogData(checkTheseElectrodes(iElec),allBadTrialsTime{checkTheseElectrodes(iElec)},checkPeriodIndices))); hold on;  %end; hold on;
    end
end
 figure; k = 0;
for iElec = 1:numElectrodes
     if find(checkTheseElectrodes==iElec)
        k= k+1;
    subplot(5,2,k);
    plot(freqVals,log10(tDplusFreq),'lineWidth',3); hold on; plot(freqVals,log10(tDminusFreq),'lineWidth',3);
    title(['electrode:',EEGelectrodeLabels{checkTheseElectrodes(iElec)+1,2}]); hold on;% text(80,50,['Threshold used: ',num2str(thresholdFreq),'sigma']);
    xlabel('Frequency in Hz'); ylabel('Power in dB'); xlim([0 100]); axis tight;
    if isempty(allBadTrialsFreq{checkTheseElectrodes(iElec)}); continue;end
    plot(freqVals,log10(powerVsFreq(allBadTrialsFreq{checkTheseElectrodes(iElec)},:))); %end; hold on;
    end
end
%}
% end

