function [TonalBurstData, TonalBurstBinData, removed, removed_ind] = ASSR_avg_LFP(LFP, EEG, ASSRblock, order, freq, events, cycle_num, timestamps, TS)
%Description: Filters and segments all probe and EEG data 
%Inputs: 
%LFP - NxM data segment with N corresponding to the number of
%data points of the recording and M corresponding to the number of LFP
%channels.
%EEG - Nx1 data segment with N corresponding to the number of
%data points of the recording. This will need to be modified when
%collecting multiple EEG channels.
%ASSRblock - The isolated parameter file for the ASSR condition.
%order - The LFP channel order corresponding to the probe used for the
%subject.
%freq - The ASSR frequency desired to be analyzed.
%events - The time markers (preferably the analog time marker)
%corresponding to the start of the ASSR event and the end of the ASSR
%event.
%cycle_num - The number of cycles desired to 

%Start by applying a notch filter to the data to remove 60Hz noise, EEG
%is fed in seperately as FP20 but is combined later in TonalBurst
%variable 
removed = [];
removed_ind = [];
Fs = 1000;
d = designfilt('bandstopiir','FilterOrder',2, ...
           'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
           'DesignMethod','butter','SampleRate',Fs);
block = find(ASSRblock == freq);
l_eeg = numel(EEG(1,:));
l_lfp = numel(LFP(1,:));
for i = 1:l_lfp
    no60LFP = filtfilt(d,LFP(:,order(i)));
    [tradfilt, ~, ~] = Butterworth_Hilbert_LR(no60LFP, Fs, [freq-5,freq+5]);
    [AvgData_LFP(:,i), LFP_bin(:,:,i)] = AnalyzeTonalBurst_CycleAverage(tradfilt, [events(block*2-1), events(block*2)], freq, 0, [0 0], cycle_num, Fs, timestamps, TS);
end
for j = 1:l_eeg
    no60EEG = filtfilt(d, EEG(:,j));
    [tradfilt_EEG, ~, ~] = Butterworth_Hilbert_LR(no60EEG, Fs, [freq-5,freq+5]);
    [AvgData_EEG(:,j), EEG_bin(:,:,j)] = AnalyzeTonalBurst_CycleAverage(tradfilt_EEG, [events(block*2-1), events(block*2)], freq, 0, [0 0], cycle_num, Fs, timestamps, TS);
end
TonalBurstData(:,1:l_eeg) = AvgData_EEG;
TonalBurstData(:,l_eeg+1:l_eeg+l_lfp) = AvgData_LFP;
TonalBurstBinData(:,:,1:l_eeg) = EEG_bin;
TonalBurstBinData(:,:,l_eeg+1:l_eeg+l_lfp) = LFP_bin;
end