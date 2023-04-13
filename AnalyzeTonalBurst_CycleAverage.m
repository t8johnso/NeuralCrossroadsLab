function [ASSRData, ASSRbins, removedbins, removedbin_ind] = AnalyzeTonalBurst_CycleAverage(Data, events, TBFreq, shift, edges, Cyclenum, Fs, timestamps, noise_timestamps)
%Description: Bins neural data according to onset events of stimulus.
%Seperates out noise epochs while identifying the segments that were
%removed. Requires isolateblockdata function.

%Inputs:
%Data - Nx1 array consisting of voltage values in volts where N is the 
%length of the entire dataset.
%events - 2x1 array consisting of a starting point and stopping point of
%what segment of data to trim around the stimulus onset and offset.
%TBFreq - Scalar value of the frequency of the tonalburst used as an
%auditory stimulus in Hz.
%shift - Scalar value that shifts the data as is desired in ms.
%edges - Scalar value that includes additional data within the desired
%numbed of cycles in ms.
%Cyclenum - Scalar value describing the desired number of cycles to be
%segmented together.
%Fs - Scalar value sampling frequency in Hz.
%timestamps - Nx1 array consisting of timing value that pairs with the data
%in ms.
%noise_timestamps - Variably sized array containing time stamp values to be
%excluded from data averaging in ms.

%Outputs:
%ASSRData - Array of voltage values in volts that have been averaged across
%cycle segments.
%ASSRbins - Matrix of voltage values in volts segmented by number of
%cycles.
%removedbins - Matrix of voltage values in volts segmented by number of
%cycles for the noisy data segments.
%removedbin_ind - Array containing indicies of segments that were excluded
%from ASSR analysis.

[segData, segTime] = isolateblockdata(Data, timestamps, events(1)+shift, events(2)+shift, timestamps);
segData = segData';
segTime = segTime';
cycleperiod = round(Fs/TBFreq);
binperiod = cycleperiod*Cyclenum;
s = size(segData);
numdatabins = s(2)/binperiod;
if round(numdatabins) ~= numdatabins
    numdatabins = floor(numdatabins);
    trimData = segData(:,1+edges(1):binperiod*numdatabins);
    trimTime = segTime(:,1+edges(1):binperiod*numdatabins);
else
    trimData = segData;
    trimTime = segTime;
end
Databin =[];
removedbins = [];
removedbin_ind = [];
for i = 1:numdatabins-3
    bin_TS = trimTime(1+(i-1)*binperiod-edges(1):1+i*binperiod+edges(2));
    [check, ~] = ismember(bin_TS, noise_timestamps);
    check2 = find(check ~= 0);
    if isempty(check2) == 1
        Databin(i,:) = trimData(:,1+(i-1)*binperiod-edges(1):1+i*binperiod+edges(2));
    else
        removedbins = [removedbins; trimData(:,1+(i-1)*binperiod-edges(1):1+i*binperiod+edges(2))];
        removedbin_ind = [removedbin_ind, i];
    end
end
ASSRbins = Databin(any(Databin,2),:);
ASSRData = mean(ASSRbins, 1);
end