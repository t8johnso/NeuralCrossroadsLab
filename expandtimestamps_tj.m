function timestamps = expandtimestamps_tj(Data, Data_ind, Data_ts, Data_ts_step)
%Creates the timestamp variable based upon the data collected from the
%electrophysiological system.
%Inputs:
%Data - Nx1 Data array where N is the number of samples of the recording.
%Data_ind - Mx1 array that represents the indicies at which the recording
%was started or paused, usually exported with the data from plexon systems.
%Data_ts - Mx1 array that represents the corresponding time points since
%system acquisition began that data recording started or was paused.
%Data_ts_step - Scalar representing the sampling rate of the system (in ms)

%Outputs:
%timestamps - Nx1 array of the corresponding timestamps for the Data array.
    timestamps = ones(size(Data));
    if numel(Data_ind) > 1
        for i = 1:numel(Data_ind)-1
            nonmodtimearray = 0:Data_ind(i+1)-Data_ind(i)-1;
            timestamps(Data_ind(i):Data_ind(i+1)-1) = Data_ts_step*nonmodtimearray+Data_ts(i);
        end
    end
    nonmodtimearray = 0:numel(Data)-Data_ind(end);
    timestamps(Data_ind(end):numel(timestamps)) = Data_ts_step*nonmodtimearray+Data_ts(end);  
end