function itpc = intertrialphasecoupling(PhaseData)
%Calculates the inter-trial phase coupling to determine whether phase 
%(of a particular oscillation) is being preserved between trials.
%Inputs:
%PhaseData - NxMxO data array where N is the number of trials being
%evaluated, M is the number of timepoints within each trial, and O is the
%number of channels being evaluated. In the case of LFP for our study, this
%value is normally 16, for EEG, this value is either 1 (making PhaseData a
%NxM matrix) or 2.

%Outputs:
%itpc - Inter-trial phase coupling values as a MxO array or matrix
%(depending on the size of O). Each value representing the phase coupling
%across trials for each time stamp within trials.

s = size(PhaseData);
if numel(s) == 3
    for i = 1:s(3)
        for j = 1:s(2)
            itpc(j,i) = abs(mean(exp(1i*PhaseData(:,j,i))));
        end
    end
elseif numel(s) == 2
    for j = 1:s(2)
        itpc(j) = abs(mean(exp(1i*PhaseData(:,j))));
    end
end
end