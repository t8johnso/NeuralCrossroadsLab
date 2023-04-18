function sig_P = P_Correction(p_val,sig_thresh,Corr_flag)
%Function used to make p-value corrections based upon different correction
%methods.
%Inputs: 
%p-val - MxN matrix containing p-values.
%sig_thresh - scalar value of the significance threshold to be corrected.
%Corr_flag - Keyword string either "BH" for Benjamini-Hochberg correction,
%"BF" for Bonferonni correction, or "NO" for no correction.

%Outputs:
%sig_P - Binary MxN matrix where 1s indicate indicies that are significant
%based upon the desired correction method.

try
    p_val = cell2mat(p_val);
end
[M, N] = size(p_val);
O = numel(p_val);
sig_P = zeros(M,N);
if Corr_flag == 'BH'
    A = reshape(p_val, numel(p_val), 1);
    B = sort(A);
    for i = 1:numel(B)
        C(i) = B(i).*(O)./i.*0.25;
    end
    C = C';
    trig = find(B>C);
    sig_P(find(p_val<B(trig(1)))) = 1;
elseif Corr_flag == 'BF'
    thresh_corr = 0.05./O;
    sig_P(find(p_val<thresh_corr)) = 1;
elseif Corr_flag == 'NO'
    sig_P(find(p_val<sig_thresh)) = 1;
end
end
    