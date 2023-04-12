function CIF = AmpPhase_CIF(amp, phase, lambda, savename)
%Description: Function designed to determine conditional intensity 
%relationships between amplitude and phase within periodic data sets (e.g. local 
%field potential phase). Additionally, provides figures for optimal 
%regressor selection (penalized negative log likelihood versus number of 
%regressors used), Q-Q plot for determining model uncertainty, and histogram 
%showing a rough indication of when spikes occurred with relation to phase. 
%Please note that use of this function requires ksdiscrete.m available at: 
%http://www.neurostat.mit.edu/software/time_rescaling

%Inputs:
%amp - Nx1 array consisting of amplitude values in volts where N is the 
%length of the entire dataset.
%phase - Nx1 array consisting of phase values in radians where N is the 
%length of the entire dataset.
%lambda - The number of L1 penalty terms used in the original ordering of 
%regressors, with 100 being standard. Increasing this value may increase 
%the accuracy of the ultimate conditional intensity function, but will 
%increase the computational time of the model.
%Savename - A string array that will ultimately serve as the name of saved 
%workspace.

%Outputs:
%CIF - 50x1 array representation describing how the conditional intensity 
%of a neuron’s spiking changes across a single phase cycle (-pi to pi).
%Figures -
%Figure 1: The normalized negative log-likelihood (NNLL) and the penalized 
%normalized negative log-likelihood (PNNLL) versus the number of regressors 
%used to model the data.
%Figure 2: A quantile-quantile plot describing the goodness of fit of the 
%model for the data being characterized.
%Figure 3: The conditional intensity function (probability of spiking given
%the circular dataset)versus the circular dataset.
%Figure 4: Histogram of the number of spikes occurring within the dataset at each 
%particular phase.
%Figure 5: Subplot containing Figures 1-4.
%Savefile - Workspace of the results of the workflow. Saved as a .mat file 
%with the string used for Savename. Workspace is saved automatically to 
%active folder. Depending on the size of the dataset being analyzed, some 
%non-essential components of the workflow may not be saved.

savefilename = [savename '.mat'];

%Inverse variance and mean values set here. Able to modify these values to 
kappa = linspace(.01,30,20);
deltamu = linspace(0,2*pi,20);
deltamu = deltamu(1:end-1);
regressor_temp = zeros(numel(amp), length(deltamu)*length(kappa));
param_hold = zeros(length(deltamu)*length(kappa), 2);
for j = 1:length(deltamu)
    for l = 1:length(kappa)
        %In these loops, the parameters are being stored into an array for
        %easy recall later and the regressor matrix is being created using
        %the mean and variance pairs previously defined.
        param_hold((j-1)*length(kappa)+l,:) = [deltamu(j),kappa(l)];
        R_temp = exp(kappa(l).*cos(phase+deltamu(j)))./(2.*pi.*besseli(0,kappa(l))) - min(exp(kappa(l).*cos(phase+deltamu(j)))./(2.*pi.*besseli(0,kappa(l))));
        R = R_temp./max(R_temp);
        regressor_temp(:,(j-1)*length(kappa)+l) = R;
    end
end

%L1 Regularized GLM step. This is used to determine model ordering.
[B,~] = lassoglm(regressor_temp, amp', 'normal', 'Link', 'identity','NumLambda',lambda, 'Options',statset('UseParallel',true),'MaxIter',1000);

%Temporarily saves the L1 weights due to how computationally intensive this
%step is.
save('B_temp.mat','B')

%Flips the weight ordering in order to have increasing model complexity.
NNZ_temp = zeros(1,numel(B(1,:)));
C = fliplr(B);

%Clears the workspace to save RAM space.
clear B kappa deltamu

%Following steps are used to determine the number of basis functions being
%used for each gamma value and stores the mean variance pairs as well as 
%the corresponding regressors for each of these gamma values.
for i = 1:numel(C(1,:))
    NNZ_temp(i) = nnz(C(:,i));
end

[uniqueNNZ, ia] = sort(NNZ_temp, 'ascend');

sortorder = ia(find(uniqueNNZ ~= 0));
regressor = {};

for j = 1:numel(sortorder)
    regressor{j} = regressor_temp(:,find(C(:,sortorder(j))~=0));
end

%Using the ordering previously determined by the L1 regularization,
%weights are now created for each increasing model order as well as the
%NNLL and PNNLL of each model order.
x = {};
for i = 1:length(sortorder)
    VM = fitglm(regressor{i},amp');
    x{i} = VM.Coefficients.Estimate;
    NNLL(i) = (log(2*pi*std(amp)^2)/2) + sum(((amp' - [ones(length(regressor_temp),1), regressor{i}]*x{i}).^2)/(2*numel(amp)*std(amp)^2));
    [~, d] = size(regressor{i});
    PNNLL(i) = NNLL(i) + d/length(amp);
end

%Removes the 0 model order and plots the lower end of model orders to show
%the PNNLL-model order relationship.
NNZ = uniqueNNZ(2:end);
[Q, UIA] = unique(NNZ);

figure(1)
plot(Q, NNLL(UIA),'k')
hold on
plot(Q,PNNLL(UIA),'r')
xlabel('Model Order')
ylabel('Normalized Negative Log-Likelihood')

%Find the optimal model order by looking for the global minimum. The
%regressors used for this model order are the evaluated to produce the QQ
%plot and the CIF plot
optorder_temp = find(PNNLL == min(PNNLL));

optorder = optorder_temp(1);
modelorder = uniqueNNZ(optorder+1);
morder = modelorder(1);
params_hold = param_hold(find(C(:,sortorder(optorder))~=0),:);

NNLL_plot = NNLL;
NNLL = NNLL(optorder);
VM = fitglm(regressor{optorder},amp');

CIF_temp = [ones(length(regressor_temp),1), regressor{optorder}]*x{optorder};

[~,rstsort,xks,cb,~] = ksdiscrete(CIF_temp, amp','spiketrain');
figure(2)
plot(xks,rstsort,'k-');
hold on;
plot(xks,xks+cb,'k--',xks,xks-cb,'k--');
axis([0,1,0,1])

%Conditional intensity function produced from a single cycle just to
%represent the relationship for spiking-phase in a simpler form.
phasetemp = linspace(-pi,pi,50);
params_length = size(params_hold,1);
regressorprob = ones(50,params_length+1);
for l = 1:params_length
    R2_temp = exp(params_hold(l,2).*cos(phasetemp+params_hold(l,1)))./(2.*pi.*besseli(0,params_hold(l,2))) - min(exp(params_hold(l,2).*cos(phase+params_hold(l,1)))./(2.*pi.*besseli(0,params_hold(l,2))));
    R2 = R2_temp./max(R2_temp);
    regressorprob(:,l+1) = exp(params_hold(l,2).*cos(phasetemp+params_hold(l,1)))./(2.*pi.*besseli(0,params_hold(l,2)));
end
CIF = regressorprob*x{optorder};

figure(3)
plot(phasetemp,CIF','r')
xlabel('Phase (rad)')
ylabel('P(Amp|Phase)')
xlim([-pi pi])

figure(4)
histogram(phase(amp==1),50)
xlim([-pi pi])
ylabel('# of Spikes')
xlabel('Phase (rad)')

figure(5)
subplot(2,2,1)
plot(Q, NNLL_plot(UIA),'k')
hold on
plot(Q,PNNLL(UIA),'r')
xlabel('Model Order')
ylabel('Normalized Negative Log-Likelihood')

subplot(2,2,2)
plot(xks,rstsort,'k-');
hold on;
plot(xks,xks+cb,'k--',xks,xks-cb,'k--');
axis([0,1,0,1])

subplot(2,2,3)
plot(phasetemp,CIF','r')
xlabel('Phase (rad)')
ylabel('P(Amp|Phase)')
xlim([-pi pi])

subplot(2,2,4)
histogram(phase(amp==1),50)
xlim([-pi pi])
ylabel('# of Spikes')
xlabel('Phase (rad)')

save(savefilename)
end