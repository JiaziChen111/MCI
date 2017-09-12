
clear all
close all

% Number of data points per condition
N=100;

% Generate data from baseline model (yes or no)
baseline=0;

if baseline
    disp('True Model: Baseline');
else
    disp('True Model: Two conditions');
end
disp(' ');

% Two condition model
[M0,U0] = mci_ddm_struct('vconflict2');
M0.L = 'mci_ddm_like_mc';
M0.IS = 'mci_ddm_gen_mc';

% Baseline model
[M1,U1] = mci_ddm_struct('unbiased');
M1.L = 'mci_ddm_like_sc';
M1.IS = 'mci_ddm_gen_sc';

% Condition 1
ddm_true{1}.a=2.1; 
ddm_true{1}.b=0.45;
ddm_true{1}.v=0.3; 
P{1} = mci_ddm_to_mci(ddm_true{1},M0,U0);

% Condition 2
ddm_true{2}=ddm_true{1};
ddm_true{2}.v=0.6; 
P{2} = mci_ddm_to_mci(ddm_true{2},M0,U0);
      
Ptrue = [P{1}(1), P{2}]';

% Data 
U0.cond{1}=[1:N]';
U0.cond{2}=[N+1:2*N]'; 

for c=1:2,
    % Generate data from PDFs
    if baseline
        [y{c},dist]=mci_ddm_gen(P{1},M0,U0,N);
    else
        [y{c},dist]=mci_ddm_gen(P{c},M0,U0,N);
    end
    
    x=y{c}(:,1);
    t=y{c}(:,2);
    tmin(c)=min(t);
    
    disp(sprintf('Condition %d',c));
    disp(sprintf('Theoretical percent correct = %1.2f',dist.pc));
    disp(sprintf('Empirical percent correct = %1.2f', mean(x)));
end

% Set maximum b value to less than minimum RT
M0.b_max=min(tmin)-0.01;
M1.b_max=M0.b_max;    

Y=[y{1};y{2}];

mcmc=[];

disp(' ');
disp('Using Laplace method for estimation and inference');
disp(' ');

tic;
lp0 = spm_mci_laplace (mcmc,M0,U0,Y);
els0=toc;

tic;
lp1 = spm_mci_laplace (mcmc,M1,U1,Y);
els1=toc;

disp(' ');
disp('Bayes factor in favour of two condition model');
disp(sprintf('is %1.2f',lp0.logev-lp1.logev));
