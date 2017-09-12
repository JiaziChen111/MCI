
clear all
close all

disp('Annealed Hierarchical Bayes on Logistic regression models');
    
% True parameters
trueP=[1,1]';

% Number of subjects
N=16;

% Number of data points per subject
T=20;

% Ratio of between-subject precision to prior precision
true_lambda=16;

% ----------- Generate Data -----------------------
[M0,U0] = mci_logistic_struct ('dct2',T);
for n=1:N,
    
    M{n}=M0; U{n}=U0;
    
    % Synthetic data, variable pars for all subjects
    pE=trueP;
    Np=length(pE);
    pC=M{n}.pC/true_lambda;
    MCI.trueP(:,n)=spm_normrnd(pE,pC,1);
    [g,y] = mci_logistic_gen(MCI.trueP(:,n),M0,U0);
    Y{n}.y=y;
end

% MCI
disp('MCI ...');
MCI.verbose=0;
MCI.M=M; MCI.U=U; MCI.Y=Y;

% -------------  Fit Model -------------------------

mcmc.inference='ais';
mcmc.anneal='power';
mcmc.prop='lmc';
mcmc.nprop=1;
mcmc.J=512;
mcmc.maxits=3;

% Second-Level Model
mlm.X = ones(N,1);

% Prior mean and precision of second level paramaters (mean over subjects)
mlm.m0 = M{1}.pE;
mlm.Psi0 = inv(M{1}.pC);

mlm.infer_Lambda=1;
if mlm.infer_Lambda
    mlm.a0 = 1;
    Np = size(M{1}.pC,1);
    mlm.B0 = eye(Np);
else
    % Between subject variance is 1/svr second level prior variance
    % (as in page 8 of Friston et al. Bayesian model reduction ... Neuroimage, 2015
    svr = true_lambda;
    mlm.Lambda=svr*mlm.Psi0;
end
MCI.mlm=mlm;

tic;
[post,MCI] = spm_mci_ahb (MCI,mcmc);
els=toc
post.els=els;
