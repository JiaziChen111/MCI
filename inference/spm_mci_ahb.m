function [post,MCI] = spm_mci_ahb (MCI,mcmc)
% Annealed Hierarchical Bayes
% FORMAT [post,MCI] = spm_mci_ahb (MCI,mcmc)
%
% MCI               Data structure containing fields:
%
% .M{n}             Model for nth of N replications (e.g. subjects)
% .U{n}             Inputs for nth replication
% .Y{n}             Data for nth replication
% .mlm              Multivariate linear model describing group effects
%                   .infer_Lambda=0 to use fixed Lambda
%                   See spm_mci_mlm.m
%
% mcmc      Optimisation parameters  eg.
%
% .J        number of temperatures
% .anneal   annealing schedule:
%           'sigmoid', 'linear', 'nonlinear', 'log', 'power', 'frozen'
% .prop     type of proposal: 'lmc' or 'mh' (default)
% .nprop    number of proposals at each temperature
% .maxits   number of independent samples to produce
% .rec_traj set to 1 to record trajectories (default is 0)
%
% The output fields are: 
%
% POSTERIOR SAMPLES (post)
%
% .m               [(K x B) x Nsamples] group effects, m
% .w               [K x N x Nsamples] subject effects, w
% .Ce               [Ny x Ny x N x Nsamples] Obs noise covariance samples
% .acc             acceptance values (are proposals accepted ?)
% .E               Energy (neg log joint prob)
% .Lambda          Between-subject precision matrix (posterior samples)
%
% where K is the number of first level random effects (eg. within subject)
% and B is the number of second level effects (e.g. between subject)
%
% POSTERIOR MOMENTS (MCI)
%
% .mlm.Ep          Posterior mean at group level (from post.m)
% .mlm.Cp          Posterior covariance at group level (from post.m)
% .mlm.Lambda      Posterior mean precision at group level (from
%                  post.Lambda)
% .M{n}.Ep         Posterior mean for subject n (from post.w)
% .M{n}.Cp         Posterior covariance for subject n (from post.w)
%
% The following fields are only returned if mcmc.rec_traj=1
% .traj     traj(p,j,n,i) is value of parameter p for subject n 
%           at temperature j for ith trajectory.
% .trajm    trajm(p,j,i) is value of parameter p for group 
%           at temperature j for ith trajectory.
%
% BLOCKED GIBBS SAMPLING
%
% The annealing aspect of this algorithm can be bypassed by using a 
% 'frozen' annealing trajectory (inverse temperature fixed at 1 for all j).
% For each trajectory (chain), i, mcmc.J samples will be returned in the
% above .traj and .trajm variables (you need to set mcmc.rec_traj=1).
% Samples from a single chain will be returned by setting mcmc.maxits=1.
%
% W. Penny, M. Klein-Flugge and G Ziegler. Annealed Hierarchical Bayes,
% Submitted, 2017
% 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

try J=mcmc.J; catch J=512; end
try anneal=mcmc.anneal; catch anneal='power'; end
try prop=mcmc.prop; catch prop='mh'; end
try maxits=mcmc.maxits; catch maxits=32; end
try nprop=mcmc.nprop; catch nprop=1; end
try verbose=mcmc.verbose; catch verbose=1; end
try vl=vl; catch vl=[]; end;
try mcmc.rec_traj; catch mcmc.rec_traj=0; end

if isfield(MCI.M{1},'Ce')
    try
        mcmc.update_obs_noise=mcmc.update_obs_noise;
    catch
        % Update observation noise by default
        mcmc.update_obs_noise=1;
        mcmc.update_obs_step=0;
    end
end

nprop=nprop+1;

N=length(MCI.M);
for n=1:N,
    MCI.M{n} = spm_mci_minit (MCI.M{n});
end

if verbose
    disp(sprintf('Using %s annealing schedule with %d temperatures',anneal,J));
end

beta = spm_mci_schedule (anneal,J);

if isfield(MCI.M{1},'Ce')
    try
        mcmc.update_obs_noise=mcmc.update_obs_noise;
    catch
        % Update observation noise by default
        mcmc.update_obs_noise=1;
        mcmc.update_obs_step=0;
    end
end

mcmc.scale=fliplr(10.^-[0:1:nprop-1]);
mcmc.beta=beta;
mcmc.nprop=nprop;

% Set up second-level model
mlm = MCI.mlm;
mlm.B = size(mlm.X,2);

disp('spm_mci_ahb: assuming that all first-level effects are random effects');

mlm.K = length(MCI.M{1}.pE);
mlm.N = length(MCI.M);
mlm.U = kron(mlm.X,eye(mlm.K));

try mlm.Psi0=mlm.Psi0; catch mlm.Psi0 = 0.1*eye(mlm.K*mlm.B); end
try mlm.m0 = mlm.m0; catch mlm.m0 = zeros(mlm.K*mlm.B,1); end
try mlm.a0=mlm.a0; catch mlm.a0 = 1; end
try mlm.B0=mlm.B0; catch mlm.B0 = 0.1*eye(mlm.K); end

mlm.m = spm_normrnd(mlm.m0,inv(mlm.Psi0),1);
if mlm.infer_Lambda
    mlm.Lambda = spm_wishrnd(mlm.B0,mlm.a0);
end
MCI.mlm = mlm;

% Loop over multiple trajectories
parfor i=1:maxits,
%for i=1:maxits,
    if verbose
        disp(sprintf('Acquiring %d out of %d IID samples',i,maxits));
    end
    
    % A single trajectory
    ais_samp = spm_mci_ahb_single (MCI,mcmc);
    
    w(:,:,i) = ais_samp.w;
    E(:,i) = ais_samp.E(:);
    logw(i) = ais_samp.logw;
    acc(:,:,i) = ais_samp.acc;
    Ce(:,:,:,i) = ais_samp.Ce;
    m(:,i) = ais_samp.m;
    Lambda(:,:,i)=ais_samp.Lambda;
    if mcmc.rec_traj
        traj(:,:,:,i)=ais_samp.traj;
        trajm(:,:,i)=ais_samp.trajm;
    end
end

% Model evidence
[logev,q] = spm_mci_ais_evidence (logw);

% Bootstrapped confidence intervals on evidence
boot = spm_mci_ais_bootstrap (logw);

% Effective number of samples (page 129, Neal 2001)
vq=std(q)^2;
Nq=length(q);
post.Neff=Nq/(1+Nq*vq);
    
post.w=w;
post.Ce=Ce;
post.logev=logev;
post.logev_se=boot.logev_se;
post.logev_lower=boot.logev_low;
post.logev_upper=boot.logev_high;
post.acc=acc;
post.q=q;
post.E=E;
post.beta=beta;
post.ind=[1:maxits];
post.logw=logw;
post.m=m;
post.Lambda=Lambda;

if mcmc.rec_traj,
    post.traj=traj;
    post.trajm=trajm;
end

% Posterior mean and covariances for each subject
for n=1:N,
    w=squeeze(post.w(:,n,:))';
    MCI.M{n}.Ep = mean(w)';
    MCI.M{n}.Cp = cov(w);
end

% Posterior mean and covariance at group level
MCI.mlm.Ep=mean(post.m')';
MCI.mlm.Cp=cov(post.m');
MCI.mlm.Lambda=mean(post.Lambda,3);
MCI.mcmc=mcmc;

end

%-------------------------------------------------------

function [logev,q] = ais_evidence (logw)

r=max(logw);
u=exp(logw-r);
S=mean(u);
logev=log(S)+r;
q=u/sum(u);

end