function [post] = spm_mci_ais (mcmc,M,U,Y,vl)
% Annealed Importance Sampling
% FORMAT [post] = spm_mci_ais (mcmc,M,U,Y,vl)
% 
% mcmc      Optimisation parameters  eg.
%
% .J        number of temperatures
% .anneal   annealing schedule:
%           'sigmoid', 'linear', 'nonlinear', 'log' or 'power'
% .prop     type of proposal: 'lmc' (default) or 'mh' 
% .nprop    number of proposals at each temperature
% .maxits   number of independent samples to produce
%
% M         Model structure 
% U         Input structure
% Y         Data 
% vl        Variational Laplace solution
%               .Ep                 Posterior Mean
%               .Cp                 Posterior Covariance
%           If this field is specified then AIS starts sampling
%           from the VL posterior. Otherwise from the model prior.
%
% The function returns data structure 'post' with fields
%
% .P                P(:,i) is ith posterior sample
% .logev            approximation to log evidence
% .logev_se         standard error thereof
% .logev_lower      5th percentile thereof
% .logev_upper      95th percentile thereof
% .logev_resample   resampled log evidences
% .traj             individual trajectories
% .acc              acceptances
% .logw             log of (unnormalised) importance weights
% .q                normalised importance weights
% .E                E(:,i) is energy (negative log joint) of ith smaple
% .beta             set of inverse temperatures
%
% R Neal (2001) Annealed Importance Sampling. Statistics and
% Computing, 11, 125-139.
%
% This implementation uses the Matlab Parallel Computing toolbox
% (see use of parfor instead of for below).
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_ais.m 6548 2015-09-11 12:39:47Z will $
            
try J=mcmc.J; catch J=512; end
try anneal=mcmc.anneal; catch anneal='power'; end
try prop=mcmc.prop; catch mcmc.prop='lmc'; end
try maxits=mcmc.maxits; catch maxits=32; end
try nprop=mcmc.nprop; catch nprop=1; end
try verbose=mcmc.verbose; catch verbose=1; end
try mcmc.rec_traj; catch mcmc.rec_traj=0; end
if nargin < 5
    vl=[];
else
    vl=vl;
end

if isfield(M,'Ce')
    try
        mcmc.update_obs_noise=mcmc.update_obs_noise;
    catch
        % Update observation noise by default
        mcmc.update_obs_noise=1;
        mcmc.update_obs_step=0;
    end
else
    mcmc.update_obs_noise=0;
    mcmc.update_obs_step=0;
end

nprop=nprop+1;

M = spm_mci_minit (M);
V  = M.V;
Np = size(V,1);
Nsub = size(V,2);

if verbose
    disp(sprintf('Using %s annealing schedule with %d temperatures',anneal,J));
end

beta = spm_mci_schedule (anneal,J);

mcmc.scale=fliplr(10.^-[0:1:nprop-1]);
mcmc.beta=beta;
mcmc.nprop=nprop;

if ~isempty(vl)
    vl.Cp=full(vl.Cp);
    vl.Lambdap=inv(vl.Cp);
    vl.const=-0.5*Np*log(2*pi)-0.5*spm_logdet(vl.Cp);
    % Mean and cov in subspace:
    vl.mr=M.V'*vl.Ep-M.vpE
    vl.Cr=M.V'*vl.Cp*M.V;
end

parfor i=1:maxits,
%for i=1:maxits,
    if verbose
        disp(sprintf('Acquiring %d out of %d IID samples',i,maxits));
    end
    
    if ~isempty(vl)
        [P(:,i),E(i),logw(i),acc(i,:),traj(i,:,:)] = spm_mci_ais_single_vl (mcmc,M,U,Y,vl);
    else
        ais_samp = spm_mci_ais_single (mcmc,M,U,Y);
        P(:,i) = ais_samp.P;
        E(i) = ais_samp.E;
        logw(i) = ais_samp.logw;
        acc(i,:) = ais_samp.acc;
        traj(i,:,:) = ais_samp.traj;
        if isfield(M,'Ce')
            Ce(:,:,i) = ais_samp.Ce;
        end
        tels(i) = ais_samp.els;
    end
end

% Model evidence and error bars
[logev,q] = spm_mci_ais_evidence (logw);
boot = spm_mci_ais_bootstrap (logw);

% Get Posterior Mean using importance weights
post.wEp=P*q';

post.P=P;
if isfield(M,'Ce')
    post.Ce=Ce;
end
post.logev=logev;
post.logev_se=boot.logev_se;
post.logev_lower=boot.logev_low;
post.logev_upper=boot.logev_high;
post.acc=acc;
post.traj=traj;
post.q=q;
post.E=E;
post.beta=beta;
post.ind=[1:maxits];
post.logw=logw;
post.tels=tels;

