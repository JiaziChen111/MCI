function [samp] = spm_mci_ahb_single (MCI,mcmc)
% A Single Trajectory of Annealed Hierarchical Bayes
% FORMAT [samp] = spm_mci_ahb_single (MCI,mcmc)
%
% MCI       Data structure containing fields:
%
% .M{n}     Model for nth of N replications (e.g. subjects)
% .U{n}     Inputs for nth replication
% .Y{n}     Data for nth replication
% .mlm      Multivariate linear model describing group effects
%           See spm_mci_mlm.m
%
% mcmc      Optimisation parameters  
%
% samp      Estimates from single trajectory
% .traj     traj(p,j,n) is value of parameter p for subject n 
%           at temperature j (only returned if mcmc.rec_traj=1)
%
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

M=MCI.M; U=MCI.U; Y=MCI.Y; mlm=MCI.mlm;
beta=mcmc.beta;
nprop=mcmc.nprop;
prop=mcmc.prop;
scale=mcmc.scale;

J=length(beta);
N=length(M);

infer_Lambda=mlm.infer_Lambda;
if ~infer_Lambda
    final_Lambda=mlm.Lambda;
end

% Initialisation
for n=1:N,
    if isstruct(M{n}.pC)
        pC = full(diag(spm_vec(M{n}.pC)));
    else
        pC = M{n}.pC;
    end
    % prior cov in subspace
    pC=full(M{n}.V'*pC*M{n}.V);
    Np=size(pC,1);
    
    % Sample parameters from prior (in subspace)
    x(:,1,n) = spm_mci_proposal (M{n},U{n},zeros(Np,1),pC);
    
    if isfield(M{n},'Ce')
        % Sample obs noise covar from prior
        Ny=size(Y{n},2); noise.D0=eye(Ny); noise.c0=Ny;
        Lprec=spm_wishrnd(noise.D0,noise.c0);
        M{n}.Ce=inv(Lprec);
    end
    
    [L(1,n),L2(n)] = spm_mci_joint (x(:,1,n),M{n},U{n},Y{n});
    Lsum(n) = (beta(2)-beta(1))*L2(n);
    acc(1,n) = 1;
end

lgv_mcmc.maxits=nprop;
if isfield(MCI.M{1},'Ce')
    lgv_mcmc.update_obs_noise=mcmc.update_obs_noise;
    lgv_mcmc.update_obs_step=mcmc.update_obs_step;
end
    
for j=2:J,
    
    % Subject Models
    for n=1:N,
        xs(:,n)=x(:,j-1,n); Ls(n)=L(j-1,n); acc(j,n)=0;
        
        % Assume, for now, that all subject effects are random effects
        if j>2
            % Set prior mean and covariance from 2nd level model
            M{n}.pE=Rhat(:,n);
            M{n}.vpE=M{n}.pE;
            M{n}.pC=pC;
        end
        
        % Generate sample at next temperature using LMC
        M{n}.beta=beta(j);
        if j==2
            % Need param in original space when calling spm_mci_lgv.m
            xs(:,n)=M{n}.V*xs(:,n)+M{n}.vpE;
        end
        lgv_mcmc.init=xs(:,n);
        [Mtmp,stats] = spm_mci_lgv (lgv_mcmc,M{n},U{n},Y{n});
        
        % Update obs noise estimate for Gauss noise models
        if isfield(M{n},'Ce')
            M{n}.Ce=Mtmp.Ce;
        end
        x(:,j,n) = stats.P(:,end);
        L(j,n) = -stats.E(:,end);
        L2 = stats.L2(end);
        acc (j,n) = any(stats.acc(2:end));
        if j < J
            Lsum(n) = Lsum(n) + (beta(j+1)-beta(j))*L2;
        end
    end
    
    if j < 0.75*J,
        mlm.infer_Lambda=0;
        mlm.Lambda=mlm.Psi0;
    else
        mlm.infer_Lambda=infer_Lambda;
        if ~infer_Lambda
            mlm.Lambda=final_Lambda;
        end
    end
        
    % Group Model
    mlm_Y=squeeze(x(:,j,:));
    y=mlm_Y(:);
    [mlm,post] = spm_mci_mlm (mlm,y,1);
    if mcmc.rec_traj
        m(:,j)=post.m;
    end
    
    % Predictions of subject effects
    rhat=mlm.U*post.m;
    Rhat=reshape(rhat,mlm.K,mlm.N);
    pC=inv(post.Lambda);
    
end

samp.logw=sum(Lsum);
samp.w=squeeze(x(:,J,:));
samp.m=post.m;

if mcmc.rec_traj
    samp.traj=x;
    samp.trajm=m;
end
samp.Lambda=post.Lambda;
samp.E=-L(J,:);
if isfield(M{1},'Ce')
    for n=1:N,
        samp.Ce(:,:,n)=M{n}.Ce;
    end
else
    samp.Ce=[];
end
samp.acc=acc;