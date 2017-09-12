function [samp] = spm_mci_ais_single (mcmc,M,U,Y)
% Produce a single independent sample using AIS
% FORMAT [samp] = spm_mci_ais_single (mcmc,M,U,Y)
%
% mcmc      Sampling settings
% M         Model structure
% U         Input structure
% Y         Data
%
% samp      Data structure for returned sample with fields
% .P         [Np x 1] parameter sample
% .Ce        observation noise covariance sample 
% .E         Negative log joint
% .logw      Contribution to model evidence
% .acc       acc(j) is acceptance rate at temperature j
% .traj      traj(p,j) is value of parameter p at temperature j
%            (only set if mcmc.rec_traj=1)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_ais_single.m 6548 2015-09-11 12:39:47Z will $

tstart=tic;

beta=mcmc.beta;
nprop=mcmc.nprop;
prop=mcmc.prop;
scale=mcmc.scale;

J=length(beta);

if isstruct(M.pC)
    pC = full(diag(spm_vec(M.pC)));
else
    pC = M.pC;
end
% prior cov in subspace
pC=full(M.V'*pC*M.V);

Np=size(pC,1);

% Sample parameters from prior (in subspace)
%x(:,1) = spm_mci_proposal (M,U,zeros(Np,1),pC);
% NOTE pE not nec. zero !!!
x1 = spm_mci_proposal (M,U,zeros(Np,1),pC);

if isfield(M,'Ce') & mcmc.update_obs_noise
    % Sample obs noise covar from prior
    Ny=size(Y,2); noise.D0=eye(Ny); noise.c0=Ny; 
    Lprec=spm_wishrnd(noise.D0,noise.c0);
    M.Ce=inv(Lprec);
end

%[L(1),L2] = spm_mci_joint (x(:,1),M,U,Y);
[L(1),L2] = spm_mci_joint (x1,M,U,Y);

% Need sample in original space when calling eg. spm_mci_lgv
x(:,1)=M.V*x1+M.vpE;

Lsum = (beta(2)-beta(1))*L2;
acc(1) = 1;
for j=2:J,
    xs=x(:,j-1); Ls=L(:,j-1);acc(j)=0;
    switch prop,
        case 'mh',
            % Generate sample at next temperature by composing
            % Metropolis moves at different scales
            for s=1:nprop,
                dx=spm_normrnd(zeros(Np,1),scale(s)^2*pC,1);
                xcand=xs+dx;
                [Lcand,L2cand] = spm_mci_joint (xcand,M,U,Y,beta(j));
                dL = Lcand-Ls;
                r = exp(dL);
                alpha = min(1,r);
                test_prob = rand(1);
                if alpha > test_prob
                    % Accept
                    xs = xcand;
                    Ls = Lcand;
                    L2 = L2cand;
                    acc (j) = 1;
                end
            end
            
        case 'lmc',
            % Generate sample at next temperature using
            % Langevin Monte Carlo
            M.beta=beta(j);
            %if j==2
                % Need param in original space when 
                % calling spm_mci_lgv.m
                %xs=M.V*xs+M.vpE;
            %end
            lgv_mcmc.init=xs;
            lgv_mcmc.maxits=nprop;
            if isfield(mcmc,'update_obs_noise')
                lgv_mcmc.update_obs_noise=mcmc.update_obs_noise;
            end
            if isfield(mcmc,'h')
                lgv_mcmc.h=mcmc.h;
            end
            
            [Mtmp,stats] = spm_mci_lgv (lgv_mcmc,M,U,Y);
            
            % Update obs noise estimate for Gauss noise models
            if isfield(M,'Ce')
                M.Ce=Mtmp.Ce;
            end
            
            xs = stats.P(:,end);
            Ls = -stats.E(:,end);
            L2 = stats.L2(end);
            acc (j) = any(stats.acc(2:end));
            
        otherwise
            disp('Unknown proposal type in spm_mci_ais_single.m');
    end
    x(:,j)=xs;
    L(j)=Ls;
    if j < J
        Lsum = Lsum + (beta(j+1)-beta(j))*L2;
    end
    
    % Stop if this trajectory has taken too long
    te(j)=toc(tstart);
    telapsed=te(j);
    if telapsed > 200
        keyboard
    end
end
samp.logw=Lsum;
samp.P=x(:,J);
samp.E=-L(J);
if isfield(M,'Ce')
    samp.Ce=M.Ce;
else
    samp.Ce=[];
end

if mcmc.rec_traj
    %nj=size(x,2);
    %samp.traj=spm_vec(M.pE)*ones(1,nj)+M.V*x;
    samp.traj=x;
else
    samp.traj=[];
end

samp.els=telapsed;
samp.acc=acc;