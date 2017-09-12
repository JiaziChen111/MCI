function [post,M] = spm_mci_laplace (mcmc,M,U,Y)
% Laplace approximation with posterior mode from Matlab optimiser
% FORMAT [post,M] = spm_mci_laplace (mcmc,M,U,Y)
%
% mcmc  Optimisation parameters
%       .opt = 'fminsearch' (default) or 'fminunc'
%       from Matlab's optimisation toolbox
% M     Model structure
% U     Inputs
% Y     Data
%
% post  Structure with fields:
%
% .Ep       Posterior mean
% .Cp       Posterior covariance
% .logev    Log model evidence
% .L        Log Joint of .Ep
%
% M     Updated model structure
% 
% The default optimisation 'fminsearch' uses the Nelder-Mead Simplex method
% The 'fminunc' algorithm uses a Trust Region method if you supply the
% gradient but, if not, will use a Quasi-Newton method.
%

% Will Penny
% $Id$

try opt=mcmc.opt; catch opt='fminsearch'; end
try sample=mcmc.sample; catch sample=0; end

if isstruct(Y)
    Y=Y.y;
end

% Compute eigen-parameterisation
M = spm_mci_minit (M);
V  = M.V;

% Initial parameter values
if sample
    % Sample from prior
    if isstruct(M.pC)
        pC = full(diag(spm_vec(M.pC)));
    else
        pC = M.pC;
    end
    % prior cov in subspace
    pC=full(M.V'*pC*M.V);
    % prior mean in subspace
    pE=V'*M.vpE;
    xinit = spm_mci_proposal (M,U,pE,pC);
    xtypical = sqrt(diag(pC));
else
    % User specified or prior mean
    try init=mcmc.init; catch init=M.vpE; end
    
    % Initial param in eigenspace
    xinit = M.V'*(init-M.vpE);
end

options = optimset(opt);
options.GradObj='on';
options.TolFun=1e-6;
%options.Display='iter';

% Default number is 200*Np
Np=length(xinit);
options.MaxFunEvals=512*Np;

switch opt
    case 'fminunc';
        [x,E]=fminunc('spm_mci_energy',xinit,options,M,U,Y);
    case 'fminsearch';
        [x,E]=fminsearch('spm_mci_energy',xinit,options,M,U,Y);
    otherwise
        disp('Error: Unknown optimiser in spm_mci_laplace.m');
        return
end

% Posterior covariance
M.hess_outer=0;
[j,iCpY] = spm_mci_joint_grad (x,M,U,Y);
Cp = inv(iCpY+M.ipC);

% Project parameters back from eigenspace into original space
post.Ep = M.vpE+V*x;
post.Cp = V*Cp*V';
post.E = E;
post.L = -E;

% Model Evidence 
% logev = log p(Y,w) + 0.5 Np log 2pi + 0.5 log |Cp|
post.logev = -E + 0.5*M.Np*log(2*pi)+0.5*spm_logdet(Cp);

% Model Evidence for Gaussian Likelihoods - for debugging only
% if isfield(M,'IS')
%     Yhat = feval(M.IS,post.Ep,M,U);
% else
%     Yhat = spm_mci_fwd (post.Ep,M,U);
% end
% ey = Y - Yhat;
% sse = trace(M.iCe*ey'*ey);
% acc_gauss = -0.5*sse-0.5*M.N*M.logdet_Ce-0.5*M.N*log(2*pi);
% ew = x;
% comp_gauss = 0.5*ew'*M.ipC*ew-0.5*spm_logdet(M.ipC)-0.5*spm_logdet(Cp);
% logev_gauss = acc_gauss - comp_gauss;


