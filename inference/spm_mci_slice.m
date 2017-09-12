function [post,M] = spm_mci_slice (mcmc,M,U,Y)
% Slice Sampling
% FORMAT [post,M] = spm_mci_slice (mcmc,M,U,Y)
%
% mcmc  Optimisation parameters
%       .maxits = number of samples
%       .burnin = the number to throw away
%       .w = vector of widths (default is prior SD)
%
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
% Uses the slicesample.m function from Matlab's stats toolbox

% Will Penny
% $Id$

try maxits=mcmc.maxit; catch maxits=4096; end
try burnin=mcmc.burnin; catch burnin=2048; end
%try w=mcmc.w; catch w=sqrt(diag(M.pC))'; end
try w=mcmc.w; catch w=0.5; end
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
else
    % User specified or prior mean
    try init=mcmc.init; catch init=M.vpE; end
    
    % Initial param in eigenspace
    xinit = M.V'*(init-M.vpE);
end

% Slice sampling
f = @(x) spm_mci_joint (x,M,U,Y);
x = slicesample(xinit',maxits,'logpdf',f,'width',w,'burnin',burnin);

% Project parameters back from eigenspace into original space
x=ones(maxits,1)*M.vpE'+x*V;

post.P = x';
post.Ep = mean(x)';
post.Cp = cov(x);

Epsub=M.V'*(post.Ep-M.vpE);
post.L = spm_mci_joint(Epsub,M,U,Y);
post.E = -post.L;

% Model Evidence 
post.logev = [];

