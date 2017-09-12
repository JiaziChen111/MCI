function [logev,logpost] = spm_mci_chib (post,M,U,Y)
% Compute Log Evidence using Chib's method
% FORMAT [logev,logpost] = spm_mci_chib (post,M,U,Y)
%
% post      data structure from spm_mci_post.m
%           .P   posterior samples
%           .C   covariance of proposal density
% M         model structure
% U         input structure
% Y         data structure
% 
% logev     Chib approximation to log model evidence
% logpost   Chib approximation to log posterior mean 
%
% This implementation assumes posterior samples were created using
% Metropolis-Hastings with a fixed proposal density being 
% a Gaussian kernel with covariance C. The single point used to
% implement the method is chosen to be the posterior mean.
%
% [1] S Chib and I Jeliazkov (2001) Marginal Likelihood From
% the Metropolis-Hastings Output. JASA, 96(453), pp 270-281.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

% Covariance of MH proposal density
C = post.C;
Np = size(C,1);
C = C + eps*eye(Np);

% Indices of samples from posterior
pind = post.ind;
Npost = length(pind);

% Single point used
P = post.Ep;
Lp = post.L_est;

% Posterior density at P from equation 9 in [1]
Lsamp = -post.E(pind);

log_alpha = min([zeros(1,Npost); Lp-Lsamp]);

z = post.P(:,pind);
z = z - P(:)*ones(1,size(z,2));
iC = spm_inv(C);
for i=1:Npost,
    logq (i) = -0.5*z(:,i)'*iC*z(:,i);
end
logq = logq - ones(1,Npost)*(0.5*Np*log(2*pi)+0.5*spm_logdet(C));
%numerator = mean(exp(log_alpha+logq))+eps;
k = max(log_alpha+logq);
log_numerator = log(sum(exp(log_alpha+logq-k)))+k-log(Npost);

try V=M.V; catch M = spm_mci_minit (M); end
wdenom = spm_normrnd(P,C,Npost);
for j = 1:Npost,
    Pr = M.V'*(wdenom(:,j)-M.vpE);
    Ldsamp(j) = spm_mci_joint (Pr,M,U,Y);
end
%alpha = min([ones(1,Npost);exp(Ldsamp-Lp)]);
%denominator = mean(alpha)+eps;
%logpost = log(numerator)-log(denominator);

log_alpha = min([zeros(1,Npost);Ldsamp-Lp]);
k = max(log_alpha);
log_denominator = log(mean(exp(log_alpha-k)))+k;

logpost = log_numerator-log_denominator;

% Log evidence from equation 4 in [1]
% (logjoint is loglike + logprior)
logev = Lp - logpost;

