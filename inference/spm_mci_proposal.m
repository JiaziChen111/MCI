function [x] = spm_mci_proposal (M,U,m,C,shrink)
% Sample from Gaussian proposal checking for stability
% FORMAT [x] = spm_mci_proposal (M,U,m,C,shrink)
%
% M         model
% U         inputs
% m         proposal mean
% C         proposal covariance
% shrink    1 or 0 (see below)
%
% x         sample
%
% For shrink=1 we sequentially halve the difference 
% between initial sample and prior mean, until stability.
%
% For shrink=0 we generate up to stabmax random samples until
% model produces stable (within-range) predictions. 
%
% Only check for stability if M.stabfun is specified.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

if nargin < 5, shrink=1; end

if ~isfield(M,'stabfun')
    x = spm_normrnd(m,C,1);
    return
end

stab=0;
nstab=0;
stabmax=16;

if shrink
    % Shrinkage method, assumes m is stable
    xinit = spm_normrnd(m,C,1);
    dx = xinit-m;
    while stab==0 & nstab < stabmax
        x = m+dx;
        stab = feval(M.stabfun,x,M,U);
        dx = dx/2;
        nstab = nstab+1;
    end
else
    % Sampling method
    while stab==0 & nstab < stabmax
        x = spm_normrnd(m,C,1);
        stab = feval(M.stabfun,x,M,U);
        nstab = nstab+1;
    end
end

if nstab == stabmax
    %disp('Cannot produce within-range predictions');
    %disp('Increase stabmax in spm_mci_proposal.m');
    x=m;
end