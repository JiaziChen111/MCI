function [M] = spm_mci_priors (M)
% Precompute quantities for evaluating log prior 
% FORMAT [M] = spm_mci_priors (M)
%
% M.subspace        set to 1 for subspace to be eigenvectors of prior cov
%                   Otherwise (default) use original parameterisation (M.V=I)
%
% M.V               projection matrix
% M.ipC             Inverse prior cov in reduced space
% M.log_prior_t2    second term of log prior 
% M.Np              dimension of reduced space
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_priors.m 6548 2015-09-11 12:39:47Z will $

try subspace=M.subspace; catch subspace=0; end
    
if isstruct(M.pC)
    pC=full(diag(spm_vec(M.pC)));
else
    pC = M.pC;
end

if subspace
    V  = spm_svd(pC,exp(-32));
    V  = full(spm_svd(pC,exp(-32)));
else
    V=eye(size(pC));
end

Np = size(V,2);
pC = V'*pC*V;
ipC = inv(pC);
log_prior_t2 = -spm_logdet(pC)/2-0.5*Np*log(2*pi);

M.ipC=ipC;
M.V=V;
M.log_prior_t2=log_prior_t2;
M.Np=Np;
