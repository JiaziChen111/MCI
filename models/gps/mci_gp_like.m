function [L,yhat,st] = mci_gp_like (P,M,U,Y)
% Log-likelihood for Gaussian Process model 
% FORMAT [L,yhat,st] = mci_gp_like (P,M,U,Y)
%
% P         parameters
% M,U,Y     as usual
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

% Status flag (only used for dynamic systems)
st = [];

[yhat,C,K] = mci_gp_gen (P,M,U);

iC = spm_inv(C);
N = length(Y);
L = -0.5*spm_logdet(C) - 0.5*Y'*iC*Y-(N/2)*log(2*pi);

