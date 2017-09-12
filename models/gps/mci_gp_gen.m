function [y,C,K] = mci_gp_gen (P,M,U)
% Gaussian Process generative model
% FORMAT [y,C,K] = mci_gp_gen (P,M,U)
%
% P         parameters
% M,U       as usual
%
% y         sample from GP model
% C         covariance function evaluated at data points
% K         kernel function evaluated at data points
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

switch M.kernel
    case 'Gaussian',
        K = mci_gp_kernel (P,M,U);
    otherwise
        disp('Unknown kernel function in mci_gp_gen.m');
end

% Covariance of observations
C = K + M.Ce * eye(M.N);

m = zeros(M.N,1);
y = spm_normrnd(m,C,1);

