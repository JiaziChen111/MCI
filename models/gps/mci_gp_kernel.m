function [K] = mci_gp_kernel (P,M,U)
% Gaussian Process Kernel
% FORMAT [K] = mci_gp_kernel (P,M,U)
%
% P         parameters
% M,U       as usual
%
% K         N x N kernel matrix where N is number of data points
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

switch M.kernel
    case 'Gaussian',
        sigma = mci_gp_params(P,M,U);
        K = exp(-U.E/(2*sigma^2)); % Eq 6.23 Bishop, 2006
        
    otherwise
        disp('Unknown kernel function in mci_gp_kernel.m');
end

