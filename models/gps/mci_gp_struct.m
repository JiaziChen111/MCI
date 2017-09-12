function [M,U] = mci_gp_struct (M,U)
% Gaussian Process data structure
% FORMAT [M,U] = mci_gp_struct (M,U)
%
% M         .kernel   'Gaussian'
%           .Ce  Error covariance
%           .T   Number of data points
% U         Input structure
%
% M         Model structure
% U         Input structure
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

M.l=1; % Single output variable

% Get Euclidean distances amongst inputs
for i=1:M.N,
    for j=1:M.N,
        U.E(i,j) = norm(U.X(i,:)-U.X(j,:));
    end
end

M.L='mci_gp_like';
M.IS='mci_gp_gen';
%M.hess_outer = 1;

switch M.kernel
    case 'Gaussian',
        M.pE=0';
        M.pC=1;
    otherwise
        disp('Unknown kernel type in mci_gp_struct.');
end

M.logdet_Ce=spm_logdet(M.Ce);
M.iCe=inv(M.Ce);
