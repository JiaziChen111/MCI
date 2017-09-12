function [E,dEdp] = spm_mci_energy (Pr,M,U,Y,beta)
% Compute energy of model
% FORMAT [E,dEdp] = spm_mci_energy (Pr,M,U,Y,beta)
%
% Pr    parameters (vectorised and in M.V subspace)
% M     model structure
% U     inputs
% Y     data
% beta  inverse temperature
%
% E     energy (negative log joint)
% dEdp  gradient of Energy
%
% A default beta=1 gives usual energy and gradient
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

if nargin < 5 | isempty(beta)
    beta=1;
end

L = spm_mci_joint (Pr,M,U,Y,beta);

M.beta = beta;
dLdp = spm_mci_joint_grad (Pr,M,U,Y);
dLdp = dLdp(:);

E = -L;
dEdp = -dLdp;
