function [E] = spm_mci_energyP (P,M,U,Y)
% Compute energy of model
% FORMAT [E] = spm_mci_energyP (P,M,U,Y)
%
% P     parameters 
% M     model structure
% U     inputs
% Y     data
%
% E     energy (negative log joint)
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

Pr = M.V'*(p-M.vpE);
E = spm_mci_energy(Pr,M,U,Y);


