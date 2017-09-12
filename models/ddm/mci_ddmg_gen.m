function [y] = mci_ddmg_gen (P,M,U)
% Generate data from DDM using Gaussian approx  
% FORMAT [y] = mci_ddmg_gen (P,M,U)
%
% P         parameters
% M         model structure
% U         input structure (empty)
%
% y         [x; t] where x is binary decision, t is reaction time
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

[m,v,p] = mci_ddm_moments (P,M,U);

x = p > rand(1,1);
t = spm_normrnd(m,v,1);
y = [x; t];