function [L,E,st] = mci_ddm_like (P,M,U,Y)
% Compute log likelihood of unbiased DDM
% FORMAT [L,E,st] = mci_ddm_like (P,M,U,Y)
%
% P         parameters
% M         model
% U         inputs
% Y         data
% 
% L         Log likelihood
% E         Errors
% st        Status flag (0 for OK, -1 for problem)
%
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

pr = mci_ddm_wfpt(P,M,U,Y);
L = log(pr);
E = -L;
st = 0;