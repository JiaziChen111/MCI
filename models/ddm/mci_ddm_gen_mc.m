function [L,E,st] = mci_ddm_like_subject (P,M,U,Y)
% Compute log likelihood of within-subject DDM 
% FORMAT [L,E,st] = mci_ddm_like_subject (P,M,U,Y)
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
% Allow for multiple trials and multiple trial types (conditions)
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

st=0;

% Loop over trial types
Nc=length(U.cond);
for c=1:Nc,
    % Parameters for each trial type are predefined
    % linear combination of overall parameters
    Ptr=U.con{c}*P;
    ddm = mci_ddm_from_mci (Ptr,M,U);
    p = mci_ddm_wfpt_vec (ddm,M,U,Y(U.cond{c},:));
    Lc(c) = sum(log(p));
end
L = sum(Lc);
E = -L;
