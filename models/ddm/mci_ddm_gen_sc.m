function [L,E,st] = mci_ddm_like_fixed (P,M,U,Y)
% Compute log likelihood of fixed-effects DDM 
% FORMAT [L,E,st] = mci_ddm_like_fixed (P,M,U,Y)
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
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

st=0;

% tic
% N=length(Y);
% for n=1:N,
%     Lfirst(n) = mci_ddm_like (P,M,U,Y(n,:));
% end
% L = sum(Lfirst);
% %toc
% E = -L;

%tic
ddm = mci_ddm_from_mci (P,M,U);
p = mci_ddm_wfpt_vec (ddm,M,U,Y);
L = sum(log(p));
E = -L;
%toc

% L 
% Lchk
% keyboard