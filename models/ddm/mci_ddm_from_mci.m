function [ddm] = mci_ddm_from_mci (P,M,U)
% DDM from MCI parameters 
% FORMAT [ddm] = mci_ddm_from_mci (P,M,U)
%
% P         parameters
% M         model structure
% U         input structure (empty)
%
% ddm   structure with fields
%   .v         slope
%   .a         decision interval
%   .b         non-decision time 
%   .r         bias (relative)
%
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

% Standard Gaussian CDF
u = 0.5+erf(P/sqrt(2))/2;

ddm.v = u(1)*(M.v_max-M.v_min)+M.v_min;
ddm.a = u(2)*(M.a_max-M.a_min)+M.a_min;
if M.bexp
    % Exponential transform for b parameter
    ddm.b = exp(P(3));
else
    ddm.b = u(3)*(M.b_max-M.b_min)+M.b_min;
end

if length(P) > 3
    ddm.r = u(4)*(M.r_max-M.r_min)+M.r_min;
end



