function [P] = mci_ddm_to_mci (ddm,M,U)
% DDM to MCI parameters 
% FORMAT [P] = mci_ddm_to_mci (ddm,M,U)
%
% ddm   structure with fields
%   .v         slope
%   .a         decision interval
%   .b         non-decision time 
%   .r         bias (relative)
%
% P         parameters
% M         model structure
% U         input structure (empty)
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

zv = (ddm.v-M.v_min)/(M.v_max-M.v_min);
P(1) = spm_invNcdf(zv,0,1);

za = (ddm.a-M.a_min)/(M.a_max-M.a_min);
P(2) = spm_invNcdf(za,0,1);

if M.bexp
    % Exponential transform for b parameter
    P(3) = log(ddm.b);
else
    zb = (ddm.b-M.b_min)/(M.b_max-M.b_min);
    P(3) = spm_invNcdf(zb,0,1);
end

if isfield(ddm,'r')
    zr = (ddm.r-M.r_min)/(M.r_max-M.r_min);
    P(4) = spm_invNcdf(zr,0,1);
end

