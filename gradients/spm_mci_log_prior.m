function [L1] = spm_mci_log_prior (Pr,M)
% Log Prior Prob of model M
% FORMAT [L1] = spm_mci_log_prior (Pr,M)
% Pr        parameters (vectorised and in M.V subspace)
% M         model structure. 
%
% L1        log prior, L1 = log p(Pr|M)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

e = Pr;
L1 = - e'*M.ipC*e/2 + M.log_prior_t2;