function [logev,q] = spm_mci_ais_evidence (logw)
% Compute model evidence from importance weights
% FORMAT [logev,q] = spm_mci_ais_evidence (logw)
%
% logw      log importance weights over trajectories
% 
% logev     log evidence
% q         normalised importance weights
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

r=max(logw);
u=exp(logw-r);
S=mean(u);
logev=log(S)+r;
q=u/sum(u);

