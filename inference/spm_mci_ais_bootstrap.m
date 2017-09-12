function [boot] = spm_mci_ais_bootstrap (logw)
% Confidence interval on log evidence from bootstrapping
% FORMAT [boot] = spm_mci_ais_bootstrap (logw)
%
% logw      Log importance weights
%
% boot      .logev_low  5th percentile
%           .logev_high 95th percentile
%           .logev_se   standard deviation
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

Nboot=1000;
maxits=length(logw);

for b=1:Nboot,
    ind=ceil(rand(1,maxits)*maxits);
    logw_resample=logw(ind);
    logev_resample(b)=spm_mci_ais_evidence (logw_resample);
end

boot.logev_se=std(logev_resample);

lind=ceil(0.05*Nboot);
uind=floor(0.95*Nboot);
lsort=sort(logev_resample);

boot.logev_low=lsort(lind);
boot.logev_high=lsort(uind);
