function [pout] = mci_gp_params (pin,M,U,rev)
% Convert parameters to quantities of interest
% FORMAT [pout] = mci_gp_params (pin,M,U,rev)
%
% pin       input parameters 
% M,U       as usual
% rev       if rev==1 do reverse transform (default=0)
%
% pout      output parameters
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

if nargin < 4 rev=0; end

switch M.prior
    case 'LogNormal',
        if rev
            pout = log(pin);
        else
            pout = exp(pin);
        end
        
    case 'Uniform',
        if rev
            r = M.sigma_max - M.sigma_min;
            pp = (pin-M.sigma_min)/r;
            pout = spm_invNcdf(pp,0,1);
        else
            pp = spm_Ncdf(pin,0,1);
            pout = M.sigma_min + pp*(M.sigma_max-M.sigma_min);
        end
        
    otherwise
        disp('Unknown prior in mci_gp_params.m');
end



