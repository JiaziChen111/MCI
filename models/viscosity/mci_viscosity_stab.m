function [stab,ymin,ymax] = mci_viscosity_stab (Pr,M,U)
% Do parameters produce within-range (sensible) predictions
% FORMAT [stab,ymin,ymax] = mci_viscosity_stab (Pr,M,U)
%
% Pr         parameters (in reduced space)
% M,U       as usual
%
% stab      1 for yes, 0 for no
% ymin      minimum predicted value
% ymax      maximum predicted value
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

% Parameters in original space
P = M.V*Pr+M.vpE;

y = mci_viscosity_gen(P,M,U);
ymax = max(y);
ymin = min(y);
stab = (ymax < 20) & (ymin > 2);
