function [mt,vt,pc] = mci_ddm_moments (P,M,U)
% Generate moments for unbiased drift diffusion model 
% FORMAT [mt,vt,pc] = mci_ddm_moments (P,M,U)
%
% P         parameters
% M         model structure
% U         input structure (empty)
%
% mt        mean reaction time
% vt        variance of reaction time
% pc        probability of correct decision
%
% See Wagenmakers et al (2005) On the relation between the
% mean and the variance of a diffusion model response time
% distribution. Journal of Mathematical Psychology, 
% 49:195-204.
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

ddm = mci_ddm_from_mci (P,M,U);
v=ddm.v;
a=ddm.a;
b=ddm.b;

% Assumed noise variance of diffusion process
s2 = M.s2;

k = -v*a/s2;

ek = exp(k);
e2k = exp(2*k);
t1 = a/(2*v);
t2 = s2/(v^2);
denom = (1+ek)^2;

pc = 1/(1+exp(k));
mt = b + t1*((1-ek)/(1+ek));
vt = t1 * t2 * (2*k*ek-e2k+1)/denom;



