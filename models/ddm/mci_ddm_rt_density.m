function [p] = mci_ddm_rt_density (P,M,U,y)
% Probability density of unbiased DDM 
% FORMAT [p] = mci_ddm_rt_density (P,M,U,y)
%
% P         parameters
% M         model
% U         inputs
% y         y(1)=x outcome, y(2)=t reaction time
% 
% p         joint probability density, p(y)=p(x,t)
%
% See equation 31 in Tuerlinckx (2004) The efficient
% computation of the cumulative distribution and 
% probability density functions in the diffusion model.
% Behaviour Research Methods, Instruments and Computers,
% 36:702-716.
%
% or equation 2 
% in Wagenmakers et al (2005)On the relation between the
% mean and the variance of a diffusion model response time
% distribution. Journal of Mathematical Psychology, 
% 49:195-204.
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

st=0;

% Number of terms in series expansion
% For t >> b we may need fewer terms
try nterms = M.nterms; catch nterms=32; end

ddm = mci_ddm_from_mci (P,M,U);
v=ddm.v;
a=ddm.a;
b=ddm.b;

t = y(2);

if ~(t > b)
    disp('Error in ddm_rt_density: time must be larger than non-decision time');
    disp(sprintf('t=%1.2f, b=%1.2f',t,b));
    keyboard
end

% Assumed noise variance of diffusion
s2 = M.s2;

r1 = s2/(a^2);
r2 = (v^2)/s2;
k = -v*a/s2;
L0 = log(pi) + log (r1) + 0.5*k;

ex=zeros(nterms+1,1);
for n=0:nterms,
    p1 = (2*n+1)*(-1)^n;
    pp2 = r2 + (pi^2)*((2*n+1)^2)*r1;
    p2 = exp(-0.5*pp2*(t-b));
    tmp = p1*p2;
    ex(n+1) = tmp;
    if ~isreal(tmp)
        keyboard
    end
end

sx = sum(ex);
if ~(sx > 0)
    L = L0 + log(eps);
else
    L = L0 + log(sx);
end

p = exp(L);

