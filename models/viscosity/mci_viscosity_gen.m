function [y] = mci_viscosity_gen (P,M,U)
% Viscosity model
% FORMAT [y] = mci_viscosity_gen (P,M,U)
%
% P         parameters
% M,U       as usual
%
% See section 6.2 in 
% DiCiccio at al, JASA, 92(439):903-915.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

x1 = U.X(:,1);
x2 = U.X(:,2);
x2s = x2.^2;

Np=length(P);
if Np==9
    % Full model
    t1 = P(1)./(P(2)+x1);
    t2 = P(3)*x2;
    t3 = P(4)*x2s;
    t4 = P(5)*x2.^3;
    m = (P(6)+P(7)*x2s).*x2;
    a = -x1./(P(8)+P(9)*x2s);
    t5 = m.*exp(a);
    y = t1+t2+t3+t4+t5;
else
    % Np==7: Reduced Model
    t1 = P(1)./(P(2)+x1);
    t2 = P(3)*x2;
    t3 = P(4)*x2s;
    m = P(5)*x2;
    a = -x1./(P(6)+P(7)*x2s);
    t5 = m.*exp(a);
    y = t1+t2+t3+t5;
end

