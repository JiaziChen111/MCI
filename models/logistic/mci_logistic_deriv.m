function [dLdp,iCpY,L] = mci_logistic_deriv (P,M,U,Y)
% Gradient of likelihood for logistic model
% FORMAT [dLdp,iCpY,L] = mci_logistic_deriv (P,M,U,Y)
%
% P         parameters
% M         model
% U         inputs
% Y         data
%
% dLdp      gradient of log joint
% iCpY      curvature (Fisher Information)
% L         log joint
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_logistic_deriv.m 6548 2015-09-11 12:39:47Z will $

g = mci_logistic_gen (P,M,U);
X = U.X(1:M.N,:);

if isstruct(Y)
    y=Y.y;
else
    y=Y;
end

% Gradient
dLdp = X'*(y-g);

if nargout > 1
    % Curvature
    Lambda=diag(g.*(1-g));
    iCpY = X'*Lambda*X;
end

if nargout > 2
    % Log Likelihood
    L = mci_logistic_like (P,M,U,Y);
end