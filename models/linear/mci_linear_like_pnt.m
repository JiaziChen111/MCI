function [L,E,st] = mci_linear_like_pnt (theta,M,U,Y,n)
% Compute log likelihood of single data point for linear model 
% FORMAT [L,E,st] = mci_linear_like_pnt (theta,M,U,Y,n)
%
% theta     regression coefficients
% M         model
% U         inputs
% Y         data
% n         nth data point
% 
% L         Log likelihood
% E         Errors
% st        Status flag (0 for OK, -1 for problem)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

st=0;

yhat=U.X(n,:)*theta(:);
if isstruct(Y)
    E=sum(sum((Y.y(n)-yhat).^2));
else
    E=sum(sum((Y(n)-yhat).^2));
end

L = -0.5*M.logdet_Ce - 0.5*log(2*pi);
L = L - 0.5*M.iCe*E;
