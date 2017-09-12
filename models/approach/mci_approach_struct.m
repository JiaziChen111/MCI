function [M,U] = mci_approach_struct (Nobs,model)
% Approach model structure
% FORMAT [M,U] = mci_approach_struct (Nobs,model)
%
% Nobs      Number of observations
% model     'full' (default) or 'reduced'
%
% M         Model structure
% U         Input structure
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_approach_struct.m 6548 2015-09-11 12:39:47Z will $

if nargin < 2
    model='full';
end

M.l=1; % Single output variable

M.T=40;
dt=M.T/(Nobs-1);
U.X=[0:dt:M.T]';
M.N=length(U.X);
M.t=U.X;

M.L='mci_approach_like';
M.IS='mci_approach_gen';
M.dL='mci_approach_deriv';

switch model
    case 'full',
        M.pE=[log(20),log(5)]';
        M.pC=[1/16 0; 0 1/16];
    case 'reduced',
        M.pE=log(20);
        M.pC=1/16;
    otherwise
        disp('Unknown model type in mci_approach_struct.');
end




