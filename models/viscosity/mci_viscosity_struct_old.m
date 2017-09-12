function [M,U] = mci_viscosity_struct (X,model,stab)
% Viscosity model structure
% FORMAT [M,U] = mci_viscosity_struct (X,model,stab)
%
% X         X(:,1) temperature
%           X(:,2) pressure
% model     'full' or 'reduced'
% stab      use priors that produce within-range predictions
%           (found using mci_viscosity_priors.m)
%
% M         Model structure
% U         Input structure
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

try model=model; catch model='full'; end
try stab=stab; catch stab=0; end

M.l=1; % Single output variable

M.T=size(X,1);
U.X=X;
M.N=M.T;
M.t=[1:M.T];

% Indices of different temperatures
U.i{1}=[1:11];
U.i{2}=[12:23];
U.i{3}=[24:38];
U.i{4}=[39:53];

M.L='mci_viscosity_like';
M.IS='mci_viscosity_gen';
M.stabfun='mci_viscosity_stab';

% Prior is set to what we know about parameters before
% nonlinear estimation - see page 88/89 Bates and Watts
M.pE(1)=983;
M.pE(2)=192;
M.pE(3)=1.35;
M.pE(4)=0;

slarge=10;
ssmall=0.03;
switch model
    case 'full',
        Np=9;
        s2=ones(Np,1);
        s2(1)=slarge;
        s2(2)=slarge;
        s2(8)=slarge;
        s2(3)=0.5;
        s2(4)=ssmall;
        s2(5)=ssmall/2;
        s2(7)=ssmall;
        
        M.pE(5)=0;
        M.pE(6)=0.22;
        M.pE(7)=0;
        M.pE(8)=35.9;
        M.pE(9)=0;
        
    case 'reduced',
        % Cubic terms (params 5 and 7) removed
        Np=7;
        s2=ones(Np,1);
        s2(1)=slarge;
        s2(2)=slarge;
        s2(3)=0.5;
        s2(6)=slarge;
        s2(4)=ssmall;
        
        M.pE(5)=0.22;
        M.pE(6)=35.9;
        M.pE(7)=0;
end

M.pE=M.pE(:);
M.pC=diag(s2.^2);

if stab
    load viscosity_stable_priors
    %M.pE=stable.pE;
    M.pC=stable.pC;
    if strcmp(model,'reduced')
        keep=[1,2,3,4,6,8,9];
        %M.pE=M.pE(keep);
        M.pC=M.pC(keep,keep);
    end
end

% page 89
%rss=0.08996;
%sigma_e=sqrt(rss/M.T);

%sigma_e=0.01;

%sigma_e=0.08;
sigma_e=0.2;

M.Ce=sigma_e^2;
M.logdet_Ce=spm_logdet(M.Ce);
M.iCe=inv(M.Ce);



