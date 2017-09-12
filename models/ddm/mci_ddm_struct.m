function [M,U] = mci_ddm_struct (ddm,bexp)
% Set up data structure for drift diffusion model
% FORMAT [M,U] = mci_ddm_struct (ddm,bexp)
%
% ddm       'biased' or 'unbiased' (default)
% bexp      exponential parameterisation for b (default=0)
%
% M         model structure
% U         input structure (empty)
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

try ddm=ddm; catch ddm='unbiased'; end
if nargin < 2, bexp=0; end

% Assumed noise variance of diffusion process
M.s2=1;

M.v_min=0; M.v_max=1.5;
M.a_min=1; M.a_max=3;
M.b_min=0.3; M.b_max=0.7;

M.L='mci_ddm_like';
M.hess_outer=0;

% Acceptable truncation error in Navarro-Fuss expansion
M.truncation_error=10^(-3);


v=0.25; % Prior variance
%v=1; % Prior variance
%r=-0.05; % Prior correlation
r=0; % Prior correlation
        
U=[];
switch ddm
    case 'biased',
        M.r_min=0.4; M.r_max=0.6;
        M.Np=4;
    case 'vconflict',
        % v_LL, v_WW, v_WL, a, b 
        U.con{1}=[1 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1]; % LL
        U.con{2}=[0 1 1 0 0; 0 0 0 1 0; 0 0 0 0 1]; % WW
        U.con{3}=[0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1]; % WL
        M.Np=5;
    case 'vconflict2',
        % v_1, v_2, a, b 
        U.con{1}=[1 0 0 0; 0 0 1 0; 0 0 0 1]; % C1
        U.con{2}=[1 1 0 0; 0 0 1 0; 0 0 0 1]; % C2
        M.Np=4;
        M.pC=r*v*ones(4,4);
        for i=1:4; M.pC(i,i)=v; end
        M.pE = zeros(M.Np,1);
    case 'aconflict',
        % v, a_LL, a_WW, a_WL, b 
        U.con{1}=[1 0 0 0 0; 0 1 0 1 0; 0 0 0 0 1]; % LL
        U.con{2}=[1 0 0 0 0; 0 0 1 1 0; 0 0 0 0 1]; % WW
        U.con{3}=[1 0 0 0 0; 0 0 0 1 0; 0 0 0 0 1]; % WL
        M.Np=5;
    case 'theta',
        % v, a0, beta, b
        M.Np=4;
        M.pE=zeros(M.Np,1);
        M.pC=r*v*ones(M.Np,M.Np);
        for i=1:M.Np; M.pC(i,i)=v; end
        %M0.pC(3,3)=1;
    case 'thetaDBS',
        % v, a0, betaOFF, betaON, b
        M.Np=5;
        M.pE=zeros(M.Np,1);
        M.pC=r*v*ones(M.Np,M.Np);
        for i=1:M.Np; M.pC(i,i)=v; end
        %M0.pC(3,3)=1;
        %M0.pC(4,4)=1;
    otherwise
        M.Np=3;
        M.pC=r*v*ones(3,3);
        M.pC(1,1)=v;
        M.pC(2,2)=v;
        M.pE = zeros(M.Np,1);
        if bexp
            M.pE(3)=log(0.3);
            M.pC(3,3)=1/32;
        else
            M.pC(3,3)=v;
        end
end

M.bexp=bexp;

% Number of outputs (decision and reaction time)
M.l=2;



