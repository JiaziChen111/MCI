function [M,U] = mci_ddm_struct_old (ddm)
% Set up data structure for drift diffusion model
% FORMAT [M,U] = mci_ddm_struct_old (ddm)
%
% ddm       'biased' or 'unbiased' (default)
%
% M         model structure
% U         input structure (empty)
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

try ddm=ddm; catch ddm='unbiased'; end

% Assumed noise variance of diffusion
% Navarro and Fuss, Wiecki et al assume s2=1
% Ratcliff, Wagenmackers et al assume s2=0.01
M.s2=1; 
scale=sqrt(M.s2)/0.1;

% See Wagenmakers et al (2005) for typical values
% when noise variance, s2=0.01, noise SD, s=0.1
wagenmakers=0;

% Multiply priors over v and a by scale if assuming 
% different noise SD. See footnote 1 in Navarro and Fuss, 2009
if wagenmakers==1
    M.v_min=0.1*scale; M.v_max=0.5*scale;
    M.a_min=0.07*scale; M.a_max=0.17*scale;
else
    M.v_min=0; M.v_max=1.5;
    M.a_min=1; M.a_max=4;
end
M.b_min=0.3; M.b_max=0.7;


M.L='mci_ddm_like';
%M.hess_outer=1;

% Acceptable truncation error in Navarro-Fuss expansion
M.truncation_error=10^(-3);

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
    case 'aconflict',
        % v, a_LL, a_WW, a_WL, b 
        U.con{1}=[1 0 0 0 0; 0 1 0 1 0; 0 0 0 0 1]; % LL
        U.con{2}=[1 0 0 0 0; 0 0 1 1 0; 0 0 0 0 1]; % WW
        U.con{3}=[1 0 0 0 0; 0 0 0 1 0; 0 0 0 0 1]; % WL
        M.Np=5;
    case 'theta',
        % v, a0, beta, b
        M.Np=4;
    case 'thetaDBS',
        % v, a0, betaOFF, betaON, b
        M.Np=5;
    otherwise
        M.Np=3;
end
M.pE = zeros(M.Np,1);
M.pC = eye(M.Np);

% Number of outputs (decision and reaction time)
M.l=2;



