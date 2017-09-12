function [M,stats] = spm_mci_amc (mcmc,M,U,Y)
% Sampling using Adaptive Monte Carlo
% FORMAT [M,stats] = spm_mci_amc (mcmc,M,U,Y)
%
% mcmc   Sampling parameters
%        .verbose            display progress
%        .maxits             maximum number of total samples 
%        .init               initial sample values (start of chain)
%        .h                  step size
%        .C                  covariance of proposal density 
%                            (default is prior covariance)
%        .Nscale,  .Ntune    number of samples in adaptive phase
%   
% M      Model Structure
% U      Inputs
% Y      Data
%
% M      Updated model structure
% stats  Structure with fields:
%
% .P     Samples, [maxits x M.Np] 
% .E     Negative log joint prob, [maxits x 1]
%
% Uses Algorithm 4 from [1]. See also equations 4 and 5 in [2].
%
% If .Nscale=0, .Ntune=0 this will run MH (non-adaptive) with 
% proposal density mcmc.C 
% 
% [1] C Andrieu and J Thoms (2008). A tutorial on adaptive MCMC. Statistical 
% computing, 180(4), pp 343-373.
%
% [2] B Sengupta, K Friston and W Penny (2015) Gradient-free MCMC methods 
% for dynamic causal mdoelling. Neuroimage, 112, pp 375-381.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

% Defaults
try verbose=mcmc.verbose; catch verbose=0; end
try maxits=mcmc.maxits; catch maxits=2048; end
try h=mcmc.h; catch h=0.5; end 
try plot_int=mcmc.plot_int; catch plot_int=1; end
try init=mcmc.init; catch init=spm_vec(M.pE); end
try Ntune=mcmc.Ntune; catch Ntune=ceil(maxits/4); end
try Nscale=mcmc.Nscale; catch Nscale=ceil(maxits/4); end

% Compute eigen-parameterisation
M = spm_mci_minit (M);
V  = M.V;

% Initial param in eigenspace
xinit = M.V'*(init-M.vpE);

% Initialise proposal covariance
try Cp=mcmc.C; catch Cp = M.pC; end
Cp = M.V' * Cp * M.V;

% Sample matrix
x = zeros(maxits,M.Np);
x(1,:) = xinit';         

if verbose figure; end

% Initialise Robbins-Monro - see page 3 of [2]
logs = 0;
mu = zeros(M.Np,1);
C = Cp;
alpha_target=0.23;

% Main Loop
i=1; 
acc=zeros(maxits,1);
while (i <= maxits),
    
    if verbose
        if mod(i,plot_int) == 0 & i > 2
            spm_mci_progress (x,E,i);
        end
    end
    
    % Proposal (first proposal always accepted)
    if i==1
        pos=x(1,:)';
        curr=[];
    else
        pos=spm_normrnd(curr.mu,curr.Cp, 1);
    end
        
    % Quantities re proposal
    prop.pos=pos;
    [prop.L,prop.L2,st] = spm_mci_joint(pos,M,U,Y);
    if st==-1
        disp('Integration problem in spm_mci_amc.m');
        keyboard
    end
    
    prop.mu = pos;
    if i==1,
        prop.Cp = Cp;
    else
        prop.Cp = curr.Cp;
    end
    prop.iCp = zeros(M.Np); % Not used
    prop.logdetCp = 0; % Not used
    
    % Accept proposal ?
    [curr,accepted,bayes_fb(i),dL(i),alpha] = spm_mci_mh_update(curr,prop,verbose);
    acc(i)=accepted;
    E(i) = -curr.L;
    L2(i) = curr.L2;
    if i > 1
        dEdit(i-1)=100*(E(i)-E(i-1))/E(i-1);
    end
    x(i,:)=curr.pos;
    
    if i < (Ntune + Nscale)
        % Equation 31 in [1]
        h = 1/(i+1);
        logs = logs + h * (alpha-alpha_target);
        err = (curr.pos - mu);
        mu = mu + h * err;
        C = C + h * (err*err'-C);
        
        % Always update proposal covariance
        % whether or not sample was accepted
        curr.Cp = exp(logs)*C;
        
        % Record proposal variances
        % sp(:,i) = diag(curr.Cp);
    end
    
    i=i+1;
end

if verbose
    disp(sprintf('Total accepted samples = %d', sum(acc)));
end

% Project parameters back from eigenspace into original space
x=x(1:i-1,:);
nj=size(x,1);
stats.P=M.vpE*ones(1,nj)+V*x';
Cp = M.V * curr.Cp * M.V';

stats.E=E;
stats.dEdit=dEdit;
stats.acc=acc;
stats.bayes_fb=bayes_fb;
stats.dL=dL;
stats.L2=L2;
stats.C=Cp;
