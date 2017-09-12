function [post] = spm_mci_post (mcmc,M,U,Y,true_P)
% Estimate posterior density
% FORMAT [post] = spm_mci_post (mcmc,M,U,Y,true_P)
%
% mcmc          .inference = 'amc','ais','vl' or 'langevin' 
%               .verbose = 0 or 1 to plot progress (default 0)
%               .maxits = max number of iterations for sampling
%               .init = init parameter values (default is prior mean)
%               .data_fit = 0 or 1 to get fit to data (default 1)
%
% M             model structure
% U             inputs (shouldn't be empty)
% Y             data
% true_P        true parameters (if known)
%
% post          structure containing posterior (mean, samples etc)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_post.m 6548 2015-09-11 12:39:47Z will $

try verbose=mcmc.verbose; catch verbose=0; end
try data_fit=mcmc.data_fit; catch data_fit=0; end
    
if nargin < 5 | isempty(true_P)
    tp=0;
else
    tp=1;
end

tic;
switch mcmc.inference,
    
    case 'laplace',
        disp('Laplace Approximation');
        
        tic;
        [post,M] = spm_mci_laplace (mcmc,M,U,Y);
        
    case 'ais',
        disp('Annealed Importance Sampling');
        disp(' ');
        
        tic;
        post = spm_mci_ais (mcmc,M,U,Y);
        toc
        
        Nsamp=size(post.P,2);
        post.ind=[1:Nsamp];
        
        % Using weighted mean ********
        post.Ep=post.wEp;
        %post.Ep=mean(post.P(:,post.ind)')';
        
        post.mcmc=mcmc;

        % Eigenrep needed for spm_mci_joint later
        M = spm_mci_minit (M);
        
    case 'multi-amc',
        disp('Multiple chains of Adaptive Monte Carlo');
        disp(' ');
        Nsamp=ceil(0.5*mcmc.J);
        Nscale=0;
        Ntune=Nsamp;
        mc = spm_mci_popdef (Nscale,Ntune,Nsamp);
        mc.verbose=verbose;
        
        MM{1}=M;
        UU{1}=U;
        
        for it=1:mcmc.maxits,
            % Loop over chains
            mcmc.init{1}=spm_normrnd(M.pE,M.pC,1);
            Psamp = spm_mci_pop (mc,MM,UU,Y);
            if it ==1
                P=Psamp{1}.theta;
            else
                P=[P,Psamp{1}.theta];
            end
        end
        post.ind=[1:size(P,2)];
        post.Ep=mean(P(:,post.ind)')';
        post.P=P;
        post.mcmc=mcmc;
        
        % Eigenrep needed for spm_mci_joint later
        M = spm_mci_minit (M);
        
    case 'amc4',
        
        [M,stats] = spm_mci_amc (mcmc,M,U,Y);
        
        Psamp=stats.P';
        Nsamp=size(Psamp,1);
        post=stats;
        post.ind=[mcmc.Ntune+1:Nsamp];
        post.targ=Psamp(post.ind,:);
        post.Ep=mean(post.targ)';
        post.C=stats.C;
        
        logev.phm = spm_mci_phm (post.L2(post.ind));
        post.logev = logev;
        post.ar = sum(post.acc(post.ind))/length(post.ind);
        
    case 'amc',
        disp('Adaptive Monte Carlo');
        disp(' ');
        
        N1 = ceil(0.25*mcmc.maxits);
        N2 = ceil(0.5*mcmc.maxits);
        try Nscale = mcmc.Nscale; catch Nscale = N1; end
        try Ntune = mcmc.Ntune; catch Ntune = N1; end
        try Nsamp = mcmc.Nsamp; catch Nsamp = N2; end
        mc = spm_mci_popdef (Nscale,Ntune,Nsamp);
        mc.verbose=verbose;
        mc.anneal='frozen';
        
        % Draw initialisation point from prior ?
        %mcmc.init{1}=spm_normrnd(M.pE,M.pC,1);
        try mc.init{1}=mcmc.init; catch mc.init{1}=spm_vec(M.pE); end
        
        MM{1}=M;
        UU{1}=U;
        tic;
        [Psamp,logev,D,MM] = spm_mci_pop (mc,MM,UU,Y);
        toc
        M=MM{1};
        
        if tp
            [post.Ep,post.SDp]=spm_mci_report (Psamp,mc,true_P);
        else
            disp('Initial params:');
            disp(mc.init{1});
            [post.Ep,post.SDp]=spm_mci_report (Psamp,mc);
        end
        
        post.Ep=post.Ep';
        post.P=Psamp{1}.theta;
        post.E=-Psamp{1}.logq;
        post.dL=Psamp{1}.dL;
        post.bayes_fb=zeros(1,length(post.E));
        post.acc=Psamp{1}.acc;
        post.D=D;
        post.C=Psamp{1}.S; % Covariance of Gaussian proposal density
        post.mcmc=mcmc;
        post.ind=mc.ind_samp;
        post.logev=logev;
        post.ar = sum(post.acc(post.ind))/length(post.ind);
        
    case 'vl',
        disp('Variational Laplace');
        disp(' ');
        
        D.y=Y;
        if ~verbose
            M.nograph=1;
        end
        
        if isstruct(U)
            UI=U;
        else
            UI.u=U';
            UI.dt=M.T/M.N;
        end
        
        % Change VL defaults
        if isfield(mcmc,'maxits'), M.Nmax=mcmc.maxits; end
        if isfield(mcmc,'init'), M.P=mcmc.init; end
        if isfield(mcmc,'hE'), M.hE=mcmc.hE; end
        if isfield(mcmc,'hC'), M.hC=mcmc.hC; end
        
        %[Ep,Cp,Eh,F,tmp1,tmp2,tmp3,k] = spm_nlsi_GN (M,UI,D);
        vl = spm_mci_vl(M,UI,D);
        
        post.Ep=spm_vec(vl.Ep);
        post.Cp=vl.Cp;
        post.hE=vl.hE; post.hC=vl.hC;
        post.Eh=vl.Eh; post.Ch=vl.Ch;
        post.Ce=diag(1./exp(vl.Eh));
        post.logev=vl.F;
        post.its=vl.k;
        if isfield(M,'P')
            post.init=M.P;
        end
        
        % Eigenrep needed for spm_mci_joint later
        M = spm_mci_minit (M);

    case 'langevin',
        disp('Langevin Monte Carlo');
        disp(' ');
        
        [M,stats]=spm_mci_lgv(mcmc,M,U,Y);
        Psamp=stats.P';
        Nsamp=size(Psamp,1);
        try burn_in=mcmc.burn_in; catch burn_in=round(0.3*Nsamp); end
        post=stats;
        post.ind=[burn_in+1:Nsamp];
        post.targ=Psamp(post.ind,:);
        post.Ep=mean(post.targ)';
        
        % 2.5%, 50% and 97.5% quantiles
        q = [.025 .5 .975];
        for j=1:M.Np,
            sx = post.P(j,post.ind);
            post.quantiles(j,:) = quantile(sx,q);
        end
        
    otherwise
        disp('Unknown inference method');
end

if strcmp(mcmc.inference,'VL')
    Up=UI;
else
    Up=U;
end

if data_fit
    % Generate data fit from posterior mean
    if isfield(M,'IS')
        if strcmp(M.IS,'spm_gen_erp')
            Ps=spm_unvec(post.Ep,M.pE);
            post.Yhat=feval('spm_gen_erp',Ps,M,Up);
        else
            post.Yhat = feval(M.IS,post.Ep,M,Up);
        end
    else
        post.Yhat = spm_mci_fwd (post.Ep,M,Up);
    end
end

post.els=toc;
disp(sprintf('Optimisation time = %1.2f seconds',post.els));

lw=2;
if isfield(M,'t') & isfield(post,'Yhat') & verbose
    % Plot time series
    figure
    rm=ceil(sqrt(M.l));
    for i=1:M.l,
        if M.l>3
            subplot(rm,rm,i);
        else
            subplot(M.l,1,i);
        end
        plot(M.t,Y(:,i),'LineWidth',lw);
        hold on
        plot(M.t,post.Yhat(:,i),'r','LineWidth',lw);
        grid on
        set(gca,'FontSize',16);
        legend('Data','Fit');
        xlabel('Time');
        ylabel(sprintf('y(%d)',i));
    end
end

% get parameters in reduced space
Pr=M.V'*(post.Ep-M.vpE);
post.L_est = spm_mci_joint (Pr,M,U,Y);
disp(sprintf('Estimated Log Joint=%1.2f',post.L_est));

if tp    
    % get parameters in reduced space
    Pr=M.V'*(spm_vec(true_P)-M.vpE);
    post.L_true = spm_mci_joint (Pr,M,U,Y);
    disp(sprintf('True Log Joint=%1.2f',post.L_true));
    disp(' ');
end

pt=4;
if M.Np > pt & verbose
    hp=figure;
    set(hp,'Name','Parameters');
    plot(post.Ep,'r','LineWidth',lw);
    xlabel('Parameter');
    set(gca,'FontSize',16);
    grid on
else
    disp('Estimated (latent) params:');
    disp(post.Ep);
end

if tp & verbose
    if M.Np > pt
        hold on
        plot(spm_vec(true_P),'LineWidth',lw);
        legend('Estimated','True');
    else
        disp('True (latent) params:');
        disp(spm_vec(true_P));
    end
end

switch mcmc.inference,
    case {'amc','langevin'},
        post.type='sample';
    otherwise
        post.type='gaussian';
end

post.M=M;
post.U=U;