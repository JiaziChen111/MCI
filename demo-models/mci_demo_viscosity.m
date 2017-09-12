
clear all
close all

% Load 'viscosity' data from Table A1.8, Logarithm of lubricant viscosity
% versus pressure and temperature, of [1] 
%
% [1] Bates and Watts, Nonlinear Regression
% and its Applications, 1988. 
% 
% See also section 6.2 in 
% [2] DiCiccio at al, JASA, 92(439):903-915.

load('viscosity_data_ord.txt','-ascii');

% Which inference algorithms ?
do_vl=0;
do_lgv=1;
do_ais=0;

X=viscosity_data_ord(:,2:3);
X(:,2)=X(:,2)/1000;

[M,U] = mci_viscosity_struct (X,'full',0);
[Malt,Ualt] = mci_viscosity_struct (X,'reduced',0);

Y = viscosity_data_ord(:,1);
yhat = mci_viscosity_gen (M.pE,M,U);

mcmc.data_fit=1;

if do_vl
    disp('Fitting model using Variational Laplace (VL)');
    mcmc.inference='vl';
    
    % These next 3 lines let VL use true noise variance
    noise_prec=1/(M.Ce(1,1));
    mcmc.hE=log(noise_prec);
    mcmc.hC=1e-8;
end

if do_lgv
    disp('Fitting model using Langevin Monte Carlo (LMC)');
    % Use outer product approximation for Hessian
    % M.hess_outer=1;
    % Malt.hess_outer=1;
    % M.hess_reg=0.01;
    % Malt.hess_reg=0.01;
    
    M.subspace=1;
    
    mcmc.inference='langevin';
    mcmc.lgv_update='alt';
    mcmc.maxits=2048;
    mcmc.verbose=0;
    %mcmc.h=0.1;
end

post1 = spm_mci_post (mcmc,M,U,Y);
post2 = spm_mci_post (mcmc,Malt,Ualt,Y);

if ~strcmp(mcmc.inference,'vl')
    diag.traceplot=1;
    diag.eplot=1;
    if strcmp(mcmc.inference,'langevin')
        diag.bplot=1;
    else
        diag.bplot=0;
    end
    spm_mci_diag(post1,diag);
    
    stats = spm_mci_mvnpost (post1,'ESS')
    stats = spm_mci_mvnpost (post1,'thinning')
    
    %spm_mci_quantiles (post1,2);
end

sse1=mci_plot_viscosity (U,X,Y,yhat);
title('Prior Mean Fit (Full/Reduced Model)');

sse2=mci_plot_viscosity (U,X,Y,post1.Yhat);
title('Posterior Mean Fit (Full Model)');

sse3=mci_plot_viscosity (U,X,Y,post2.Yhat);
title('Posterior Mean Fit (Reduced Model)');

disp(' ');
disp('Full Model:');
disp('Prior Mean  Posterior Mean');
for p=1:length(M.pE);
    disp(sprintf('%1.2f         %1.2f',M.pE(p),post1.Ep(p)));
end

disp(' ');
disp(sprintf('SSE for prior mean = %1.2f',sse1));
disp(sprintf('SSE for posterior mean (full) = %1.2f',sse2));
disp(sprintf('SSE for posterior mean (reduced) = %1.2f',sse3));

if strcmp(mcmc.inference,'vl')
    disp(' ');
    disp(sprintf('VL Log evidence (full) = %1.2f', post1.logev));
    disp(sprintf('VL Log evidence (reduced) = %1.2f', post2.logev));
    disp(sprintf('VL Log BF = %1.2f', post1.logev-post2.logev));
end

if do_ais
    disp('Fitting model using Annealed Importance Sampling');
    lmc.inference='ais';
    lmc.anneal='power';
    lmc.prop='lmc';
    lmc.nprop=1;
    lmc.maxits=32;
    lmc.J=512;
    %lmc.h=0.1;
    
    ais_post1 = spm_mci_post (lmc,M,U,Y);
    ais_post2 = spm_mci_post (lmc,Malt,Ualt,Y);
    
    disp(' ');
    disp(sprintf('AIS Log evidence (full) = %1.2f', ais_post1.logev));
    disp(sprintf('AIS Log evidence (reduced) = %1.2f', ais_post2.logev));
    disp(sprintf('AIS Log BF = %1.2f', ais_post1.logev-ais_post2.logev));

    ais_sse2=mci_plot_viscosity (U,X,Y,ais_post1.Yhat);
    title('AIS Posterior Mean Fit (Full Model)');
    
    ais_sse3=mci_plot_viscosity (U,X,Y,ais_post2.Yhat);
    title('AIS Posterior Mean Fit (Reduced Model)');
    
    disp(' ');
    disp(sprintf('SSE for prior mean = %1.2f',sse1));
    disp(sprintf('SSE for posterior mean (full) = %1.2f',ais_sse2));
    disp(sprintf('SSE for posterior mean (reduced) = %1.2f',ais_sse3));
end
