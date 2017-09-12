function [log_prob,S,E] = mci_plot_surface_slices (M,U,Y,S,dist)
% Plot log probability surface for trivariate parameter space
% FORMAT [log_prob,S,E] = mci_plot_surface_slices (M,U,Y,S,dist)
%
% M         Model structure
% U         Inputs
% Y         Data
% S         Surface data structure
% dist      'prior', 'like', 'post'
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$


% Number of bins defining surface in each dimension
Nbins=S.Nbins;
pxy=S.pxy;

[bup,bdown] = meshgrid(pxy(1,:),pxy(2,:));
        
% For computing log prior term
M = spm_mci_priors (M);
M = spm_mci_minit (M);

% Which params to vary
P1str=[S.param{1},'=P1;'];
P2str=[S.param{2},'=P2;'];
P3str=[S.param{3},'=P3;'];

[I,J]=size(bup);
K=length(S.slice);
rk=ceil(sqrt(K));

h=figure;
set(h,'Name',['log ',dist]);

lpmin=Inf;
lpmax=-Inf;

for k=1:K,
    P3=S.slice(k);
    eval(P3str);
    for i=1:I,
        for j=1:J,
            P1=bup(i,j);
            P2=bdown(i,j);
            eval(P1str);
            eval(P2str);
            
            % Get parameters in eigenspace of prior
            Pv = spm_vec(P);
            M.vpE=spm_vec(M.pE);
            p = M.V'*(Pv-M.vpE);
            switch lower(dist),
                case 'post',
                    log_prob{k}(i,j)= spm_mci_joint (p,M,U,Y);
                case 'prior',
                    log_prob{k}(i,j)=-p'*M.ipC*p/2 + M.log_prior_t2;
                case 'like',
                    log_prob{k}(i,j)= feval(M.L,Pv,M,U,Y);
                otherwise
                    disp('Unknown distribution type');
                    return
            end
        end
    end
    % Work out common scale for all plots
    lpmin=min([min(min(log_prob{k})),lpmin]);
    lpmax=max([max(max(log_prob{k})),lpmax]);
end

for k=1:K,
    subplot(rk,rk,k);
    imagesc(S.x1val(:),S.x2val(:),log_prob{k},[lpmin lpmax]);
    set(gca,'FontSize',18);
    axis xy
    hold on
    xlabel(S.name{1});
    ylabel(S.name{2});
    title(sprintf('%s = %1.2f',S.name{3},S.x3val(k)));
end


