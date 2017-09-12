
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

disp('Setting viscosity priors so that predictions are');
disp('in reasonable range');

% Run Annealed Importance Sampling as well as VL ?
do_ais=1;

X=viscosity_data_ord(:,2:3);
X(:,2)=X(:,2)/1000;

[M,U] = mci_viscosity_struct (X);

y_thresh=20;

N=1000;
k=1;m=1;
w=spm_normrnd(M.pE,M.pC,N);
for n=1:N,
    yhat = mci_viscosity_gen (w(:,n),M,U);
    ymax(n)=max(yhat);
    if ymax(n) < y_thresh
        % Keep the good samples
        theta(:,k)=w(:,n);
        k=k+1;
    else
        theta_bad(:,m)=w(:,n);
        m=m+1;
    end
end
Ngood=size(theta,2);
disp(sprintf('%d good samples out of %d',Ngood,N));

stable.pE=mean(theta')';
stable.pC=cov(theta');
save viscosity_stable_priors stable

figure
for p=1:9,
    sbad=length(find(theta_bad(p,:)<0))/m;
    sgood=length(find(theta(p,:)<0))/k;
    disp(sprintf('Var %d, bad(-ve)=%1.2f, good(-ve)=%1.2f',p,sbad,sgood));
    subplot(3,3,p);
    [nb,xb]=hist(theta_bad(p,:));
    bar(xb,nb,'r');
    hold on
    [ng,xg]=hist(theta(p,:));
    bar(xg,ng,'b');
    title(sprintf('w(%d)',p));
end

return
% Check - are all samples from prior good ?

N=1000;
k=1;
w=spm_normrnd(stable.pE,stable.pC,N);
for n=1:N,
    yhat = mci_viscosity_gen (w(:,n),M,U);
    ymax(n)=max(yhat);
    if ymax(n) < y_thresh
        % Keep the good samples
        theta_w(:,k)=w(:,n);
        k=k+1;
    end
end
Ngood=size(theta_w,2);
disp(sprintf('New priors: %d good samples out of %d',Ngood,N));