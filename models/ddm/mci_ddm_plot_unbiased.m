
clear all
close all

do_rand = 1;
disp('Unbiased drift diffusion model');
disp('Plot RT density and generate data');
disp('Compare Gaussian with series approx');

% Same params for all trials
[M,U] = mci_ddm_struct();
if do_rand
    P = spm_normrnd(M.pE,M.pC,1);
    ddm = mci_ddm_from_mci (P,M,U);
else
    ddm.v=0.32; ddm.a=0.22; ddm.b=0.47;
    P = mci_ddm_to_mci(ddm,M,U);
end

[mt,vt,pc] = mci_ddm_moments (P,M,U);
disp(sprintf('Mean time %1.2f',mt));
disp(sprintf('SD time %1.2f',sqrt(vt)));
disp(sprintf('Prob correct %1.2f',pc));

t=linspace(ddm.b*1.1,mt+3*sqrt(vt),100);
for i=1:length(t),
    pr(i) = mci_ddm_rt_density (P,M,U,[0 t(i)]);
    pr_nf(i) = mci_ddm_wfpt (P,M,U,[0 t(i)]);
    pr_norm(i) = spm_Npdf (t(i),mt,vt);
end

% Both ddm_rt_density.m and ddm_wfpt.m return 
% unnormalised pdfs. So normalise here:
pr = pr/(1-pc);
pr_nf = pr_nf/(1-pc);

figure
plot(t,pr);
hold on
plot(t,pr_nf,'r');
plot(t,pr_norm,'k');
legend('Tuerlinckx-Wagenmakers','Navarro-Fuss','Gauss');
xlabel('t')
ylabel('p(t)');
grid on
hold on
plot([ddm.b ddm.b],ylim,'r');

% Generate data from PDFs
[y,dist]=mci_ddm_gen(P,M,U,1000); 
figure
plot(dist.t,dist.pdf_correct);
hold on
plot(dist.t,dist.pdf_error,'r');
legend('x=Correct','x=Error');
xlabel('t');
ylabel('p(t|x)');

figure
hist(y(:,2));
dist
smt=mean(y(:,2));
sst=std(y(:,2));
disp(sprintf('Sample Mean RT = %1.2f',smt));
disp(sprintf('Sample SD RT = %1.2f', sst));
disp(sprintf('Sample pc = %1.2f',mean(y(:,1))));
