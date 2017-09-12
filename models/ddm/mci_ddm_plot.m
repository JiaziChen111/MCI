
clear all
close all

do_rand = 1;
disp('Drift diffusion model');
disp('Plot RT densities and generate data');

% Same params for all trials
[M,U] = mci_ddm_struct('biased');
if do_rand
    P = spm_normrnd(M.pE,M.pC,1);
    ddm = mci_ddm_from_mci (P,M,U);
else
    ddm.v=3.66; ddm.a=2.88; ddm.b=0.43; ddm.r=0.58;
    P = mci_ddm_to_mci(ddm,M,U);
end

disp(ddm);

% Generate data from PDFs
[y,dist]=mci_ddm_gen(P,M,U,1000); 

sx(2)=mean(y(:,1));
sx(1)=1-sx(2);
disp(sprintf('Sample pc = %1.4f',sx(2)));
ind{1}=find(y(:,1)==0);
ind{2}=find(y(:,1)==1);
for i=1:2,
    if sx(i)==0
        continue
    end
    ii=ind{i};
    figure
    hist(y(ii,2));
    if i==1
        disp('Error:');
        title('Error');
    else
        disp('Correct:');
        title('correct');
    end
    smt=mean(y(ii,2));
    sst=std(y(ii,2));
    disp(sprintf('Sample Mean RT = %1.2f',smt));
    disp(sprintf('Sample SD RT = %1.2f', sst));
end

figure
plot(dist.t,dist.pdf_correct);
hold on
plot(dist.t,dist.pdf_error,'r');
legend('x=Correct','x=Error');
xlabel('t');
ylabel('p(t|x)');