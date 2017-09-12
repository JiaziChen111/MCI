function [mlm,post] = spm_mci_mlm (mlm,Y,S)
% Blocked Gibbs sampling for Multivariate Linear Model
% FORMAT [mlm,post] = spm_mci_mlm (mlm,Y,S)
%
% mlm       Multivariate Linear Model
%           .infer_Lambda=0 to use known precision Lambda
% Y         Data
% S         Number of Samples
%
% mlm       Data structure
% post      Samples 
%           .m(:,s) and .Lambda(:,:,s)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

Lambda=mlm.Lambda;

for s=1:S,
    
    % Group Effects
    Lambda_tilde = kron(eye(mlm.N),Lambda);
    ml = mlm.U'*Lambda_tilde;
    PhiN = ml*mlm.U+mlm.Psi0;
    iPhiN = inv(PhiN);
    mN = iPhiN *(ml*Y+mlm.Psi0*mlm.m0);
    m = spm_normrnd(mN,iPhiN,1);
    
    if mlm.infer_Lambda
        % Group Precision
        z = Y - mlm.U*m;
        Z = reshape(z,mlm.K,mlm.N);
        Cz = cov(Z');
        
        aN = mlm.a0+mlm.N/2;
        BN = mlm.B0+Cz*mlm.N/2;
        Lambda = spm_wishrnd(BN,aN);
    end
    
    post.m(:,s)=m;
    post.Lambda(:,:,s)=Lambda;
end

