function [prop] = spm_mci_compute_proposal (prop,M,U,Y,beta)
% Compute proposal mean and covariance
% FORMAT [prop] = spm_mci_compute_proposal (prop,M,U,Y,beta)
%
% prop      Proposal containing 
%           .h step size, 
%           .pos the prev value
%           .update 'default' or 'alt'
% M         Model (struct for single model, cell for two models)
% U         Inputs
% Y         Data
% beta      temperature
% 
% prop      Proposal Data Structure
%           .Cp     covariance
%           .dpos   mean = pos+dpos 

pos=prop.pos;
h=prop.h;

old_code=1;

if ~iscell(M)
    % Single Model
    M.beta=beta;
    [j,iCpY,st,prop.L,prop.L2] = spm_mci_joint_grad (pos,M,U,Y);
    if st==-1
        disp('Integration problem in spm_mci_lgv.m');
        keyboard
    end
    % Posterior covariance under local linear approximation
    if strcmp(prop.update,'alt')
        prop.iCp = (iCpY+M.ipC)/h; % Old code
    else
        prop.iCp = (iCpY+M.ipC)/(h^2); % As in paper
    end
    prop.Cp = inv(prop.iCp);
    prop.logdetCp = -spm_logdet(prop.iCp);
    if strcmp(prop.update,'alt')
        prop.dpos = 0.5*h*prop.Cp*j(:); % Old code
    else
        prop.dpos = 0.5*prop.Cp*j(:); % As in paper
    end
else
    % Two models
    M1=M{1}; M2=M{2};
    M1.beta=1; M2.beta=1;
    [jF,iCpYF,st,LogJoint1] = spm_mci_joint_grad (pos,M1,U,Y);
    jR = spm_mci_joint_grad (pos,M2,U,Y);
    if st==-1
        disp('Integration problem in spm_mci_lgv.m');
        keyboard
    end
    
    % Compute gradient and curvature gC and FC
    gC = (1-beta)*jF+beta*jR;
    FC = iCpYF;
    
    LogPrior1 = spm_mci_log_prior (pos,M1);
    LogPrior2 = spm_mci_log_prior (pos,M2);
    Like = LogJoint1-LogPrior1; % Log likelihood of params (under either model)
    LogJoint2 = Like + LogPrior2; % Joint for model 2
    
    prop.L = (1-beta)*LogJoint1+beta*LogJoint2;
    prop.Lj2 = LogJoint2;
    prop.Ldiff = LogPrior2-LogPrior1;
    
    % Posterior covariance under local linear approximation
    prop.iCp = (FC+(1-beta)*M1.ipC+beta*M2.ipC)/h;
    prop.Cp = inv(prop.iCp);
    prop.logdetCp = -spm_logdet(prop.iCp);
    prop.dpos = 0.5*h*prop.Cp*gC(:);
end