function [beta] = spm_mci_schedule (anneal,J)
% Create annealing schedule
% FORMAT [beta] = spm_mci_schedule (anneal,J)
%
% anneal        type of schedule:
%               'sigmoid', 'linear', 'nonlinear', 'log', 'power', 'frozen'
% J             Number of discretisations (temperatures)
%
% beta          vector of inverse temperatures
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id$

% Set annealing schedule
switch anneal
    case 'sigmoid',
        % Good for model switching integration
        x=linspace(-10,10,J);
        beta=1./(1+exp(-x));
    case 'linear',
        beta=linspace(0,1,J);
    case 'geometric',
        beta=logspace(log10(1/J),0,J-1);
        beta=[0 beta];
    case 'nonlinear',
        eta=0.2;
        j=1:J;
        beta=(eta*j/J)./(1-j/J+eta);
    case 'power',
        % From Calderhead & Girolami
        beta=linspace(0,1,J);
        beta=beta.^5;
    case 'power2',
        beta=linspace(0,1,J);
        beta=beta.^2;
    case 'power3',
        beta=linspace(0,1,J);
        beta=beta.^3;
    case 'power8',
        beta=linspace(0,1,J);
        beta=beta.^8;
    case 'power-rev',
        % As power but reversed ie. start at low temp
        beta=linspace(0,1,J);
        beta=beta.^5;
        beta=fliplr(beta);
    case 'frozen',
        beta=ones(1,J);
    otherwise
        disp('Unknown type of annealing schedule in spm_mci_schedule.m');
        return
end