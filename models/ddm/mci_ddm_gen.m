function [y,dist] = mci_ddm_gen (P,M,U,N)
% Generate data from DDM   
% FORMAT [y,dist] = mci_ddm_gen (P,M,U,N)
%
% P         parameters
% M         model structure
% U         input structure (empty)
% N         Number of data points to generate
%
% y         N x 2 matrix with columns [x, t] where x 
%           is binary decision, t is reaction time
% dist      distributions that are sampled from
%
%           .t              RT interval
%           .pdf_correct    pdf over correct RTs
%           .cdf_correct    cdf over correct RTs
%           .pdf_error      pdf over error RTs
%           .cdf_error      cdf over error RTs
%           .rmsq           RMS quantization error 
%           .pc             probability of correct
%
% Sample from series approximation to PDF for DDM.
% Invert CDF using numerical quantization with number
% of bins M.Nq and then linearly interpolate.
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id$

try Nq=M.Nq; catch Nq=1000; end
    
ddm = mci_ddm_from_mci (P,M,U);
v=ddm.v;
a=ddm.a;
b=ddm.b;

tmin = ddm.b*1.01;
if isfield(ddm,'r')
    % Biased DDM
    
    % Compute p(x=1)
    z=ddm.r*a;
    num = exp(-2*z*v/M.s2)-1;
    denom = exp(-2*a*v/M.s2)-1;
    pc = num/denom;
    
    % Find value of t required so that right tail of 
    % p(t) falls to less than pdf_thresh
    pdf_thresh=1e-4;
    tmax = 2*ddm.b;
    p1 = mci_ddm_wfpt (P,M,U,[1 tmax])/pc;
    p0 = mci_ddm_wfpt (P,M,U,[0 tmax])/(1-pc);
    pmax = max([p0,p1]);
    while pmax > pdf_thresh;
        tmax = tmax+ddm.b;    
        p1 = mci_ddm_wfpt (P,M,U,[1 tmax])/pc;
        p0 = mci_ddm_wfpt (P,M,U,[0 tmax])/(1-pc);
        pmax = max([p0,p1]);
    end
else
    % Unbiased DDM
    [mt,vt,pc] = mci_ddm_moments (P,M,U);
    tmax = mt+5*sqrt(vt);
end

% Generate decisions
x = pc > rand(N,1);

% RT CDF quantized over tmin to tmax
dt = (tmax-tmin)/(Nq-1);
t = [tmin:dt:tmax];
for i=1:Nq,
    % get p(x=1,t)
    pdf_correct(i) = mci_ddm_wfpt (P,M,U,[1 t(i)]);
    % get p(x=0,t)
    pdf_error(i) = mci_ddm_wfpt (P,M,U,[0 t(i)]);
end
% get p(t|x=1)
pdf_correct = pdf_correct/pc;
% get p(t|x=0)
pdf_error = pdf_error/(1-pc);
cdf_error = cumsum(pdf_error)*dt;
cdf_correct = cumsum(pdf_correct)*dt;

% Sample from RT density
interpol=1;
p = rand(N,1);
for j=1:N,
    % Distance to quantized cdf values
    if x(j)==0
        d=cdf_error-p(j);
    else
        d=cdf_correct-p(j);
    end
    
    if interpol
        % Interpolate between nearest two bins
        [dd,ind]=sort(d.^2);
        id=1./(dd(1:2));
        w=id/sum(id);
        rt(j)=t(ind(1))*w(1)+t(ind(2))*w(2);
        pq(j)=dd(1);
    else
        % Pick closest
        [pq(j),ind]=min(d.^2);
        rt(j)=t(ind);
    end
        
end

do_plot=0;
if do_plot
    figure
    plot(t,cdf_error,'r');
    hold on
    plot(t,cdf_correct);
    grid on
    xlabel('t');
end

y=[x(:),rt(:)];

dist.t=t;
dist.pdf_correct=pdf_correct;
dist.pdf_error=pdf_error;
dist.cdf_correct=cdf_correct;
dist.cdf_error=cdf_error;
dist.pc=pc;
dist.rmsq=sqrt(mean(pq.^2));