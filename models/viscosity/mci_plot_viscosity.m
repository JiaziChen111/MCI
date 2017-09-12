function [sse] = mci_plot_viscosity (U,X,Y,yhat)
% Plot viscosity versus pressure at four temperatures
% FORMAT [sse] = mci_plot_viscosity (U,X,Y,yhat)
%
% See e.g. Figure 3.4 in [1]

figure
sym={'kd','ko','k*','kx'};
sse=0;
for i=1:4,
    j=U.i{i};
    hold on
    plot(X(j,2),Y(j),sym{i},'MarkerSize',10);
    plot(X(j,2),yhat(j),'r-','LineWidth',2);
    sse=sse+sum((Y(j)-yhat(j)).^2);
end
set(gca,'FontSize',16);
xlabel('Pressure');
ylabel('Log Viscosity');
grid on
ylim([2 14]);