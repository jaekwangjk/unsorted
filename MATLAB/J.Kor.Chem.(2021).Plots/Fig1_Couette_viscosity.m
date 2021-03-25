
close all 

eta_str=40;
eta_infty=1.0;
kd=2.0; 
ka=0.1; 
tc=kd/ka; 


x=logspace(-4,4); 
y=eta_infty + eta_str./(tc*x+1); 

h=figure; 

loglog(x,y,'r-','linewidth',2); 



xlabel('$t_c \dot{\gamma}$','fontsize',30,'interpreter','latex')
ylabel('$\eta_{ss}/(\eta_\infty+\eta_{str})$','fontsize',30,'interpreter','latex')
axis([10^(-4) 10^(4) 0.9 100])

text(1,40, '$\eta(\dot\gamma)=\eta_{\infty}+\frac{\eta_{str}}{t_c \dot\gamma+1}$', ...
'interpreter','latex','fontsize',24)

set(gca,'fontsize',30,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,600,450])


set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'Couette_viscosity','-dtiff','-r0')