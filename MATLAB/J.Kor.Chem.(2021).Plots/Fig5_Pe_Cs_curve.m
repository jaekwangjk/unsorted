
clear all; 
close all; 

data=importdata('Pe_d2a_4.txt'); 



D1=data(:,2); 
D2=data(:,3); 


%%% Separate analysis 
fig=figure; 

a=0.025;
K=1.9789; 
kd=2; ka=0.1; 

d= 8 * a; 
Pe=data(:,1); 

semilogx(Pe,D1,'ro-','linewidth',2); hold on; 
semilogx(Pe,D2,'bo-','linewidth',2); 
%axis([min(x) max(x) 0.0 1.2]); 

axis([10^(-3) 10^2 -0.1 1.2])

legend('Front Sphere','Back sphere')

xlabel('Pe','fontsize',20,'interpreter','latex','fontsize',20)
ylabel('Drag Coefficient Cs','fontsize',20,'interpreter','latex','fontsize',20)
set(gca,'fontsize',23,'ticklabelinterpreter','latex')



set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,'Figure 5','-dtiff','-r0')


