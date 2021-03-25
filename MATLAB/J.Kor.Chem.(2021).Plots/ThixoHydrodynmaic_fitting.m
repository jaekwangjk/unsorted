clear all 
close all 

data4=importdata("dats/Pe_d2a_4.txt"); 


for i=1:30
   fitdata4(i,1)=log10( data4(i,1)) ; %Pe
   fitdata4(i,2)=data4(i,3);
   
end

%h=figure; 

x1=fitdata4(:,1); 
y=fitdata4(:,2);
x2=data4(:,1); 

y2=1.536./(1+exp(7.131*x2)) + 0.08148; 
y3=-0.7682 * tanh(3.564*x2) + 0.8496 ; 
y4=-0.5095 * atan(6.937 *x2) + 0.8589 ; 


semilogx(data4(:,1),data4(:,3),'ko','markersize',10); hold on ; 
semilogx(x2,y2,'m-.','linewidth',2); 
semilogx(x2,y3,'b--','linewidth',2); 
semilogx(x2,y4,'r-','linewidth',3); 

 
% 
% 
xlabel('$Pe$','fontsize',20,'interpreter','latex','fontsize',20)
ylabel('$\phi$','fontsize',20,'interpreter','latex','fontsize',20)
set(gca,'fontsize',30,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,750,750])

% plot3(plotdata10(:,1),plotdata10(:,2),plotdata10(:,3));