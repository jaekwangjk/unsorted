clear all 
close all 

data4=importdata("dats/Pe_d2a_4.txt"); 


for i=1:30
  
   plotdata(i,1)=data4(i,1) ; %Pe
   plotdata(i,2)=data4(i,2);
   plotdata(i,3)=data4(i,3);
   
end

figure; 

semilogx( plotdata(:,1),plotdata(:,2),'bo-','linewidth',2); hold on 
semilogx( plotdata(:,1),plotdata(:,3),'ro-','linewidth',2); 
axis([10^(-3) 10^2 -0.1 1.2])

legend('Front Sphere','Back sphere')

xlabel('Pe','fontsize',20,'interpreter','latex','fontsize',20)
ylabel('Interaction Factor $\phi$','fontsize',20,'interpreter','latex','fontsize',20)
set(gca,'fontsize',25,'ticklabelinterpreter','latex')
%set(gcf,'units','points','position',[30,30,750,750])

