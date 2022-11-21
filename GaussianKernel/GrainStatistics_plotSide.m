%% Data values
clear 

close all 


sideBin = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5]; 
x=zeros(length(sideBin)-1,1);

for k=1:length(sideBin) -1
  x(k) = (sideBin(k)+sideBin(k+1))* 0.5; 
end

y_iso = [0,0.0010,0.0088,0.0867,0.2729,0.3207,0.1940,0.0799,0.0283,0.0068,0.0010,0]; 

y_poly =[0, 0.0010,0.0152,0.1013,0.2918,0.2634,0.1854,0.0892,0.0426,0.0101,0,0 ];


h=figure; 


plot(x,y_iso,'ro-','linewidth',2,'DisplayName','Isotropic'); hold on; 
plot(x,y_poly,'bo-','linewidth',2,'DisplayName','Anisotropy'); hold on; 





%% Setting Plot 
legend()
xlabel('Side','fontsize',15)
ylabel('Probability density','fontsize',15)

%axis([0,5,0,0.35]); 

set(gca,'fontsize',25,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,800,600])


set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'Poly3','-dpdf','-r0')


