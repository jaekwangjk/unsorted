
clear all; 
close all; 

data=importdata('U_infty_set_new_mesh.txt'); 


U=data(:,1); 
D1=data(:,2); 
D2=data(:,3); 

data(:,4) = 6*pi* 50 *data(:,1)*0.025; 

d=0.2; % distance between center 
a=0.025; % radius 
ka=0.1; % ka model value 

Processed_data=data; 

for i=1:length(data)
   
    Processed_data(i,2)=Processed_data(i,2)/Processed_data(i,4); 
    Processed_data(i,3)=Processed_data(i,3)/Processed_data(i,4); 
    Processed_data(i,1)=Processed_data(i,1)/(ka*d); 
end


semilogx(Processed_data(:,1),Processed_data(:,2),'ro-','linewidth',2); 
hold on; 
semilogx(Processed_data(:,1),Processed_data(:,3),'bo-','linewidth',2); 

axis([min(Processed_data(:,1)) max(Processed_data(:,1)), -0.1,1.3]); 

%add horizontal line 
hori=[0.0001 0.8441; 10000 0.8441]; 
plot(hori(:,1),hori(:,2),'k--','linewidth',2); 

 
xlabel('$\left(\frac{1}{k_a d} \right)U$','interpreter','latex'); 
ylabel('$\phi$','interpreter','latex'); 
title('Interaction Factor $\phi$','interpreter','latex'); 
set(gca,'fontsize',20,'ticklabelinterpreter','latex')
% 
% % 
% % figure; 
% % loglog(x,data(:,2),'ro-'); hold on; 
% % loglog(x,data(:,3),'bo-'); 

