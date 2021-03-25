clear all; 
close all; 

data=importdata("Moseley_Data.txt"); 

plot(data(:,1),data(:,2),'ro-'); 

xlabel('$z/r$','interpreter','latex'); 
ylabel('$V_F/V_S$','interpreter','latex'); 
title('Moseley Experiment Data','interpreter','latex'); 
set(gca,'fontsize',20,'ticklabelinterpreter','latex')