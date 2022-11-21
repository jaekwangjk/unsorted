clear all 
close all 


%This file load grain data and output statistics. 


%% Isotropic data 

% grains = load("grain_later_iso_120.mat").grains; 
% [x,y] = getAreaStatistics(grains); 
% plot(x,y,'bo-','linewidth',3,'DisplayName','Iso 120 step'); hold on; 
% 
% grains = load("grain_later_iso_180.mat").grains; 
% [x,y] = getAreaStatistics(grains); 
% plot(x,y,'go-','linewidth',2,'DisplayName','Iso 200 step'); hold on; 
% 
% grains = load("grain_later_iso_240.mat").grains; 
% [x,y] = getAreaStatistics(grains); 
% plot(x,y,'ko-','linewidth',2,'DisplayName','Iso 240 step'); hold on; 
% 
% 
grains = load("grainFiles/grain_later_iso_240.mat").grains; 
dims =[2048, 2048]; 

[x,y] = getAreaStatistics(grains,dims); 
plot(x,y,'ro--','linewidth',2,'DisplayName','Iso 240 step'); hold on; 


%% Test data
% grains = load("grain_later_Test2_120.mat").grains; 
% [x,y] = getAreaStatistics(grains); 
% plot(x,y,'bo-','linewidth',3,'DisplayName','Toy 120 step'); hold on; 
% 
% grains = load("grain_later_Test2_180.mat").grains; 
% [x,y] = getAreaStatistics(grains); 
% plot(x,y,'go-','linewidth',2,'DisplayName','Toy 120 step'); hold on; 
% 
% grains = load("grain_later_Test2_240.mat").grains; 
% [x,y] = getAreaStatistics(grains); 
% plot(x,y,'ko-','linewidth',2,'DisplayName','Toy 240 step'); hold on; 

grains = load("grainFiles/grain_later_Test2_240.mat").grains; 
[x,y] = getAreaStatistics(grains,dims); 
plot(x,y,'ro-','linewidth',2,'DisplayName','Toy 240 step'); hold on; 


%% Setting Plot 
legend()
set(gca,'fontsize' ,24)
xlabel('A/A_{avg}','fontsize',24)
ylabel('Probability density','fontsize',24)
