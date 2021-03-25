clear all 
close all 

data2=importdata("dats/Pe_d2a_2.txt"); 
data4=importdata("dats/Pe_d2a_4.txt"); 
data6=importdata("dats/Pe_d2a_6.txt"); 
data8=importdata("dats/Pe_d2a_8.txt"); 
data10=importdata("dats/Pe_d2a_10.txt");

data3=importdata("dats/Pe_d2a_3.txt"); 
data5=importdata("dats/Pe_d2a_5.txt"); 
data7=importdata("dats/Pe_d2a_7.txt"); 
data9=importdata("dats/Pe_d2a_9.txt"); 


ka=0.1; 
a=0.025; 

for i=1:30
   plotdata2(i,1)=2; %d/2a
   plotdata3(i,1)=3; 
   plotdata4(i,1)=4; 
   plotdata5(i,1)=5; 
   plotdata6(i,1)=6; 
   plotdata7(i,1)=7; 
   plotdata8(i,1)=8; 
   plotdata9(i,1)=9; 
   plotdata10(i,1)=10; 
   
   plotdata2(i,2)=log10( data2(i,1)*0.1*0.2 / ka/ (2*a *2 ) )  ; %Pe
   plotdata2(i,3)=data2(i,3);
   plotdata2(i,4)=data2(i,2);
   
   plotdata3(i,2)=log10( data3(i,1)*0.1*0.2 / ka/ (2*a *3 )) ; %Pe
   plotdata3(i,3)=data3(i,3);
   plotdata3(i,4)=data3(i,2);
   
   plotdata4(i,2)=log10( data4(i,1)*0.1*0.2 / ka/ (2*a *4 )) ; %Pe
   plotdata4(i,3)=data4(i,3);
   plotdata4(i,4)=data4(i,2);
   
   plotdata5(i,2)=log10( data5(i,1)*0.1*0.2 / ka/ (2*a *5 )) ; %Pe
   plotdata5(i,3)=data5(i,3);
   plotdata5(i,4)=data5(i,2);
   
   plotdata6(i,2)=log10( data6(i,1)*0.1*0.2 / ka/ (2*a *6 )) ; %Pe
   plotdata6(i,3)=data6(i,3);
   plotdata6(i,4)=data6(i,2);
   
   plotdata7(i,2)=log10( data7(i,1)*0.1*0.2 / ka/ (2*a *7 )) ; %Pe
   plotdata7(i,3)=data7(i,3);
   plotdata7(i,4)=data7(i,2);
   
   plotdata8(i,2)=log10( data8(i,1)*0.1*0.2 / ka/ (2*a *8 )) ; %Pe
   plotdata8(i,3)=data8(i,3);
   plotdata8(i,4)=data8(i,2);
   
   plotdata9(i,2)=log10( data9(i,1)*0.1*0.2 / ka/ (2*a *9 )) ; %Pe
   plotdata9(i,3)=data9(i,3);
   plotdata9(i,4)=data9(i,2);
   
   plotdata10(i,2)=log10( data10(i,1)*0.1*0.2 / ka/ (2*a *10)) ; %Pe
   plotdata10(i,3)=data10(i,3);
   plotdata10(i,4)=data10(i,2);
   
end


x=[2,3,4,5,6,7,8,9,10]; 
y=plotdata2(:,2); 
[X,Y]=meshgrid(x,y); 

Z=zeros(length(y),length(x)); 
Zf=zeros(length(y),length(x)); 

for i=1:length(y)
   Z(i,1) = plotdata2(i,3); 
   Z(i,2) = plotdata3(i,3); 
   Z(i,3) = plotdata4(i,3); 
   Z(i,4) = plotdata5(i,3); 
   Z(i,5) = plotdata6(i,3); 
   Z(i,6) = plotdata7(i,3);
   Z(i,7) = plotdata8(i,3);
   Z(i,8) = plotdata9(i,3);
   Z(i,9) = plotdata10(i,3);
   
   
   Zf(i,1) = plotdata2(i,4); 
   Zf(i,2) = plotdata3(i,4); 
   Zf(i,3) = plotdata4(i,4); 
   Zf(i,4) = plotdata5(i,4); 
   Zf(i,5) = plotdata6(i,4); 
   Zf(i,6) = plotdata7(i,4);
   Zf(i,7) = plotdata8(i,4);
   Zf(i,8) = plotdata9(i,4);
   Zf(i,9) = plotdata10(i,4);
end



%%

figure; 
[M,c]=contourf(X,Y,Z,[0.01,0.1,0.2,0.3,0.5,0.7,0.8,0.85,0.87,0.88,0.9,0.91],'Showtext','on'); 


colormap parula
% % 
c.LineWidth=1.5;
clabel(M,c,'FontSize',18)

title('Back Sphere','fontsize',20);
xlabel('$d/2a$','fontsize',20,'interpreter','latex','fontsize',20)
ylabel('$\log{\mathrm{Pe}}$','fontsize',20,'interpreter','latex','fontsize',20)
set(gca,'fontsize',30,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,750,750])

%%

figure; 
[M,c]=contourf(X,Y,Zf,[0.01,0.1,0.2,0.3,0.5,0.7,0.8,0.85,0.87,0.88,0.9,0.91],'Showtext','on'); 


colormap parula
% % 
c.LineWidth=1.5;
clabel(M,c,'FontSize',18)

title('Front Sphere','fontsize',20);
xlabel('$d/2a$','fontsize',20,'interpreter','latex','fontsize',20)
ylabel('$\log{\mathrm{Pe}}$','fontsize',20,'interpreter','latex','fontsize',20)
set(gca,'fontsize',30,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,750,750])



% %%
figure; 
[M,c]=contourf(X,Y,abs(Zf-Z),[0, 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1],'Showtext','on' ); 


colormap parula
% % 
c.LineWidth=1.5;
clabel(M,c,'FontSize',18)

%surf(X,Y,Zf-Z); 
title('Difference','fontsize',20);
xlabel('$d/2a$','fontsize',20,'interpreter','latex','fontsize',20)
ylabel('$\log{\mathrm{Pe}}$','fontsize',20,'interpreter','latex','fontsize',20)
set(gca,'fontsize',30,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,750,750])
