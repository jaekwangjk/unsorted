
clear all 
close all 

%% Setting up grid point  
L=1.0; 
N=1024; %Number of discretization  (This will make x=0,1,2,3....100 exactly)
h=L/N;

x=zeros(1,N); 
for i=1:N %so that you have N data in total in the grid, 
   x(i)= (i-1)*h;  % 'xb' include the boundary..
end

y=x; 
 
[X,Y] = meshgrid(x,y);

data=importdata("chi_199.txt"); 

fig=figure; 


h=surf(X,Y,data);
set(h,'edgecolor','none')
pbaspect([1 1 1])
view(2)

set(gca,'fontsize' ,20)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(gca,'XColor', 'none','YColor','none')

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,'Kernel','-dpng','-r0')