close all 
clear all 

%Stokes Correction for 

nP=40; 
alpha_vlist=linspace(0.5,3.0,nP); 
phi=zeros(nP,1); 


for i=1:nP
    
    alpha = alpha_vlist(i); 
    
    for n=1:20

        top = 4 * ( sinh (n*alpha + 0.5*alpha ) )^2  - (2*n+1)^2 * ( sinh(alpha) )^2;          
        bottom = 2*sinh( 2*n*alpha+ alpha)   +  (2*n+1) * ( sinh(2*alpha) ); 

        value = top/ bottom; 
        
        
        phi(i) = phi(i) + (n^2+n)/(2*n-1)/(2*n+3)  *...
                (1- value);
    end

    % Here, the front fractor 4/3 is correct by Brenner, 1961
    phi(i) = (4/3) * sinh(alpha) * phi(i); 
end



x=cosh(alpha_vlist); 

h=figure; 

plot(real(x),phi,'r-','linewidth',2);
axis([0.5 11 0.6 1]); 

ylabel('$\phi(d/2a)$','fontsize',24,'interpreter','latex');
xlabel('$d/2a$','fontsize',24,'interpreter','latex'); 
set(gca,'fontsize',25,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,600,450])

hold on

x0=4; y0=0.8441; 
plot(x0,y0,'bd','markerfacecolor','blue','markersize',10); 

% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'Newtonian_Interaction_Factor','-dtif','-r0')
% 

zz(:,1)=real(x);
zz(:,2)=phi;
