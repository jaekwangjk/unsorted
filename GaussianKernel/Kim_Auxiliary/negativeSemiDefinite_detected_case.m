

CVE = importdata("FCC_110STGB.txt");

x=linspace(0,70.5,100); 
y=interp1(CVE(:,1), CVE(:,2),x); 



N=4;



caseN=4; 

ori = detectedCase(caseN,:); 

%misorienList 
%energyList 



S=zeros(N,N); 

marker =1; 
for i = 1:N
  for j=i+1:N
      
  oriA=ori(i);
  oriB=ori(j);
  
  misorient= abs(oriA-oriB);
  energy = interp1(CVE(:,1),CVE(:,2),misorient);
    
  S(i,j) = energy; 
  S(j,i) = S(i,j); 
  
  
  misorientList(marker,1) = misorient; 
  energyList(marker,1) = energy; 
  marker=marker+1; 
  
  end
    
end



subplot(2,2,caseN)
plot(x,y,'r-'); hold on; 
plot(misorientList,energyList, 'bo', 'markersize',10,'MarkerFaceColor','b' );  

xlabel('Misorientation samples','fontsize',24)
ylabel('Energy','fontsize',24)
set(gca,'fontsize',20,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,800,600])