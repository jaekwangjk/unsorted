clear; 
close all 

N= 13;
CVE = importdata("FCC_110STGB.txt");

testNumber = 1000; 
% Size of sigma N 
largeEigenvalues = zeros(testNumber,2) ; 


countPositive =0; 

for t =1:testNumber    
S=zeros(N,N); 

% for compasrion add read schockely energy 
SRead=zeros(N,N); 
angBrandon = 30; 
angBrandonRad = angBrandon * pi/180; 

% construct S(i,j) such that j>i and symmetrize them. 
% construct random (N-1) orientation 
ori =zeros(N,1); 
energy = zeros(N,1); 


for i = 2:N
    ori(i) = rand * 65; 
end


for i = 1:N
  for j=i+1:N
      
  oriA=ori(i);
  oriB=ori(j);
  
  misorient= abs(oriA-oriB);
  energy = interp1(CVE(:,1),CVE(:,2),misorient);
    
  S(i,j) = energy; 
  S(j,i) = S(i,j); 
  
  end
end

%% Eigen Value analysis 
J = eye(N)-1/N*ones(N,1)*ones(1,N);
eigS = eig(J*S*J);

e_real = real(eigS); 
e_img = imag(eigS); 


sort_ereal = sort(e_real); 

% save all the largest eigven value 
largeEigenvalues(t,1) = max(e_real); 
largeEigenvalues(t,2) = sort_ereal(N-1);


if( max(e_real) > 0.000001) 
   countPositive =countPositive +1;
   detectedEig(countPositive) =  largeEigenvalues(t,1);
   for q=1:N
   detectedCase(countPositive,q) = ori(q);
   end
end


% for q=1:N
% allcase(k,q)=ori(q); 
% end

end



Probability = countPositive/testNumber; 


%plot(e_real,e_img,'bo','markersize',10,'MarkerFaceColor','b'); hold on; 
%xlabel('Real','fontsize',24)
%ylabel('Complex','fontsize',24)
%set(gca,'fontsize',20,'ticklabelinterpreter','latex')
%set(gcf,'units','points','position',[30,30,800,600])


