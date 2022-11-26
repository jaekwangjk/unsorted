clear; 
close all 

N= 10;

CVE = importdata("FCC_110STGB.txt");

% Size of sigma N 

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
    ori(i) = rand * 70.5; 
    %energy(i)=interp1(CVE(:,1),CVE(:,2),ori(i)); 
end

figure; 
%plot(ori,'ro', 'markersize',10,'MarkerFaceColor','r');
 

for i = 1:N
  for j=i+1:N
      
  oriA=ori(i);
  oriB=ori(j);
  
  misorient= abs(oriA-oriB);
  energy = interp1(CVE(:,1),CVE(:,2),misorient);
    
  S(i,j) = energy; 
  S(j,i) = S(i,j); 
  
  if(misorient < angBrandon)
      misorientRad = misorient*pi/180; 
      energy = misorientRad/angBrandonRad.*(1-log(misorientRad/angBrandonRad));
  else
      energy =1; 
  end
      
   SRead(i,j) = energy; 
   SRead(j,i) = SRead(i,j); 
  
  
  
  
  end
    
end

x=ori; 

plot(ori,S(1,:),'ro', 'markersize',10,'MarkerFaceColor','r'); hold on; 
plot(ori,SRead(1,:), 'bo', 'markersize',10,'MarkerFaceColor','b'); 

% 
xlabel('Misorientation','fontsize',24)
ylabel('Energy','fontsize',24)
set(gca,'fontsize',20,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,800,600])
% 



% Construct the basis matrix (the column is the basis) 
A = zeros(N,N); 

for j=1:N
   
    for i=1:N
       
        if(j<2) % first column 
            A(i,j)=1;
        end
        
        if(j>1 && i==j)
            A(i,j)=1; 
        end
    end
    
end


% The colums of Q are the orthogonal vectors
Q=zeros(N,N);
R=zeros(N,N); 

% % Do gram-schmidt algorithm

for j=1:N 
    v=A(:,j); % v begins as j colum of A
    for i=1:j-1
        R(i,j)=Q(:,i)'* A(:,j); % modify(A:,j) to v for more accurcay
        v=v-R(i,j)*Q(:,i); % substract projection 
    end % v is now perpendicular to all q1,...q_j-1
    R(j,j)=norm(v); 
    Q(:,j)=v/R(j,j); % normalized v to be the next unit vector q_j
end

% Transformation matrix 

% T=Q*inv(A); % This is the matrix that T*A = Q


% new Sigma 
Snew = (inv(Q)*S)*Q; 

% and its submatrix 
Snew_sub = Snew(2:N,2:N);
e = eig(Snew_sub);

e_real = real(e);
e_img = imag(e);


figure;
plot(e_real,e_img,'ro','markersize',10,'MarkerFaceColor','r'); hold on; 


% new Sigma for Read
SReadnew = (inv(Q)*SRead)*Q; 

% and its submatrix 
SReadnew_sub = SReadnew(2:N,2:N);
e = eig(SReadnew_sub);

e_real = real(e);
e_img = imag(e);

plot(e_real,e_img,'bo','markersize',10,'MarkerFaceColor','b'); hold on; 




xlabel('Real','fontsize',24)
ylabel('Complex','fontsize',24)
set(gca,'fontsize',20,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,800,600])
