N=200;
dims=[512 512]; 

energy_curve = importdata("FCC_110STGB.txt");

S = zeros(N,N); 
ori = zeros(N,1); 


%% set orientaiton 
for i=1:length(ori)
    

 if(mod(i,10)==0 || mod(i,10)==1 || mod(i,10)==2 || mod(i,10)==3 )
    ori(i) = 5*rand(1); % from 0 to 5
 elseif(mod(i,10)==4)
    ori(i) = 15 + 5 * rand(1); %from 15 to 20
 elseif(mod(i,10)==5)
    ori(i) = 25 + 5 * rand(1); % from 25 to 30
 elseif(mod(i,10)==6)
    ori(i) = 35 + 5 * rand(1); % from 35 to 40
 elseif(mod(i,10)==7)
    ori(i) = 45 + 5 * rand(1); % from 45 to 50
 elseif(mod(i,10)==8)
    ori(i) = 55 + 5 * rand(1); % from 55 to 60
 elseif(mod(i,10)==9)
    ori(i) = 65 + 5 * rand(1); % from 65 to 70
 end
       
       
end


oriprime = 2*pi*rand(N,1);

ori=oriprime; 



angBrandon = 30; 
angBrandonRad = angBrandon * pi/180; 

%% Alternative construction

N = length(ori);
ori1 = repmat(ori*pi/180,1,N);
ori2 = repmat((ori*pi/180)',N,1);

ang1 = (ori2>ori1).*(2*pi - ori2 + ori1) + (ori2<=ori1).*(2*pi - ori1 + ori2);
minang = min(ang1,abs(ori1-ori2));

Sprime = ones(N,N)-eye(N,N);

% Read-Shockley with Bradon angle %
select = minang<=angBrandonRad;
Sprime(select) = minang(select)/angBrandonRad.*(1-log(minang(select)/angBrandonRad));
Sprime(1:N+1:N^2) = 0;



for i=1:N
    for j=i+1:N
      oriA=ori(i); oriB=ori(j);
      misorient= abs(oriA-oriB);
      
      if(misorient < angBrandon)
      misorientRad = misorient*pi/180; 
      energy = misorientRad/angBrandonRad.*(1-log(misorientRad/angBrandonRad)); 
      else
      energy =1; 
      end
      
      S(i,j) = energy; 
      S(j,i) = S(i,j); 
    end
end

%Mobility Matrix  
M = ones(N)-eye(N);


S=Sprime;

J = eye(N)-1/N*ones(N,1)*ones(1,N); % ortogonal projection 
eigS = eig(J*S*J);
[~,i] = min(abs(eigS));
eigS(i) = [];
mineigS = min(eigS);
maxeigS = max(eigS);

% option 1 
eigreciprocalM = -ones(N-1,1);
maxeigreciprocalM = -1;
mineigreciprocalM = -1;


alpha = mineigS/maxeigreciprocalM;
beta = maxeigS/mineigreciprocalM;

disp(mineigS)
disp(maxeigS)

aux = S.*M;
aux(1:N+1:N^2) = -Inf;
alpha = max(alpha,max(aux(:)));
aux(1:N+1:N^2) = Inf;
beta = min(beta,min(aux(:)));

if sum(eigS > 0)
    error('The matrix \sigma is not conditionally negative semidefinite.')
end
if sum(eigreciprocalM > 0)
    error('The matrix 1/\mu is not conditionally negative semidefinite.')
end





% % for i=1:N
% %     for j=i+1:N
% %       oriA=ori(i); oriB=ori(j);
% %       misorient= abs(oriA-oriB);
% %       energy = interp1(energy_curve(:,1),energy_curve(:,2),misorient);
% %       S(i,j) = energy; 
% %       S(j,i) = S(i,j); 
% %     end
% % end