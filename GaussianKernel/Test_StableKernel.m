clear 

N=300; 
ori=zeros(N,1); 

%% reset orientaiton % It should run with this...
for i=1:length(ori)    
 if(mod(i,7)==0 )
 ori(i) = 2 * randi([0,5]); % from 0 to 5
 elseif(mod(i,7)==1)
 ori(i) = 12 +  2 * randi([0,5]); %from 15 to 20
 elseif(mod(i,7)==2)
 ori(i) = 22 +  2 * randi([0,5]); % from 25 to 30
 elseif(mod(i,7)==3)
 ori(i) = 32 +  2 * randi([0,5]); % from 35 to 40
 elseif(mod(i,7)==4)
 ori(i) = 42 +  2 * randi([0,5]); % from 45 to 50
 elseif(mod(i,7)==5)
 ori(i) = 52 +  2 * randi([0,5]); % from 55 to 60
 elseif(mod(i,7)==6)
 ori(i) = 60 +  2 * randi([0,3]); % from 65 to 70
 end
      
end

ori = ori * pi/180 ; % make radian 
distinctOri= unique(ori); 

%% Kernel Generation 


energy_curve= importdata("FCC_110STGB.txt");
option =2; 

N = length(distinctOri); % returns the number of grain

if(option > 0)    
  S = ones(N,N)-eye(N,N); % Surface tension matrix 
  
    for i=1:N
     for j=i+1:N
      % convert to degree...   
      oriA=distinctOri(i)* 180/pi; oriB=distinctOri(j)*180/pi;
      
      misorient= abs(oriA-oriB);
      energy = interp1(energy_curve(:,1),energy_curve(:,2),misorient);
      S(i,j) = energy; 
      S(j,i) = S(i,j); 
     end
    end
  
  % High mobility for high angle grain boundaries 
  % Low mobility for low angle grain boundaries 
  if(option==2)
    iM = zeros(N,N) ; % inverseMobility  
    % Mobility function parameters 
    ThetaMax = 10; 
    Mmax = 1; 
    
    for i=1:N
     for j=i+1:N
      % convert to degree...   
      oriA=distinctOri(i)* 180/pi; oriB=distinctOri(j)*180/pi;
      misorient= abs(oriA-oriB); 
      
       if(misorient < ThetaMax)
        mobility = Mmax * 0.905; 
       else
        mobility = Mmax; 
       end
       
      iM(i,j) = 1.0/mobility;
      iM(j,i) = iM(i,j); 
     end
    end
    
  end    
    
    
elseif(option==-1)
  minang = misorientation_angle2d(distinctOri);
  S = ones(N,N)-eye(N,N); % Surface tension matrix 
  angBrandonRad = 30 * pi/180;  
  select = minang<=angBrandonRad;
  S(select) = minang(select)/angBrandonRad.*(1-log(minang(select)/angBrandonRad));
  S(1:N+1:N^2) = 0;
end


% Read-Shockley with Bradon angle %
% sufrace_tension2d function looks as follow %

J = eye(N)-1/N*ones(N,1)*ones(1,N);
eigS = eig(J*S*J);
[~,i] = min(abs(eigS)); % This get rid of unstable value 
eigS(i) = [];


% Get max and min of reciprocal M 
% All mobilities = 1
if (option==1 || option==-1) 
 eigreciprocalM = -ones(N-1,1);
 maxeigreciprocalM = -1;
 mineigreciprocalM = -1;
elseif(option ==2) 
 eigiM = eig(J*iM*J);
 [~,i] = min(abs(eigiM)); % This get rid of unstable value 
 eigiM(i) = [];
 
 assert( max(eigiM) < 0 ) ; % for numberical stability 
 absoluteEigeniM=-abs(eigiM);  % These are all negative 
 maxeigreciprocalM = max(absoluteEigeniM); 
 mineigreciprocalM = min(absoluteEigeniM); 
end



absoluteEigenS=-abs(eigS);  
mineigS = min(absoluteEigenS);
maxeigS = max(absoluteEigenS);


alpha = mineigS/maxeigreciprocalM;
beta = maxeigS/mineigreciprocalM;

