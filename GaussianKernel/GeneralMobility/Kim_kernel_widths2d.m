function [alpha,beta] = Kim_kernel_widths2d(distinctOri,energy_curve,option)

% option is choice of mobility function 
% option ==1, corresponds to unit mobility. 

N = length(distinctOri); % returns the number of grain

if(option==1)    
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
[~,i] = min(abs(eigS));
eigS(i) = [];



% All mobilities = 1

eigreciprocalM = -ones(N-1,1);
maxeigreciprocalM = -1;
mineigreciprocalM = -1;

% JK correction, force them they all negative. 
absoluteEgienS=-abs(eigS); 

mineigS = min(absoluteEgienS);
maxeigS = max(absoluteEgienS);

alpha = mineigS/maxeigreciprocalM;
beta = maxeigS/mineigreciprocalM;
 
% disp('maxeigS'); disp(maxeigS); 


M = ones(N)-eye(N);
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


end

