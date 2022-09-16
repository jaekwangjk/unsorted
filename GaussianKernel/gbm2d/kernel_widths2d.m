function [alpha,beta] = kernel_widths2d(ori,angBrandon,option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [alpha,beta] = kernel_widths2d(ori,angBrandon,option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the kernel widths to use in the Algorithm 2 from
%   Salvador, T. & Esedoglu, S. J Sci Comput (2019) 79: 648.
%   https://doi.org/10.1007/s10915-018-0866-8   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(ori,1); % returns the number of grain
S = surface_tension2d(ori,angBrandon);

%
% sufrace_tension2d function looks as follow 
%
% minang = misorientation_angle2d(ori);
% S = ones(N,N)-eye(N,N); Surface tension matrix 
% Read-Shockley with Bradon angle %
%select = minang<=angBrandonRad;
%S(select) = minang(select)/angBrandonRad.*(1-log(minang(select)/angBrandonRad));
% S(1:N+1:N^2) = 0;
%
%

switch option
    case 1
        % All mobilities = 1
        M = ones(N)-eye(N);
    case 2
        % All mobilities = 1/(surface tension)
        M = 1./S;
        M(1:(N+1):end) = 0;
end

J = eye(N)-1/N*ones(N,1)*ones(1,N);

eigS = eig(J*S*J);
[~,i] = min(abs(eigS));
eigS(i) = [];

mineigS = min(eigS);
maxeigS = max(eigS);

if option == 1
    eigreciprocalM = -ones(N-1,1);
    maxeigreciprocalM = -1;
    mineigreciprocalM = -1;
else
    reciprocalM = M.^-1;
    reciprocalM((1:N+1:end)) = 0;
    eigreciprocalM = eig(J*reciprocalM*J);
    [~,i] = min(abs(eigreciprocalM));
    eigreciprocalM(i) = [];
    maxeigreciprocalM = max(eigreciprocalM);
    mineigreciprocalM = min(eigreciprocalM);
end

alpha = mineigS/maxeigreciprocalM;
beta = maxeigS/mineigreciprocalM;

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

% In the special case where the surface tensions and the reciprocal
% mobilities are the same, alpha and beta will be the same. In such event,
% using the algorithm with two Gaussian kernels is not
% recommended as one Gaussian kernel is indeed enough. Nonetheless the
% algorithm works provided we chose \alpha > \beta.
if alpha == beta
    alpha = 1.1*alpha;
    beta = 0.9*beta;
end
end

