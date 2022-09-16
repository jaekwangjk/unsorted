function S = surface_tension2d(ori,angBradon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = surface_tension2d(ori,angBradon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the surface tension matrix given by the Read-Shockley model

N = length(ori);
angBrandonRad = angBradon*pi/180; % converting to radians
% need to update updatelevelsetdata2d.c accordingly


minang = misorientation_angle2d(ori);

S = ones(N,N)-eye(N,N);

% Read-Shockley with Bradon angle %
select = minang<=angBrandonRad;
S(select) = minang(select)/angBrandonRad.*(1-log(minang(select)/angBrandonRad));
S(1:N+1:N^2) = 0;

end