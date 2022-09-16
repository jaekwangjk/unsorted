function minang = misorientation_angle2d(ori)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minang = misorientation_angle2d(ori)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the misorientation angle matrix.

N = length(ori);
ori1 = repmat(ori,1,N);
ori2 = repmat(ori',N,1);

ang1 = (ori2>ori1).*(2*pi - ori2 + ori1) + (ori2<=ori1).*(2*pi - ori1 + ori2);
minang = min(ang1,abs(ori1-ori2));


end