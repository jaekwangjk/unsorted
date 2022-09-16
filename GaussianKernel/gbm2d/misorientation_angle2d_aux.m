function minang = misorientation_angle2d_aux(ori1,ori2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minang = misorientation_angle2d_aux(ori1,ori2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the misorientation angle between ori1 and ori2.

ang1 = (ori2>ori1).*(2*pi - ori2 + ori1) + (ori2<=ori1).*(2*pi - ori1 + ori2);
minang = min(ang1,abs(ori1-ori2));


end