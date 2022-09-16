function [indN,indNE,indE,indSE,indS,indSW,indW,indNW] = find_ind_neighbors2d(I,J,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [indN,indNE,indE,indSE,indS,indSW,indW,indNW] = find_ind_neighbors2d(I,J,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns the linear indices of the 8 neighbours of (I,J). I and J can be
% vectors.
m = dims(1); n = dims(2);

Iplus = I+1;
Iplus(Iplus>m) = 1;
Iminus = I-1;
Iminus(Iminus<1) = m;
Jplus = J+1;
Jplus(Jplus>n) = 1;
Jminus = J-1;
Jminus(Jminus<1) = n;


indN = sub2ind(dims,Iminus,J);
indNE = sub2ind(dims,Iminus,Jplus);
indE = sub2ind(dims,I,Jplus);
indSE = sub2ind(dims,Iplus,Jplus);
indS = sub2ind(dims,Iplus,J);
indSW = sub2ind(dims,Iplus,Jminus);
indW = sub2ind(dims,I,Jminus);
indNW = sub2ind(dims,Iminus,Jminus);

end