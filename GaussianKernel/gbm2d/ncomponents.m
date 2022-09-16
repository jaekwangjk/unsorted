function ncomponents = ncomponents(dims,ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ncomponents = ncomponents(dims,ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the number of components in the set with linear indices ind.

[I,J] = ind2sub(dims,ind);

still_left_boundary = sparse(I,J,1,dims(1),dims(2),length(ind));

[indN,indNE,indE,indSE,indS,indSW,indW,indNW] = find_ind_neighbors2d(I,J,dims);

boundary = not(still_left_boundary(indN))+not(still_left_boundary(indNE))+not(still_left_boundary(indE))+not(still_left_boundary(indSE))+not(still_left_boundary(indS))+not(still_left_boundary(indSW))+not(still_left_boundary(indW))+not(still_left_boundary(indNW));

Ibd = I(boundary>0);
Jbd = J(boundary>0);

boundary = sparse(Ibd,Jbd,1,dims(1),dims(2),length(Ibd));

still_left_boundary = sparse(Ibd,Jbd,1,dims(1),dims(2),length(Ibd));

nleft = length(Ibd);
ncomponents = 0;
while nleft > 0
    ncomponents = ncomponents+1;
    component = false(dims);
    ind_component = find(still_left_boundary==1,1);
    still_left_boundary(ind_component) = 0;
    component(ind_component) = true;
    change = 1;
    while true
        [I,J] = ind2sub(dims,ind_component);
        [indN,indNE,indE,indSE,indS,indSW,indW,indNW] = find_ind_neighbors2d(I,J,dims);
        ind_neighbors = [indN;indNE;indE;indSE;indS;indSW;indW;indNW];
        ind_component = ind_neighbors(boundary(ind_neighbors)>0);
        ind_component = unique(ind_component(not(component(ind_component))));
        component(ind_component) = true;
        change = change+length(ind_component);
        still_left_boundary(ind_component) = 0;
        if change == 0
            break
        else
            nleft = nleft-change;
            change = 0;
        end
    end
end

end