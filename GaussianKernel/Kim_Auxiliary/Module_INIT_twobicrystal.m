function [grains,ori] = Module_INIT_twobicrystal(N,dims,Rgrains)
 
assert(N == 3); 


label = -1e10 * ones(dims); % very large value initially 
m = dims(1); n = dims(2);
W = int32(-ones(dims));
work_x = int32(zeros(prod(dims),1));
work_y = int32(zeros(prod(dims),1));
work_bx = int32(zeros(prod(dims),1)); % boundary /* List of bdry element indices. */
work_by = int32(zeros(prod(dims),1));
work_candidate_x = int32(zeros(prod(dims),1));
work_candidate_y = int32(zeros(prod(dims),1));

x=zeros(2,1); 
y=zeros(2,1); 

x(1) = m/4; x(2) = 3*m/4; 
y(1) = n/2; y(2) = n/2; 

grains_voronoi = cell(N,1);

[x2,y2] = meshgrid(1:m,1:n);

dist2 = sqrt( (x2-x(1)).^2 + (y2-y(1)).^2 ); 
dist3 = sqrt( (x2-x(2)).^2 + (y2-y(2)).^2 );
        
ind = sub2ind(dims,x2(:),y2(:));
ind2 = ind(dist2(:)< m/25);
ind3 = ind(dist3(:)< m/25); 

ind1 = setdiff(ind,ind2);
ind1 = setdiff(ind1,ind3);

grains_voronoi{1,1} = ind1;
grains_voronoi{2,1} = ind2;
grains_voronoi{3,1} = ind3;



count =0; 
for k=1:N
    count=count+1; 
    grains{k,1} = grains_voronoi{k};
    grains{count,2} = [];
    grains{count,3} = [];
    grains{count,4} = [];
end

% Dilate the grain neighborhood, and define a level set function
% on it (1 inside, -1 outside):
for k=1:N % Loop over the grains.
    ind = grains{k,1};
    [x,y] = ind2sub(dims,ind); % dims is the grid matrix, and ind is the x grid in the 1d darry 
    [x2,y2] = pgrow3(int32(x),int32(y),Rgrains,W,...
        work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y);
    ind2 = sub2ind(dims,x2,y2);
    label(ind2) = -1; % outside
    label(ind) = 1; % inside 
    grains{k,1} = ind2;
    grains{k,2} = label(ind2);
    grains{k,3} = 0*ind2;       % Convolution vals. of alpha kernel init to 0.
    grains{k,4} = 0*ind2;       % Convolution vals. of beta kernel init to 0.
	display_progress(k,N,1);
end % (for k). Loop over grains ends.



ori=zeros(N,1);