function [grains,ori] = initialvoronoidata2d(N,dims,Rgrains,th_ang,rn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [grains,ori] = initialvoronoidata2d(N,dims,Rgrains,th_ang)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates initial condition with N grains on a grid of size dims.
% Grains with misorientation angle smaller than th_ang degrees are merged.
% Each grain will have several disconnected components.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(rn); 
label = -1e10 * ones(dims); % very large value initially 
m = dims(1); n = dims(2);
W = int32(-ones(dims));
work_x = int32(zeros(prod(dims),1));
work_y = int32(zeros(prod(dims),1));
work_bx = int32(zeros(prod(dims),1)); % boundary /* List of bdry element indices. */
work_by = int32(zeros(prod(dims),1));
work_candidate_x = int32(zeros(prod(dims),1));
work_candidate_y = int32(zeros(prod(dims),1));

r = 0.1;
go = 1;
while go
    Naux = floor(N*(1+r));
    x = 1 + floor(rand(1,Naux)*(m-1));
    y = 1 + floor(rand(1,Naux)*(n-1));
    
    % We ensure that there are not repeated seeds for the voronoi region
    aux = unique([x;y]','rows','stable');
    x = aux(:,1)';
    y = aux(:,2)';
    Naux = length(x);
    
    if Naux>=N
        go = 0;
        select = randperm(Naux,N);
        x = x(select);
        y = y(select);
    else
        r = r + 0.1;
    end
end

grains_voronoi = cell(N,1);

% Parameter d controls the extent to which the distance funciton
% to each point in the dataset is constructed:
d = 1 + ceil( sqrt(m*n/N) );
d = 2*d;

start_progress(' - Constructing distance function');
% Construct distance function, label, to the union of all the points:
for k=1:N % Loop over the random points.
    [x2,y2] = meshgrid((x(k)-d):(x(k)+d),(y(k)-d):(y(k)+d));
    dist = -sqrt( (x2-x(k)).^2 + (y2-y(k)).^2 );
    x2 = 1 + mod(x2-1,m);
    y2 = 1 + mod(y2-1,n);
    ind = sub2ind(dims,x2(:),y2(:));
    label(ind) = max(dist(:),label(ind));
    display_progress(k,N,1)
end % (for k). Loop over random points ends.

% If the union of d neighborhoods of the random points do not cover
% the entire computational domain, we cannot trust the construction:
if min(label(:)) < -0.5*1e10
    error('Parameter d too small.');
end

start_progress(' - Constructing Voronoi regions');
% Associate each grid point with the random point it is closest to,
% forming the grains:
indgridQ = zeros(dims);
for k=1:N % Loop over the random points again.
    
    [x2,y2] = meshgrid((x(k)-d):(x(k)+d),(y(k)-d):(y(k)+d));
    dist = -sqrt( (x2-x(k)).^2 + (y2-y(k)).^2 );
    x2 = 1 + mod(x2-1,m);
    y2 = 1 + mod(y2-1,n);
    
    ind = sub2ind(dims,x2(:),y2(:));
    ind2 = ind( dist(:) >= label(ind) );
    ind2 = ind2(not(indgridQ(ind2)));
    indgridQ(ind2) = 1;
    grains_voronoi{k,1} = ind2;
    display_progress(k,N,1)
end % (for k). Loop over random points ends.


% Bins with width th_ang are formed for the orientation angles. Grains in
% the same bin are merged together with orientation angle determined by the
% center of the bin.
ori_voronoi = 2*pi*rand(N,1);
bin_edges = 0:th_ang:360;
center_bins = th_ang/2:th_ang:360;
nBins = length(center_bins);
y = zeros(size(ori_voronoi));
for k=1:nBins
    select =  ori_voronoi*360/(2*pi) >= bin_edges(k);
    y(select) = k;
end
% alternatively y can be computed by
% y = discretize(ori_voronoi*360/(2*pi),bin_edges);
% However, only the most recent Matlab versions have the
% function discretize
ori_voronoi = center_bins(y)/360*2*pi;
ori_voronoi = ori_voronoi';

count = 0;
for k = 1:nBins
    ind = find(ori_voronoi == center_bins(k)/360*2*pi);
    if not(isempty(ind))
        count = count+1;
        grains{count,1} = grains_voronoi{ind(1)};
        grains{count,2} = [];
        grains{count,3} = [];
        grains{count,4} = [];
        ori(count,1) = center_bins(k)/360*2*pi;
        if length(ind)>1
            for l = length(ind):-1:2
                grains{count,1} = [grains{count,1};grains_voronoi{ind(l),1}];
            end
        end
    end
end % (for k). Loop over bins ends.

N = size(grains,1);

start_progress(' - Constructing grains');
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