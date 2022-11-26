function Module_DataProcess_showgrainswithori2d(grains,dims,ori,id)
% showgrainswithori2d(grains,dims,ori,id)
% Last input variable, "id", is optional.
% The visualization color is so bad...

% id is just optional..it 
if nargin < 4
  id = (1:1:size(grains,1))';
end


N = size(grains,1); % Number of grains.
u1 = zeros(dims);
u2 = zeros(dims);

for k=1:N % Loop over grains.
  ind = grains{k,1}; % Pixels in a nhd. of the grain.
  val = grains{k,2}; % Level set values. 
  ind2 = ind(val>0); % Pixels in the interior of grain.
  ang = ori(id(k)); % read orientatio values 
  u1(ind2) = (1+cos(ang))/2;
  u2(ind2) = (1+sin(ang))/2;
end

% Colormap does not work because U is provided with RGB 
% You need to figure out nice combination of RGB. 
U(:,:,1) = u1; % R
U(:,:,2) = u2; % G
U(:,:,3) = ones(dims); %B
clf % Clear the figure 
image(U);
axis square;
%colormap gray
colorbar 
end