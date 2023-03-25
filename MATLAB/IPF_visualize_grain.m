
% This code visualizes orientation distribution maps for 2D grains,
% of which their [110] planes are aligned with the z-axis of the reference
% frame
% I refered the post 
% https://mathematica.stackexchange.com/questions/47492/how-to-create-an-inverse-pole-figure-color-map

% Load polycrystal file 
grains = load("TestGrain.mat").grains;
dims = load("TestGrain.mat").dims;
ori = load("TestGrain.mat").ori;
ori = ori * pi/180; % it should have been a radians....but it is okay...

Ng = size(grains,1); % Number of grains.
direction =zeros(length(ori),3); 
rgb =zeros(length(ori),3); 

% RGB colorkeys 
u1 = zeros(dims);
u2 = zeros(dims);
u3 = zeros(dims); 

for k=1:Ng % Loop over grains.
  ind = grains{k,1}; % Pixels in a nhd. of the grain.
  val = grains{k,2}; % Level set values. 
  ind2 = ind(val>0); % Pixels in the interior of grain.

  RotGrain = determineRotationField_rad(ori(k)); % Rotation Field 

  % What is the meaning of RotGrain?
  % RotGrain makes 110 crystal planes of sample grains (which are arbitrarily described in the reference axis)
  % align with [001]-direction (z-axis) of the reference frame; 

  % Inversely, its inverse (or transpose) maps an arbitrary direction of crystal
  % plane to [001]-drirection, (and note that this corresponds to the view from
  % transverse axis) 

  % p.s. if x=[0,0,1], then you should see a single component, as all [110] crystall is aligned in z 
  % The above point is checked 
  % Rather, you are viewing from the x-axis 
  x=[1,0,0]; x=transpose(x); 

  % 'vec' is a directional vector, i.e. size is 1 
  vec= transpose(RotGrain)*x; 

  direction(k,1)=vec(1); 
  direction(k,2)=vec(2); 
  direction(k,3)=vec(3); 

  u = direction(k,3) - direction(k,2); 
  v = direction(k,2) - direction(k,1);
  w = direction(k,1);

  rgb(k,:)=[u,v,w]; 
  rgb(k,:) = rgb(k,:)/ max(abs(rgb(k,:))); % normalize the colors  
  rgb(k,:) = abs(rgb(k,:)); 

  u1(ind2)=rgb(k,1);
  u2(ind2)=rgb(k,2);
  u2(ind2)=rgb(k,3);

end

data(:,:,1)=u1; 
data(:,:,2)=u2; 
data(:,:,3)=u3; 

fig=figure; 

image(data);
axis square;
axis off;