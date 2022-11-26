function TestVisual(grains,dims,id,pause_v)

Ng = size(grains,1); % Number of grains.

%x=linspace(1,dims(1),dims(1)); 
%y=linspace(1,dims(1),dims(1)); 
%[X,Y]=meshgrid(x,y);  

%data=zeros(dims(1)*dims(1),1); 

resolution = 1.0/ (fix(Ng/3)*2.5); 

u1 = zeros(dims);
u2 = zeros(dims);
u3 = zeros(dims); 

for k=1:Ng % Loop over grains.
  ind = grains{k,1}; % Pixels in a nhd. of the grain.
  val = grains{k,2}; % Level set values. 
  ind2 = ind(val>0); % Pixels in the interior of grain.
  
  
  if(mod(id(k),3)==0)
      u1(ind2) = 0.5 + fix(id(k)/3) * resolution; 
  elseif( mod(id(k),3)==1)
      u2(ind2) = 0.5 + fix(id(k)/3) * resolution; 
  else %(mod(id(k),3)==2)
      u3(ind2) = 0.5 + fix(id(k)/3) * resolution; 
  end
 
end


data(:,:,1)=u1; 
data(:,:,2)=u2; 
data(:,:,3)=u3; 

image(data);
axis square;
axis off;
pause(pause_v);

end


% %%%%%
% for k=1:N % Loop over grains.
%   ind = grains{k,1}; % Pixels in a nhd. of the grain.
%   val = grains{k,2}; % Level set values. 
%   ind2 = ind(val>0); % Pixels in the interior of grain.
%   ang = ori(id(k)); % read orientatio values 
%   u1(ind2) = (1+cos(ang))/2;
%   u2(ind2) = (1+sin(ang))/2;
% end
% 
% % Colormap does not work because U is provided with RGB 
% % You need to figure out nice combination of RGB. 
% U(:,:,1) = u1; % R
% U(:,:,2) = u2; % G
% U(:,:,3) = ones(dims); %B
% clf % Clear the figure 
% image(U);
% axis square;
% colormap gray
% colorbar 