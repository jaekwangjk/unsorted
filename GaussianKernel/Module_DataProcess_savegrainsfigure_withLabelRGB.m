function Module_DataProcess_savegrainsfigure_withLabelRGB(grains,dims,id,fileName, fileNumber)


Ng = size(grains,1); % Number of grains.

% matlab rgb colors takes values in between [0,1]
% the following code will well resolve large polycrystals....

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

fig=figure; 


image(data);
axis square;
axis off;

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
formatSpec='%03d'; 

FnameString = append(fileName,num2str(fileNumber, formatSpec)) ; 
print(fig,FnameString,'-dpng','-r0')

close(fig); 


end
