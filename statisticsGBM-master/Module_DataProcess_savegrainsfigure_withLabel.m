function Module_DataProcess_savegrainsfigure_withLabel(grains,dims,id,fileNumber)


Ng = size(grains,1); % Number of grains.


x=linspace(1,dims(1),dims(1)); 
y=linspace(1,dims(1),dims(1)); 
[X,Y]=meshgrid(x,y);  

data=zeros(dims(1)*dims(1),1); 


for k=1:Ng % Loop over grains.
  ind = grains{k,1}; % Pixels in a nhd. of the grain.
  val = grains{k,2}; % Level set values. 
  ind2 = ind(val>0); % Pixels in the interior of grain.
  data(ind2) = id(k); 
end


data=reshape(data,dims(1),dims(1)); 


fig=figure; 

h=surf(X,Y,data); 
set(h,'edgecolor','none')


axis([1 dims(1) 1 dims(1)])
pbaspect([1 1 1])
view(2)
ax = gca;
set(gca,'XColor', 'none','YColor','none') %ticks off


outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
formatSpec='%03d'; 

FnameString = append('grainFigureIso',num2str(fileNumber, formatSpec)) ; 
print(fig,FnameString,'-dpng','-r0')

close(fig); 


end