
%% To obtain side distribution...
close all 
%% Load grain 
grains = load("grainFiles/grain_later_Poly_300.mat").grains; 
dims=[2048, 2048]; % dimensions of grid 

% Neighbor === ADJCENCY MATRIX...
% construct domain pixel..
Ng = size(grains,1); %number of grain 
data=zeros(dims(1)*dims(1),1); 

%% Determine the area of each type of grain 


Areas = zeros(Ng,1); 


Areas_Atype = zeros(Ng,1); 
Areas_Btype = zeros(Ng,1); 
Areas_Ctype = zeros(Ng,1); 


% Extract areas of grain...
posind = cell(1,Ng); 
for k=1:Ng 
    ind = grains{k,1}; % Pixels within a nhd. of grain.
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    posind{k} = ind(val>0); % Pixels in the interior of grain.
    Areas(k) = length(posind{k})/prod(dims);
   
    type = mod(k,3); 
    
    if(type == 0)
      Areas_Atype(k) = Areas(k); 
    elseif(type ==1)
      Areas_Btype(k) = Areas(k);
    else
      Areas_Ctype(k) = Areas(k); 
    end    
        
end


Areas = nonzeros(Areas); 
Areas_Atype = nonzeros(Areas_Atype); 
Areas_Btype = nonzeros(Areas_Btype); 
Areas_Ctype = nonzeros(Areas_Ctype); 




%% Adjacency Matrix 
ADJ = zeros(Ng,Ng); 

for k=1:Ng % Loop over grains.
  ind = grains{k,1}; % Pixels in a nhd. of the grain.
  val = grains{k,2}; % Level set values. 
  ind2 = ind(val>0); % Pixels in the interior of grain.
  data(ind2) = k; 
end

data=reshape(data,dims(1),dims(1)); 

xgrid=[1, 0, -1, 0];
ygrid=[0, 1, 0, -1];


for k=1:dims(1)
   for j=1:dims(1)
       
      cx = k; cy=j; 
      currentIndex = data(cx,cy);
      
      for m = 1:4
       xn = cx+xgrid(m);  xn = mod(xn+dims(1), dims(1))+1;  
       yn = cy+ygrid(m);  yn = mod(yn+dims(1), dims(1))+1; 
       neighborIndex = data(xn,yn);
      
       if(neighborIndex ~= currentIndex)
          ADJ(neighborIndex, currentIndex)= ADJ(neighborIndex,currentIndex)+1; 
          ADJ(currentIndex, neighborIndex)= ADJ(currentIndex,neighborIndex)+1; 
       end
    
      end % end loop over neighbors...
      
   end
end


for k=1:Ng
   for j = 1:Ng
      if(ADJ(k,j)> 5)
          ADJ(k,j)=1; 
      else
          ADJ(k,j)=0;
      end
   end
end




%% Start "Extract Side" process over here....

%% ID and is mixed...

side =zeros(1,Ng);
for k=1:Ng
   side(k) = sum(ADJ(k,:));  
end

Side_Atype = zeros(Ng,1); 
Side_Btype = zeros(Ng,1); 
Side_Ctype = zeros(Ng,1); 


for k=1:Ng 
     
     type = mod(k,3); 
 
     if(type == 0)
       Side_Atype(k) = side(k);
     elseif(type ==1)
       Side_Btype(k) = side(k); 
     else
       Side_Ctype(k) = side(k); 
     end    
         
end

 
Side_Atype = nonzeros(Side_Atype); NgA = size(Side_Atype); NgA=NgA(1);
Side_Btype = nonzeros(Side_Btype); NgB = size(Side_Btype); NgB=NgB(1);
Side_Ctype = nonzeros(Side_Ctype); NgC = size(Side_Ctype); NgC=NgC(1);


sideBin = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5]; 

x=zeros(length(sideBin)-1,1);

for k=1:length(sideBin) -1
  x(k) = (sideBin(k)+sideBin(k+1))* 0.5; 
end

side=nonzeros(side); 

Ng = size(side); Ng=Ng(1);


y=histcounts(side,sideBin); 
y=y';
y=y/Ng; 
disp("y"); disp(y)
plot(x,y,'ko-','linewidth',2,'DisplayName','Overall'); hold on; 

yA=histcounts(Side_Atype,sideBin); 
yA=yA';
yA=yA/NgA; 
plot(x,yA,'ro-','linewidth',2,'DisplayName','Atype'); hold on; 

yB=histcounts(Side_Btype,sideBin); 
yB=yB';
yB=yB/NgB; 
plot(x,yB,'go-','linewidth',2,'DisplayName','Btype'); hold on; 

yC=histcounts(Side_Ctype,sideBin); 
yC=yC';
yC=yC/NgC; 
plot(x,yC,'bo-','linewidth',2,'DisplayName','Ctype'); hold on; 


%y_weighted = ((yA*NgA) * sum(Areas_Atype) + (yB*NgB) * sum(Areas_Btype) + (yC*NgC) * sum(Areas_Ctype))/Ng ; 

y_weighted = yA * sum(Areas_Atype) + yB * sum(Areas_Btype) + yC* sum(Areas_Ctype) ; 

plot(x,y_weighted,'ko--','linewidth',2,'DisplayName','Area Weighted mean'); hold on; 




%% Setting Plot 
legend()
xlabel('Side','fontsize',15)
ylabel('Probability density','fontsize',15)

%axis([0,5,0,0.35]); 

set(gca,'fontsize',25,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,800,600])



%% Check the area fraction...

