% This code creates nucleate grain from clustered grain structure 
clear all 

addpath(genpath(pwd)) % add all subfolders to the path


N=1000;
dims=[2000 2000]; % dimensions of grid 
Rgrains = 10; % number of pixel to outgrow 
rn=1; 
FactorNc=20; 

[grains,ori,cluster] = Module_INIT_initialdoubleVoronoidata2d(N,dims,Rgrains,rn, FactorNc);

ID = (1:1:size(grains,1))';

%% reset orientaiton % It should run with this...
for i=1:length(ori)    
 if(mod(cluster(i),6)==0 )
 ori(i) = 3 * rand; 
 elseif(mod(cluster(i),6)==1)
 ori(i) = 12 +  3 * rand; 
 elseif(mod(cluster(i),6)==2)
 ori(i) = 22 +  3 * rand; 
 elseif(mod(cluster(i),6)==3)
 ori(i) = 32 +  3 * rand; 
 elseif(mod(cluster(i),6)==4)
 ori(i) = 42 +  3 * rand; 
 elseif(mod(cluster(i),6)==5)
 ori(i) = 62 +  3 * rand; % from 62 to 65
 end  
end

ori = ori * pi/180 ; % make radian 



% fileName1 = 'checkup_';
% Module_DataProcess_savegrainsfigure_withLabelRGB_ori(grains,dims,ID,ori,fileName1,0);

Ng = size(grains,1); %number of grain 

fileName1 = 'original_';
Module_DataProcess_savegrainsfigure_withLabelRGB_ori(grains,dims,ID,ori,fileName1,0);





%% Adjacency Matrix 

ADJ = zeros(Ng,Ng); 
data=zeros(dims(1)*dims(1),1);

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

%%
newOri = ori; % first just copy the new ori 
nucleated = zeros(Ng,1); % first just copy the new ori 
rng(1); 
for i=1:Ng
    neighbors=find(ADJ(i,:)>0); 
    n= length(neighbors); 
    orientA=ori(i); 
    for j=1:n
        nIndex = neighbors(j); 
        if(nIndex>i) % Just to one random flip, for each relation 
            orientB= ori(nIndex);
            misorient = abs(orientA - orientB);
            if(misorient>10 * pi/180)
             coinToss = rand; 
             if(coinToss>0.8)
                type=randi([0,5],1); 
                nucleated(i)=1; 
                if(type==0 )
                    newOri(i) = (3 * rand) *pi/180; % from 0 to 5
                elseif(type==1)
                    newOri(i) = (12 +  3 * rand)*pi/180; %from 15 to 20
                elseif(type==2)
                    newOri(i) = (22 +  3 * rand)*pi/180; % from 25 to 30
                elseif(type==3)
                    newOri(i) = (32 +  3 * rand)*pi/180; % from 35 to 40
                elseif(type==4)
                    newOri(i) = (42 +  3 * rand)*pi/180; % from 45 to 50
                elseif(type==5)
                    newOri(i) = (62 +  3 * rand)*pi/180; % from 65 to 70
                end 

             end
            end
        end
    end
end

% fileName1 = 'nucleate_';
% Module_DataProcess_savegrainsfigure_withLabelRGB_ori(grains,dims,ID,newOri,fileName1,0);

% % filename='./InitialCondition/init_grain_nucleate_1';  
% % save(filename,'grains','dims','-v7.3')
% % 
% % 
% % filename='./InitialCondition/ori_case_d_1';  
% % save(filename,'ori','dims','-v7.3')
% % 
% % 
% % filename='./InitialCondition/nucleated_case_d_1';  
% % save(filename,'nucleated','dims','-v7.3')

