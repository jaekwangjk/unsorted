%% Load grain 
clear 
close all 

grains = load("grainFiles/Poly_orient_c1_init.mat").grains; 
dims=[2048 2048]; % dimensions of grid 

Ng = size(grains,1); %number of grain 
Area_and_Orient = zeros(Ng,2); 

%% 코딩!! info Oreintation를 저장한 뒤 Area와 함께 non zero element를 추출해야함..!!

posind = cell(1,Ng); 
for k=1:Ng 
    
    Area_and_Orient(k,2) = mod(k, 90); 
    
    ind = grains{k,1}; % Pixels within a nhd. of grain.
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    posind{k} = ind(val>0); % Pixels in the interior of grain.
    Area_and_Orient(k,1) = length(posind{k})/prod(dims);
    
end

% How orientation is defined from ID  
aliveID = find(Area_and_Orient(:,1)); % Grains with area > 0 


Area_and_Orient = Area_and_Orient(aliveID,:); 


meanArea = mean(Area_and_Orient(:,1)); 

for k =1:length(aliveID) 
    Area_and_Orient(k,1)= Area_and_Orient(k,1)/meanArea;     
end


% For each orientaiton 
% sum area add


Area_by_Orient =zeros(90,1); 
Orient = zeros(90,1); 
Count_by_Orient = zeros(90,1); 

for k=1:90
    Orient(k) = k-1; 
end

for k =1:length(aliveID) 
    
    psOrient = Area_and_Orient(k,2); 
    Area_by_Orient(psOrient+1)= Area_by_Orient(psOrient+1) + Area_and_Orient(k,1); 
    Count_by_Orient(psOrient+1)= Count_by_Orient(psOrient+1) +1; 
end


%for k=1:90
%    Area_by_Orient(k) = Area_by_Orient(k) / Count_by_Orient(k); 
%end


% N_alive = length(aliveID); 
% Orientations = zeros(N_alive,1); 
% 
% for k=1:N_alive
% 
%     Orientations(k) = mod(aliveID(k), 90); 
%     
% end
% 
% 
% h=histogram(Orientations,90);


%% Setting Plot 
legend()
plot(Orient, Area_by_Orient,'ro-','linewidth',2); 
xlabel('Orientation','fontsize',15)
ylabel('Reduced Area (Sum)','fontsize',15)

%axis([-1,91,0,2]); 

set(gca,'fontsize',25,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,800,600])


%set(h,'Units','Inches');
%pos = get(h,'Position');
%set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'Poly3','-dpdf','-r0')
