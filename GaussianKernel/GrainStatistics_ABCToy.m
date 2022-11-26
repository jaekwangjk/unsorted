
clear
close all

%% Load grain 
grains = load("grainFiles/grain_later_Test2_240.mat").grains; 

%%  Get the mean area for all type 
Ng = size(grains,1); %number of grain 
dims = [2048, 2048]; 

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


% Extract non zero element 
Areas = nonzeros(Areas); 
AT= Areas; 
AreaMean  = mean(Areas); 
Areas=Areas/AreaMean;

% Extract non zero element 
Areas_Atype = nonzeros(Areas_Atype); 
AA= Areas_Atype; % Temp which keeps the actual area not reduced area
Areas_Atype = Areas_Atype/AreaMean;

Areas_Btype = nonzeros(Areas_Btype); 
AB= Areas_Btype; 
Areas_Btype = Areas_Btype/AreaMean;

Areas_Ctype = nonzeros(Areas_Ctype); 
AC= Areas_Ctype; 
Areas_Ctype = Areas_Ctype/AreaMean;

Ng = size(Areas,1);
NgA = size(Areas_Atype); NgA=NgA(1);
NgB = size(Areas_Btype); NgB=NgB(1);
NgC = size(Areas_Ctype); NgC=NgC(1);


%% Draw histogram...


Nbins = 31; 


% construct x
x=zeros(Nbins,1);
areaBin = linspace(0, 5, Nbins+1); 
for k=1:Nbins
    x(k) = (areaBin(k)+areaBin(k+1))* 0.5; 
end

y=histcounts(Areas,areaBin); 
y=y';
y=y/Ng; 


yA= histcounts(Areas_Atype,areaBin); 
yA=yA';
yA=yA/NgA; 


yB= histcounts(Areas_Btype,areaBin); 
yB=yB';
yB=yB/NgB; 


yC= histcounts(Areas_Ctype,areaBin); 
yC=yC';
yC=yC/NgC; 



h=figure; 
plot(x,y,'ko-','linewidth',2,'DisplayName','Overall'); hold on; 
plot(x,yA,'ro-','linewidth',2,'DisplayName','A type'); hold on; 
plot(x,yB,'go-','linewidth',2,'DisplayName','B type'); hold on; 
plot(x,yC,'bo-','linewidth',2,'DisplayName','C type'); hold on; 



%% Setting Plot 
legend()
xlabel('A/A_{avg}','fontsize',24)
ylabel('Count','fontsize',24)


% axis([0,5,0,0.35]); 

set(gca,'fontsize',40,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,800,600])


set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'Poly3','-dpdf','-r0')



%%  Count plot
% 
% h=figure; 
% 
% X = categorical({'A type','B type','C type'});
% X = reordercats(X,{'A type','B type','C type'}); % To preserve the category order in bar function 
% y=[length(Areas_Atype)/Ng, length(Areas_Btype)/Ng, length(Areas_Ctype)/Ng];
% bar(X,y) 
% 
% %xlabel('A/A_{avg}','fontsize',24)
% ylabel('Fraction Count','fontsize',24)
% set(gca,'fontsize',40,'ticklabelinterpreter','latex')
% set(gcf,'units','points','position',[30,30,800,600])

%% Count fraction and area fraction plot 

areaFraction = [sum(AA), sum(AB), sum(AC)]; 
countFraction = [NgA/Ng, NgB/Ng, NgC/Ng]; 

xref = linspace(0,1,50); 
yref = linspace(0,1,50); 

h=figure; 
plot(areaFraction, countFraction,'ro','Markersize',10,'linewidth',2); hold on; 
plot(xref,yref,'k--'); 
xlabel('Area Fraction','fontsize',24)
ylabel('Count Fraction','fontsize',24)
set(gca,'fontsize',20,'ticklabelinterpreter','latex')
set(gcf,'units','points','position',[30,30,800,600])



