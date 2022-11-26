function [x,y]=getAreaStatistics(grains,dims)

Nbins = 21; 

% construct x
x=zeros(Nbins,1);
areaBin = linspace(0, 5, Nbins+1); % These are edges of bin 
for k=1:Nbins
    x(k) = (areaBin(k)+areaBin(k+1))* 0.5; 
end


Ng = size(grains,1); %number of grain 

posind = cell(1,Ng);
area_grain = zeros(Ng,1);

for k=1:Ng 
    ind = grains{k,1}; % Pixels within a nhd. of grain.
    val = grains{k,2}; % Lev. set. vals. at those pixels.
    posind{k} = ind(val>0); % Pixels in the interior of grain.
    area_grain(k) = length(posind{k})/prod(dims);
end

% Extract non zero element 
area_grain = nonzeros(area_grain); 
AreaMean  = mean(area_grain); 
area_grain=area_grain/AreaMean;

Ng = size(area_grain,1);

y= histcounts(area_grain,areaBin); 
y=y';
y=y/Ng; 



end