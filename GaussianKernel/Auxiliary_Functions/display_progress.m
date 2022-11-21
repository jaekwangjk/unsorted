function display_progress(k,N,x,d)

if nargin < 4
    d = 0;
end

if and(mod(floor(k/N*100*10^d),x)==0,floor(k/N*100*10^d)~=floor((k-1)/N*100*10^d))
    dispstat(sprintf(['Progress %1.',num2str(d),'f%%'],floor(k/N*100*10^d)/10^d),'timestamp');
end

if floor(k/N*100*10^d) == 100*10^d
    dispstat('Finished.','keepprev');
end
end