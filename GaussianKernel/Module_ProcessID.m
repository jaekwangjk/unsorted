function newID = Module_ProcessID(ID)


Ng = length(ID); 
newID = zeros(Ng, 1);


for k=1:Ng
    
    
    grainType = mod(k,3); 
    
    if(grainType==0) % Domiant mode 
        newID(k,1) = fix( rand(1) * 10);
    else
        newID(k,1) = fix( 10 + rand(1) *  80); 
    end
    
    
    
    
end


end