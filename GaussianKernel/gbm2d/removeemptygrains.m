function [newgrains,newid] = removeemptygrains(grains,id)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [newgrains,newid] = removeemptygrains(grains,id)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removes empty grains from the structure "grains".

N = size(grains,1); % Number of grains.

count = 0; % Number of non-empty grains, initialized.
for k=1:N % Loop over the grains.

  maxlev = max(grains{k,2}); % Max. val. of level set function.
  if maxlev > 0
    count = count + 1;
    newgrains{count,1} = grains{k,1};
    newgrains{count,2} = grains{k,2};
    newgrains{count,3} = grains{k,3};
    newgrains{count,4} = grains{k,4};
    newid(count,1) = id(k,1);
    grains{k,1} = [];
    grains{k,2} = [];
    grains{k,3} = [];
    grains{k,4} = [];
  end % end if.

end % k. Loop over grains ends.