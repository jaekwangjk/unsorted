function [newgrains,newid] = removeemptygrainsoriginal(grains,id)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [newgrains,newid] = removeemptygrainsoriginal(grains,id)
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
    newid(count,1) = id(k,1);
    grains{k,1} = [];
    grains{k,2} = [];
    grains{k,3} = [];
    grains{k,4} = [];
    % this allows the function 'gbm2doriginal_sim' to be initialized with
    % the same grain as in ''gbm2_sim' which is prepared to store
    % convolution values from two distinct kernels
  end % end if.

end % k. Loop over grains ends.