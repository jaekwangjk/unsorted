function save_variable(filename,var,name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  save_variable(filename,var,name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save_variable saves the variable var in filename
% under the name 'name'


aux.(name)= var;
save(filename,'-struct','aux',name,'-append')

end
