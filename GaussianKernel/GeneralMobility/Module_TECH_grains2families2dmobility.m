function [families,famgrains] = Module_TECH_grains2families2dmobility(grains,R,dims)

% Note that this code can easily include the single kernel code. 

global WORKSPACE;
global work_x work_y work_bx work_by work_candidate_x work_candidate_y;

N = size(grains,1); % Number of grains.
maxFamilies = 1000;
famgrains = cell(maxFamilies,1);

grainlabels = zeros(dims); % Grid describing which grain each pixel belongs to.

for k=1:N % Loop over grians.
  ind = grains{k,1};  % List of pixels in a nhd of grain k.
  val = grains{k,2}; % Level set function val. at pixels in grain k.
  posind = ind(val>0); % Pixels contained in the interior of grain k.
  grainlabels(posind)=k; % Populate grid of grain labels.
end

remaininggrains = zeros(1,N); % Grains not allocated to families so far.
maxlabel = 0;
families = cell(maxFamilies,4); % Maximum of maxFamilies families anticipated.

collect_ind = zeros(dims(1)*dims(2),1);
collect_val = zeros(dims(1)*dims(2),1);
collect_cval1 = zeros(dims(1)*dims(2),1);
collect_cval2 = zeros(dims(1)*dims(2),1);

label = 0;
while N>sum(remaininggrains>0) % Loop over families.
    label = label+1;
    intersects = zeros(1,N);
    collected = 0;
    
    listofgrains = [];  % Grains contained in this family, initialized.
    numberofgrains = 0; % Number of grains in this family, initialized.
    
    for k=1:N % Loop over the grains.
        
        if (intersects(k)==0 && remaininggrains(k)==0) % If grain k is unallocated
            % and is far from other
            % grains in this family
            
            % Remove grain k from remaining grains:
            remaininggrains(k) = 1;
            
            % Grow grain k:
            ind = grains{k,1};
            [x,y] = ind2sub(dims,ind);
            [x2,y2] = pgrow3(int32(x),int32(y),R,WORKSPACE,...
                work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y); % Dilation.
            ind2 = sub2ind(dims,x2,y2); % Pixels in R neighborhood of current grain.
            
            % Move grain to family label:
            pixingrain = size(ind,1);
            collect_ind(collected+1:collected+pixingrain) = ind;
            collect_val(collected+1:collected+pixingrain) = grains{k,2};
            collect_cval1(collected+1:collected+pixingrain) = grains{k,3};
            
            %added 
            collect_cval2(collected+1:collected+pixingrain) = grains{k,4};
            collected = collected + pixingrain;
            
            maxlabel = label;
            
            % Grains ineligible for inclusion in current family due to proximity:
            gind = grainlabels(ind2); % Labels of grains found in R neighborhood.
            if max(gind) > 0
                ind2 = gind(gind>0); % Eliminate zeros from labels.
                intersects(ind2) = 1; % Those labels marked as ineligible.
            end
            
            % Add grain k to the list of grains contained in this family:
            numberofgrains = numberofgrains + 1;
            listofgrains(numberofgrains,1)= k;
            
        end % end if.
    end % (for k) Loop over grains ends.
    
    families{label,1} = collect_ind(1:collected);
    families{label,2} = collect_val(1:collected);
    families{label,3} = collect_cval1(1:collected);
    
    %added 
    families{label,4} = collect_cval2(1:collected);
    famgrains{label} = listofgrains;
end % (for label). Loop over families ends.

families = families(1:maxlabel,1:4);
famgrains = famgrains(1:maxlabel);