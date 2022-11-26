clear
close all 

addpath(genpath(pwd)) % add all subfolders to the path

%% Hyper parameters 
N=5000; % may be number of grain 
dims=[2048 2048]; % dimensions of grid 

%N=200;
%dims=[500 500]; 
Rgrains = 10; % number of pixel to outgrow 
Rfamilies = 30; % minimum pixel distance to have two grains in same family 

th_ang_for_merge = 0.001; % th_ang = Threshold angle to merge grains.
angBrandon= 30; % not necessarily needed.....the value is not used all to be honest... 
%% Prepare workspaces used globally 

global KERNEL Z WORKSPACE;

I = sqrt(-1); % Imaginary number i.
w=exp(I*2*pi/dims(1)); % nth root of unity.
[x,y] = meshgrid(1:dims(1),1:dims(1)); x = x'; y = y';
dt = 10/(dims(1)^2); % suggested time step size 

KERNEL = exp( -dt*dims(1)*dims(1)*( 4-w.^(x-1)-w.^(1-x)-w.^(y-1)-w.^(1-y) ) );

Z = -ones(dims(1),dims(2)); % workspace needed in loc_levset_to_volfluid.c.
WORKSPACE = int32(-ones(dims)); % workspace needed for pgrow.c.

global work_x work_y work_bx work_by work_candidate_x work_candidate_y; % Needed for pgrow3.c.
work_x = int32(zeros(prod(dims),1));
work_y = int32(zeros(prod(dims),1));
work_bx = int32(zeros(prod(dims),1));
work_by = int32(zeros(prod(dims),1));
work_candidate_x = int32(zeros(prod(dims),1));
work_candidate_y = int32(zeros(prod(dims),1));

%% Initialization

rn=2; % random number assigned for voronoi tessellation 

[grains,ori] = Module_INIT_initialvoronoidata2d(N,dims,Rgrains,th_ang_for_merge,rn);
ID = (1:1:size(grains,1))';
pauseTime = 0.5;
% Module_DataProcess_showgrains_withLabel(grains,dims,ID,pauseTime); 


%% save initial grain and load grain into a 'mat'file
% filename = './TestGrain';
% save(filename,'grains','dims','-v7.3')


%% Time marching

fileName1 = './png/polycrystal_';
for t = 1:301% Main time iteration starts.
%for t = 1:2% Main time iteration starts.
    
    status = append(int2str(t),' Time step is running');
    disp(status); 
    
    %Module_DataProcess_showgrains_withLabel(grains,dims,ID,pauseTime); 
    %Module_DataProcess_savegrainsfigure_withLabelRGB(grains,dims,ID,fileName1,t); 
    
    display=0;
    [families,famgrains] = Module_TECH_grains2families2doriginal(grains,Rfamilies,dims,display); 
    Nf = size(families,1); % Number of families.
    
    % %% save grain and load grain into a 'mat'file
    if(t>299  && mod(t,30)==0)
    %if(t==1)
    filename = './grain_later_Poly_confi2_';
    filename = append(filename,int2str(t)); 
    save(filename,'grains','dims','-v7.3')
    end 
    
    % Convolving families     
    for k=1:Nf % Loop over families.
        
        cval = Module_TECH_convolvefamily2d(families{k,1},families{k,2},dims);
        % Distribute convolution values to the grains contained in this family:
        numberofgrains = size(famgrains{k},1); % Number of grains in this family.
        listofgrains = famgrains{k};           % Column vector of grain indices.
        
        % Family convolution result로부터 grain convolution을 얻는 과정..
        for ell = 1:numberofgrains % Loop over grains contained in this family.
            label = listofgrains(ell);
            ind = grains{label,1};
            grains{label,3} = cval(ind); % Read off and record the convolution vals.
        end % (for ell) Loop over grains ends.
    end % (for k) Loop over families ends.
    
    
    % REDISTRIBUTION STEP:
    presence = CModule_get_nhd_grains2doriginal(grains,dims(1)*dims(2));
    
    % Change the function here!
    %updatelevelsetdata2dToyTest2(presence,grains,ID,ori,angBrandon);
    %updatelevelsetdata2dIsotropic(presence,grains,ID,ori,angBrandon);
    updatelevelsetdata2dPolyCrystal(presence,grains,ID,ori,angBrandon);
    
    % do not use this. This is a nasty fucntion 
    %[grains,ID] = removeemptygrainsoriginal(grains,ID);
    
    
    %% REFRESH GRAIN BUFFERS:
    % Write this as a function 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level set data for each grain extends beyond the interface,
    % to a tubular neighborhood. That neighborhood has to be
    % updated once the interface has (potentially) moved.
    N = size(grains,1); % Number of grains.
    for k=1:N % Loop over grains.
        ind = grains{k,1}; % Pixels within a nhd. of grain.
        val = grains{k,2}; % Lev. set. vals. at those pixels.
        cval = grains{k,3}; % Convolution vals. at those pixels.
        Z(ind) = val;      % Lev. set. representation on grid.
        posind = ind(val>0); % Pixels in the interior of grain.
        [x,y] = ind2sub(dims,posind);
        [x2,y2] = pgrow3(int32(x),int32(y),Rgrains,WORKSPACE,...
            work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y); % Dilation.
        ind2 = sub2ind(dims,x2,y2);
        val2 = Z(ind2); % Level set vals.
        Z(ind2) = -1; % Reset Z back to all -1's. (outside)
        Z(ind) = cval - 1; % Convolution values - 1.(if it was belong to buffered grain)
        cval2 = Z(ind2); % Convolution vals - 1.
        Z(ind2) = -1;
        grains{k,1} = ind2;   % Refresh grain's data structure.
        grains{k,2} = val2;   % Ditto.
        grains{k,3} = cval2 + 1; % Ditto. 0 is substutived
    end % (for k). Loop over grains ends.
    
    

end

