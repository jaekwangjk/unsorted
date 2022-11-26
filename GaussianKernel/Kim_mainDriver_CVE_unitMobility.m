clear
close all 

addpath(genpath(pwd)) % add all subfolders to the path

% Surface energy tension must be considered in the two places. 
% When it generate a kernel 
% when it updatelevelset2D

% Is this a issue with a Kerel Width
% Module_TECH_grains2families2dmobility(grains,Rfamilies,dims,display); 
%[families,famgrains] = Module_TECH_grains2families2dmobility(grains,Rfamilies,dims,display); 
% cval = Module_TECH_convolvefamily2dmobility(families{k,1},families{k,2},dims);
% presence = CModule_get_nhd_grains2dmobility(grains,dims(1)*dims(2));
% CModule_updatelevelsetdata2d_RSenergy(presence,grains,ID,alpha, beta, ori,angBrandon,option);
    

%% Hyper parameters 
%N=5000; % may be number of grain 
%dims=[2048 2048]; % dimensions of grid 

N=1000;
dims=[2000 2000]; 
n = dims(1); % Size of computational grid.

Rgrains = 10; % number of pixel to outgrow 
Rfamilies = 30; % minimum pixel distance to have two grains in same family 


% global KERNEL Z WORKSPACE;
global KERNELalpha KERNELbeta Z WORKSPACE;


I = sqrt(-1); % Imaginary number i.
w=exp(I*2*pi/dims(1)); % nth root of unity.x
[x,y] = meshgrid(1:dims(1),1:dims(1)); x = x'; y = y';
dt = 10/(dims(1)^2); % suggested time step size 

%% Kernel Option 
% understand how alpha, beta is being chosen 
 
CVE = importdata("FCC_110STGB.txt");

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
rn=1; % random number assigned for voronoi tessellation 
[grains,ori] = Module_INIT_initialdoubleVoronoidata2d(N,dims,Rgrains,rn,40);
ID = (1:1:size(grains,1))';

% reset orientaiton % It should run with this...
for i=1:length(ori)
    
 if(mod(i,7)==0 )
 ori(i) = 1 * randi([0,5]); % from 0 to 5
 %ori(i) = 5 * rand(1);   %randi([0,5]); % from 0 to 5
 elseif(mod(i,7)==1)
 %   ori(i) = 15+ 5 * rand(1); 
 ori(i) = 15 +  1 * randi([0,5]); %from 15 to 20
 elseif(mod(i,7)==2)
 %  ori(i) = 25 * rand(1); 
 ori(i) = 25 +  1 * randi([0,5]); % from 25 to 30
 elseif(mod(i,7)==3)
 %   ori(i) = 35+ 5 * rand(1); 
 ori(i) = 35 +  1 * randi([0,5]); % from 35 to 40
 elseif(mod(i,7)==4)
 %   ori(i) = 45+ 5 * rand(1); 
 ori(i) = 45 +  1 * randi([0,5]); % from 45 to 50
 elseif(mod(i,7)==5)
 %  ori(i)=  55+ 5 * rand(1); 
 ori(i) = 55 +  1 * randi([0,5]); % from 55 to 60
 elseif(mod(i,7)==6)
 %   ori(i) = 60+ 5 * rand(1);  
 ori(i) = 62 +  1 * randi([0,3]); % from 65 to 70
 end
      
end

ori = ori * pi/180 ; % make radian 
distinctOri= unique(ori); 

% Build Kernel for unit mobility 
option =1; 
%option =-1;  % Read Sochekley. 
[alpha, beta] = Kim_kernel_widths2d(distinctOri,CVE,option);


KERNELalpha = exp(-dt*alpha*n*n*(4-w.^(x-1)-w.^(1-x)-w.^(y-1)-w.^(1-y)));
KERNELbeta = exp(-dt*beta*n*n*(4-w.^(x-1)-w.^(1-x)-w.^(y-1)-w.^(1-y)));

disp('alpha'); disp(alpha);
disp('beta'); disp(beta); 


% reset orientaiton % It should run with this...
for i=1:length(ori)
    
 if(mod(i,7)==0 )
   ori(i) = 5 * rand(1);   %randi([0,5]); % from 0 to 5
 elseif(mod(i,7)==1)
   ori(i) = 15 + 5 * rand(1); 
 elseif(mod(i,7)==2)
   ori(i) = 25 + 5 * rand(1); 
 elseif(mod(i,7)==3)
   ori(i) = 35 + 5 * rand(1); 
 elseif(mod(i,7)==4)
   ori(i) = 45 + 5 * rand(1); 
 elseif(mod(i,7)==5)
   ori(i)=  55 + 5 * rand(1); 
  elseif(mod(i,7)==6)
   ori(i) = 62 +  3 * rand(1); % from 62 to 65
 end
 
end

ori = ori * pi/180 ;


%% save initial grain and load grain into a 'mat'file
% filename = './TestGrain';
% save(filename,'grains','dims','-v7.3')
% 
% filename = './oriman';
% save(filename,'ori','-v7.3')


%% Time marching

fileName1 = './png/cvemobility_';
for t = 1:5% Main time iteration starts.
%for t = 1:1% Main time iteration starts.
    
    status = append(int2str(t),' Time step is running');
    disp(status); 
    
   
    %Module_DataProcess_showgrains_withLabel(grains,dims,ID,pauseTime); 
    Module_DataProcess_savegrainsfigure_withLabelRGB_ori(grains,dims,ID,ori,fileName1,t); 
    
    display=0;
    
    
    % These two are different, as the second prepares cell structure for
    % saving the second kernel...
	[families,famgrains] = Module_TECH_grains2families2dmobility(grains,Rfamilies,dims); 
    Nf = size(families,1); % Number of families.
   
     
    % Convolving families     
    for k=1:Nf % Loop over families.
        
        cval = Module_TECH_convolvefamily2dmobility(families{k,1},families{k,2},dims);
        % Distribute convolution values to the grains contained in this family:
        numberofgrains = size(famgrains{k},1); % Number of grains in this family.
        listofgrains = famgrains{k};           % Column vector of grain indices.
        
        % Family convolution result로부터 grain convolution을 얻는 과정..
        for ell = 1:numberofgrains % Loop over grains contained in this family.
            label = listofgrains(ell);
            ind = grains{label,1};
            cvalaux = cval{1}; 
            grains{label,3} = cvalaux(ind); % Read off and record the convolution vals.
            cvalaux = cval{2};
            grains{label,4} = cvalaux(ind); 
        end % (for ell) Loop over grains ends.
    end % (for k) Loop over families ends.
    
    
    % REDISTRIBUTION STEP:
    % presence = CModule_get_nhd_grains2doriginal(grains,dims(1)*dims(2));
    presence = CModule_get_nhd_grains2dmobility(grains,dims(1)*dims(2));

    option=1; angBrandon =30; 
    
    ori_in_deg = ori * 180/pi; % This has to be handled correctly
    CModule_updatelevelsetdata2d_CVE(presence,grains,ID,ori_in_deg,alpha,beta,angBrandon,option,CVE);
    
    
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
        cval1 = grains{k,3}; % Convolution vals. at those pixels.
        cval3 = grains{k,4}; % Convolution vals. at those pixels.
        
        Z(ind) = val;      % Lev. set. representation on grid.
        posind = ind(val>0); % Pixels in the interior of grain.
        [x,y] = ind2sub(dims,posind);
        [x2,y2] = pgrow3(int32(x),int32(y),Rgrains,WORKSPACE,...
            work_x,work_y,work_bx,work_by,work_candidate_x,work_candidate_y); % Dilation.
        ind2 = sub2ind(dims,x2,y2);
        val2 = Z(ind2); % Level set vals.
        Z(ind2) = -1; % Reset Z back to all -1's.
        Z(ind) = cval1 - 1; % Convolution values - 1.
        cval2 = Z(ind2); % Convolution vals - 1.
        Z(ind2) = -1;
        Z(ind) = cval3 - 1; % Convolution values - 1.
        cval4 = Z(ind2); % Convolution vals - 1.
        Z(ind2) = -1;
        grains{k,1} = ind2;   % Refresh grain's data structure.
        grains{k,2} = val2;   % Ditto.
        grains{k,3} = cval2 + 1; % Ditto.
        grains{k,4} = cval4 + 1; % Ditto.
    end % (for k). Loop over grains ends.
    
 
end

