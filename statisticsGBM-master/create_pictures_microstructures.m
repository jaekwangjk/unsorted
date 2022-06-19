%% 2D microstructures
clear
load('simulation2d-4096x4096-1.mat','dims')
load('simulation2d-4096x4096-1.mat','grains_t0','ori_t0')
showgrainswithori2d(grains_t0,dims,ori_t0)
set(gca,'visible','off')
print(gcf, '-depsc2', './Microstructures/microstructure2d_t0');


load('simulation2d-4096x4096-1.mat','grains_t1','ori_t1')
showgrainswithori2d(grains_t1,dims,ori_t1)
set(gca,'visible','off')
print(gcf, '-depsc2', './Microstructures/microstructure2d_t1');

load('simulation2d-4096x4096-1.mat','grains_t2','ori_t2')
showgrainswithori2d(grains_t2,dims,ori_t2)
set(gca,'visible','off')
print(gcf, '-depsc2', './Microstructures/microstructure2d_t2');

shownetwork2d(grains_t2,dims)
set(gca,'visible','off')
axis([1 1024 1 1024])
print(gcf, '-depsc2', './Microstructures/microstructure2d_t2_zoomin');

%% 2D microstructures (constant)
clear
load('simulation2dconstant-4096x4096-1.mat','dims')
load('simulation2dconstant-4096x4096-1.mat','grains_t0_constant')
shownetwork2d(grains_t0_constant,dims);
set(gca,'visible','off')
print(gcf, '-depsc2', './Microstructures/microstructure2dconstant_t0');


load('simulation2dconstant-4096x4096-1.mat','grains_t1_constant')
shownetwork2d(grains_t1_constant,dims);
set(gca,'visible','off')
print(gcf, '-depsc2', './Microstructures/microstructure2dconstant_t1');

load('simulation2dconstant-4096x4096-1.mat','grains_t2_constant')
shownetwork2d(grains_t2_constant,dims);
set(gca,'visible','off')
print(gcf, '-depsc2', './Microstructures/microstructure2dconstant_t2');

axis([1 1024 1 1024])
print(gcf, '-depsc2', './Microstructures/microstructure2dconstant_t2_zoomin');

%% 3D microstructures

clear
load('simulation3d-1-256x256x256-1.mat','dims')
load('simulation3d-1-256x256x256-1.mat','grains_t0','grains_t1','grains_t2')
load('simulation3d-1-256x256x256-1.mat','ID_t1','ID_t2')
load('simulation3d-1-256x256x256-1.mat','ori_t0')

index_reference = ID_t2(randi(size(grains_t2,1)));

showgrainswithori3d(grains_t0,dims,ori_t0,index_reference)
set(gca,'visible','off')
print(gcf, '-depsc2', './Microstructures/microstructure3d_t0');
set(gca,'visible','off')
showgrainswithori3d(grains_t1,dims,ori_t0,index_reference,ID_t1)
print(gcf, '-depsc2', './Microstructures/microstructure3d_t1');
set(gca,'visible','off')
showgrainswithori3d(grains_t2,dims,ori_t0,index_reference,ID_t2)
print(gcf, '-depsc2', './Microstructures/microstructure3d_t2');
