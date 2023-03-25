#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<numeric>
#include<limits>
#include "simple_heap.h"
#include "Matrix.h"
#include<array> 
#include<vector> 
#include <fftw3.h>

#define MAXVAL std::numeric_limits<double>::max()
#define LARGEVAL 99998


// To be done... initial profile diffusion [done]

// wall, visited, chi(grid), velocity로 변수 구분시키고 
// 데이터 아웃 포맷 만들기 (MATLAB으로 즉각즉각 확인하게 좋도록...)
// Fast Marching Paper 다시한번 읽어보기..

// [2] velocity field 
// [3] ...
// Create a grid
// Define a function chi on the grid
// Input: grid size, velocity field, chi describing a surface, step size
// Output: chi (t)
// Create a grid
// Define a function chi on the grid
// Define a surface, by initializing chi on all outer and inner grid points (w.r.t surface) to 0 and 1 respectively
// (J.K.??????)
// Tag the grid points as "alive", "close" (dead points that share a neighbor with alive points) and "far" points
// Initialize u to zero at all alive points
// loop 
//     Find the value of u at all close points by solving the quradratic equation
//     Heapify close points with back pointers to the grid
//     move the grid associated to the root to "alive"
//     identify the "neighbors of the root" that are in "far", and "move" to "close" to form the new "close"
// end

// Input: grid size, velocity field, chi describing a surface, step size
// Output: chi (t)
double eikonalUpdate_2d(double h1s, double h2s, double U1, double U2, double finverse); 
double eikonalUpdate(int nx, int ny, int nz, int x, int y, int z, double *grid, double *velocity); 

//Poisson equation FFT solver struct
typedef struct  {
  fftw_plan Forward;
  fftw_plan Backward;
  double *kernel;

  fftw_complex *workspace;
  
}poissonSolver;

//double initial_profile_new(const std::array<double,3>& x)
//{
//    double func =  20 * (x[0]-0.5)*(x[0]-0.5) + 30 * (x[1]-0.5)*(x[1]-0.5) - 1 ;
//    return func; 
//}

void initialize_wall(int nx, int ny, int nz, double *wall, double *visited, const std::array<double,3>& corner, const double& h ){
    std::array<double,3> x;
	
	int pcount = nx*ny*nz; 
		
    std::ofstream wallFile;
    wallFile.open("initialWall.txt");
	
    // Two step ---do profile first do other step later. 
    for(int iz= 0; iz< nz; iz++){ 
	  x[2] = corner[2]+ iz*h;	
      for(int iy= 0; iy< ny; iy++){
		x[1] = corner[1]+ iy*h ;  
		  for(int ix =0; ix < nx; ix++){
		    x[0] = corner[0]+ ix*h ; 
			 
			int indx = iz*nx*ny+ iy*nx + ix; 
			double func =  20 * (x[0]-0.5)*(x[0]-0.5) + 30 * (x[1]-0.5)*(x[1]-0.5) - 1 ;
			
			if(func < 0 ) 
			object[indx] = 1.0; 
			else
			object[indx] = 0.0; 
			
			wallFile << x[0] << ", " << x[1]  << ", " << chi[indx] << std::endl;
	}}}
	
    for(int i=0; i< pcount; i++){	
	  visited[i]=0.0; //initialize 
		 
	  if( fabs(wall[indx]-1.0) < tol)
	   visited[indx]=1.0; 	 
	}  
			 
    wallFile.close();
}


// this will be the main model of the software 
void initialize_velocity(int nx, int ny, int nz, double *chi, double *velocity){
	
	double smallValue = 1e-9;
    
	// smoothing obstacles
	poissonSolver fftps; /* The"workspace"-stucture */
	int pcount = nx*ny*nz; 
	fftps.workspace= new fftw_complex[pcount]();
	fftps.kernel = new double[pcount]();
	fftps.Forward=fftw_plan_dft_3d(nz,ny,nx,fftps.workspace,fftps.workspace,
		                                       -1,FFTW_MEASURE);
	fftps.Backward=fftw_plan_dft_3d(nz,ny,nx,fftps.workspace,fftps.workspace,
		                                        1, FFTW_MEASURE);
	// prepare kernel...
	for(int iz= 0; iz< nz; iz++){ 
	 for(int iy= 0; iy< ny; iy++){
	  for(int ix =0; ix < nx; ix++){
	    double xk = 2*M_PI*ix/(nx*1.0);
	    double yk = 2*M_PI*iy/(ny*1.0);
		double zk = 2*M_PI*iz/(nz*1.0);
		double Eigenvalues=2*nx*nx*(1-cos(xk))+2*ny*ny*(1-cos(yk))+2*nz*nz*(1-cos(zk));
		fftps.kernel[iz*ny*nx+ iy*nx+ix]=Eigenvalues;
    }}}
	
	std::cout<<"I am here! A4 !" << std::endl; 

    std::cout<<"I am here! B !" << std::endl; 
	
	for(int i=0;i<pcount;i++){
		fftps.workspace[i][0]=chi[i];
		fftps.workspace[i][1]=0.0;
	}	
	
	fftw_execute(fftps.Forward);
	
	std::cout<<"I am here! C !" << std::endl; 

	double dt = 0.01; 
	
	for(int i=0;i<pcount;i++){
	    fftps.workspace[i][0]*=exp(-fftps.kernel[i]*dt)/pcount;
	    fftps.workspace[i][1]*=exp(-fftps.kernel[i]*dt)/pcount;
	}	
	
	fftw_execute(fftps.Backward);
	
	std::cout<<"I am here! D !" << std::endl; 
	
	for(int i=0;i<pcount;i++){
		chi[i]=fftps.workspace[i][0]; 
	}	
	
	
	
	for(int iz= 0; iz< nz; iz++){ 
      for(int iy= 0; iy< ny; iy++){ 
		  for(int ix =0; ix < nx; ix++){
			  int indx = iz*nx*ny+ iy*nx + ix; 
			  
			  if(chi[indx] < 0.0)
				  velocity[indx] = smallValue; 
			  else
				  velocity[indx] = 1.0; 
   }}}
}



int main()
{
    typedef double data_t;

    const int grayscale = 255;
	const int grayscale_obs = 127;
    const int dim = 3;
    const int numOfIterations = 1e6;
    int nx,ny,nz,indx;
    const std::array<data_t,dim> corner{0.,0.,0.};
    
	std::array<int,dim> grid_indx;
    std::array<double,dim> x;
    data_t h;

    nx = 100; ny = 100; nz = 1;
    h  = 1.0/nx;

    // total number of grid points
    int pcount = nx*ny*nz;

	double *chi; 
	double *visited; 
	double *velocity; 
	
	
	chi = new double[pcount]();
	visited = new double[pcount]();
	velocity = new double[pcount]();
	
	
	// This is input!! 
	//Matrix<dim> velocity(nx,ny,nz);
    initialize_grid(nx,ny,nz, chi,visited,corner,h);
	initialize_velocity(nx,ny,nz,chi,velocity);
	
    // For Movie
    //std::string str = "ffmpeg -y -f rawvideo -vcodec rawvideo -pix_fmt gray -s " + std::to_string(nx) + 
	//                  "x" + std::to_string(ny) + " -r 30 -i - -f mp4 -q:v 5 -an -vcodec mpeg4 out.mp4";
	
    //FILE* pipeout = popen(str.c_str(), "w");
    unsigned char* pixels = new unsigned char[pcount];
    
	
	i_heap heap;
	heap = i_heap_create_empty_heap(pcount, pcount);
	i_heap_clear_heap(&heap);
	
    clock_t b,e;
    b=clock();
	
    // Add all boundary points to the heap
    double tol =1e-7;
	
    for(int iz= 0; iz< nz; iz++){ 
      for(int iy= 0; iy< ny; iy++){
		for(int ix =0; ix < nx; ix++){
					
        int indx = iz*nx*ny+ iy*nx + ix; 
		
        if ( fabs(chi[indx] -0.0) < tol){
	      i_heap_insert_node_with_location(&heap, indx, chi[indx], indx);
          // the latter 'indx' is for the location pointer...
          pixels[indx] = grayscale;
        }
		
		if (chi[indx] < 0.0){
		  pixels[indx] = grayscale_obs; 
		}
    }}}
	
	// At this point the heap contains all the boundary points
 
    double chi_east,chi_west,chi_north,chi_south,chi_top,chi_bottom;
    double chi_x,chi_y,chi_z;
    double val;
	
    // Main Iteration Loop begins 
    //     Extract the root, i.e remove it from heap
    //     Find its neighbors that are in the interior
    //     Compute chi at all interior neighbors of the extracted root and add them to the heap
    // 	   Repeat Above
	
	
	int xgrid[6]={1, 0, -1,  0, 0, 0};
	int ygrid[6]={0, 1, 0 , -1, 0, 0};
	int zgrid[6]={0, 0, 0 ,  0, 1,-1};
	
	
	// Mark initial point and visited pixels 
	
    
    int it;
    //for(it=0;it<numOfIterations && close.heap_size!=0 ;it++)
   
    while (!i_heap_empty(&heap))
	{   //Previous Declation MinHeap<data_t> close(pcount);
        
		//use the root  
		int *currentIndex=&heap.root[0].originalIndex;
		double *min=&heap.root[0].key;
		
	    //index(i,j,k) is descrbied as "i*n2*n1+j*n1+k"
	    //current point
	    int cx = *currentIndex % nx;
	    int cy = (*currentIndex % (nx*ny)) / nx;
	    int cz = (*currentIndex)/ (nx*ny);
		
	    for(int m=0; m<6; m++)  {
    
	      int x=cx+xgrid[m];
	      int y=cy+ygrid[m];
	      int z=cz+zgrid[m];
		  
		  // this will result in Neumann BC for chi 
	      if(x>nx-1){x=nx-1;} if(x<0){x=0;}
	      if(y>ny-1){y=ny-1;} if(y<0){y=0;}
	      if(z>nz-1){z=nz-1;} if(z<0){z=0;}
		  
		  int index= z*ny*nx + y*nx + x;
		  
		  if(visited[index] < 0.9) //not equal to 1.0 means it is not visited 
		  {
			double possible = fmax(eikonalUpdate(nx,ny,nz,x,y,z,chi,velocity),*min);
		
  	        if(possible<chi[index]) {
  	         if(chi[index]==LARGEVAL) {
  	          chi[index]=possible;
  	          i_heap_insert_node_with_location(&heap, index, chi[index], index);
  	          }else{
  	           int heapLocation=heap.locations[index];
  	           i_heap_decrease_key(&heap,heapLocation,chi[index]);
  	          }
  	  	     }//end of possible  
		  }	
		}// end of neighbor m   
		
    // mark presented position now as visited and delete the root
	visited[* currentIndex] = 1.0; 
	i_heap_delete_min(&heap);
			  
	}//end of while roop
	
	e=clock();
	std::cout << (e-b)/(CLOCKS_PER_SEC*1.0) << std::endl;;	


    std::ofstream closeFile;
    closeFile.open("closeChi.txt");
 	
    for(int iz= 0; iz< nz; iz++){ 
	  x[2] = corner[2]+ iz*h;	
      for(int iy= 0; iy< ny; iy++){
		x[1] = corner[1]+ iy*h ;  
		  for(int ix =0; ix < nx; ix++){
			 x[0] = corner[0]+ ix*h ; 	
			 int indx = iz*nx*ny+ iy*nx + ix; 
			 closeFile << x[0] << " " << x[1]  << " " << chi[indx] << std::endl;	  	  
	}}}
	
   //fflush(pipeout);
   //pclose(pipeout);
   //delete[] pixels;
    
   return 0;
}



// detailed functions 

double eikonalUpdate_2d(double h1s, double h2s, double U1, double U2, double finverse)
{
  // conduct 2d eikonalUpdate, when discrim <0
  // h1s,h2s: h1^2 where h1= 1.0/n1;
  // U1,U2 : minimum value in 1 and 2 direction
    
  double result=0.0;
  double alpha = h1s/h2s;
  
  double a = h1s + h2s;
  double b = -2.0 * (h1s*U1 + h2s*U2);
  double c = h1s*U1*U1 + h2s*U2*U2 - h1s*h2s*pow(finverse,2);      
  double discrim = b*b - 4*a*c;
    
  if(discrim>0)  {
    result = (sqrt(b*b-4*a*c)-b)/(2*a);
  }else if (U2 < U1)  {
    result = U2 + sqrt(h2s) * finverse;
  }else  {
    result = U1 + sqrt(h1s) * finverse;
  }
  
  return result;
}

double eikonalUpdate(int nx, int ny, int nz, int x, int y, int z, double *grid, double *velocity)
{
  double result =0.0;
  //eikonal update in 3d with non-equal spacing for each direction
  double hx=1/(nx*1.0); double hy=1/(ny*1.0); double hz=1/(nz*1.0);
  
  int index = z*ny*nx + y*nx + x;
  
  int x_west,x_east,y_south,y_north,z_down,z_up;
  
  x_west=fmax(x-1,0);
  x_east=fmin(nx-1,x+1);
  y_south=fmax(y-1,0);
  y_north=fmin(ny-1,y+1);
  z_down=fmax(z-1,0);
  z_up=fmin(nz-1,z+1);  
  
  int west= z*ny*nx+ y*nx + x_west;
  int east= z*ny*nx+ y*nx + x_east;
  int south = z*ny*nx + y_south * nx + x;
  int north = z*ny*nx + y_north * nx + x;
  int down = z_down*ny*nx + y*nx + x;
  int up = z_up*ny*nx + y*nx +x;

  //f inverse in slowness
  double finverse = 1.0/velocity[index]; 
  double tol = 1e-7;
  
  //testgrid is grid, true grid is indicator -- this is equivalent to upwind
  double U_H=fmin(grid[west],grid[east]);
  double U_V=fmin(grid[north],grid[south]);
  double U_Z=fmin(grid[up], grid[down]);
    
  //We solve the quadratic equation for U_ij, current cell, using up wind
  double hxs = hx*hx;
  double hys = hy*hy;
  double hzs = hz*hz;

  double a = hxs*hys + hys*hzs + hxs*hzs;
  double b = (-2.0) * (hys*hzs * U_H + hxs*hzs * U_V + hxs*hys* U_Z);
  double c = hys*hzs*pow(U_H,2)+ hxs*hzs*pow(U_V,2) +
                   hzs*hys*pow(U_Z,2)- hxs*hys*hzs*pow(finverse,2);
    
  double discrim = b*b - 4.0 *a *c;
    
  if(discrim>0) {
    result=(sqrt(b*b-4.0*a*c)-b)/(2.0*a);
  }else  {   //if discrim < 0,  then we take 2 dimensional updates in each plane
        // and choose the smallest.
    double result_xy;
    double result_yz;
    double result_zx;
    
    //in x-y plane, y-z, z-x plane separatedly
    result_xy = eikonalUpdate_2d(hxs,hys,U_H,U_V,finverse);
    result_yz = eikonalUpdate_2d(hys,hzs,U_V,U_Z,finverse);
    result_zx = eikonalUpdate_2d(hzs,hxs,U_Z,U_H,finverse);
    
    // result=fmin(result_xy,result_yz);
    // result=fmin(result, result_zx);
    // just simply choose xy results
	
	result=result_xy;
  }
     
     return result;
}