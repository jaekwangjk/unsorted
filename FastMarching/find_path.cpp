#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<numeric>
#include<limits>

#include<array> 
#include<vector> 

#define MAXVAL std::numeric_limits<double>::max()
#define LARGEVAL 99998

// Things to do : kill periodicity 

// Create a grid
// Define a function chi on the grid
// Input: grid size, velocity field, chi describing a surface, step size
// Output: chi (t)

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
  double finverse = 1.0/ velocity[index]; 
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


int initial_profile_new(const std::array<double,3>& x)
{
    double func =  20 * (x[0]-0.5)*(x[0]-0.5) + 30 * (x[1]-0.5)*(x[1]-0.5) - 1 ;
    return (func <= 0.) ? 1 : 0;
}



// this will be the main model of the software 
void read_chi(int nx, int ny, int nz, double *chi, std::string fileName){
	
	std::ifstream chiFile;
	
	chiFile.open(fileName);
	
	int pcount = nx*ny*nz;
	double trash; 
	
	int indx=0; 
	if ( chiFile.is_open() ) {
		while (indx<pcount){
		
		chiFile >> trash; // x position 
		std::cout << trash << std::endl;   
		chiFile >> trash; // y position 
		std::cout << trash << std::endl; 
		chiFile >> chi[indx];
		std::cout << chi[indx] << std::endl; 
		indx+=1;
		}
    }
	std::cout << "File has been read succesfully..." << std::endl;
	
	std::cout << "Checking data..." << std::endl;
	
    for(int iz= 0; iz< nz; iz++){ 
      for(int iy= 0; iy< ny; iy++){
		  for(int ix =0; ix < nx; ix++){
	   
		  int indx = iz*nx*ny+ iy*nx + ix; 
		  std::cout << chi[indx] << std::endl;
			 
	}}}
    
	chiFile.close();
}



int main()
{
  const int dim = 3;
  const int numOfIterations = 1e6;
  int nx,ny,nz,indx;
  double h;

  nx = 100; ny = 100; nz = 1;
  h  = 1.0/nx;

  // total number of grid points
  int pcount = nx*ny*nz;
  double *chi; 
  chi = new double[pcount]();
  std::cout << "read chi start" << std::endl; 
  read_chi(nx,ny,nz,chi,"closeChi.txt"); 
  
  // Movie
  /*
  std::string str = "ffmpeg -y -f rawvideo -vcodec rawvideo -pix_fmt gray -s " + std::to_string(nx) + 
                      "x" + std::to_string(ny) + " -r 30 -i - -f mp4 -q:v 5 -an -vcodec mpeg4 out.mp4";
	
  FILE* pipeout = popen(str.c_str(), "w");
  unsigned char* pixels = new unsigned char[pcount];
    	

  fflush(pipeout);
  pclose(pipeout);
  delete[] pixels;
  */  
  return 0;
}


