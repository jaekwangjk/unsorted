#include <iostream>

// absolute path for eigen
#include </Users/jaekwangkim/Program/eigen/Eigen/Dense>
#include <memory>
#include <fftw3.h>
#include <fstream>
//C++ Libraries
#include <algorithm>
#include <vector>

#include <cstdio> // sprintf 

#define PI 3.14159265359

using namespace std;

int main (); 
void MyGaussianKernel (int totalTimeStep);


void MyGaussianKernel (int totalTimeStep)
{
  std::cout <<"MyGaussianKernel function begins..." << std::endl;
	
	// Discretization
  int nx = 1024;
  int ny = 1024;

	// Define Surface Tension Matrix \sigma_{ij}
	// For now, it is hard coded
	
	// [a] Diagonal terms are zero, and [b] it should be symmetric
	// [c] Trigonometric inequaltiy should be considred.
	double sigma[3][3]; 
	sigma[0][0]=0.0; sigma[1][1]=0.0; sigma[2][2]=0.0; 
	
	sigma[0][1]=0.62; sigma[1][0]=0.62;
	sigma[0][2]=0.87; sigma[2][0]=0.87;
	sigma[1][2]=0.62; sigma[2][1]=0.62;
  
  //Step[2] Generate grid
	double Lx=1; double Ly=1; 		
    
	double *x_grid; double *y_grid; 
	x_grid=(double *) malloc(sizeof(double)*nx); 
	y_grid=(double *) malloc(sizeof(double)*ny); 

	for (int j = 0; j < nx; j++ ){
    x_grid[j] = 0.0+ Lx*j/(nx-1);
    y_grid[j] = 0.0+ Ly*j/(ny-1);
  }
	
	//Step[3] Generate Initial Condition
	
	fftw_complex *chi_1;
	fftw_complex *chi_2;
	fftw_complex *chi_3;
	
	chi_1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny);
	chi_2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny);
	chi_3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny);
   	
	//Convention of current 2D mesh grid is to increase x first, as a colunm,
	//then increase y as a row
	
  //Initialization
  double xPoint,yPoint;
 //loop folded
  for (int j = 0; j < ny; j++ ){
   for (int k = 0; k < nx; k++ ){
		   
     xPoint = x_grid[k]; yPoint = y_grid[j];
		 
     if(yPoint>0.8 || yPoint<0.2){
       //Characteristic Function 2
       chi_1[j*nx+k][0]=0.0;  /* complex part is zero */ chi_1[j*nx+k][1]=0.0;
       chi_2[j*nx+k][0]=1.0;  /* complex part is zero */ chi_2[j*nx+k][1]=0.0;
       chi_3[j*nx+k][0]=0.0;  /* complex part is zero */ chi_3[j*nx+k][1]=0.0;
     }else{
       if(xPoint<0.25||xPoint>0.75){
        //Characteristic Function 1
       chi_1[j*nx+k][0]=1.0;  /* complex part is zero */ chi_1[j*nx+k][1]=0.0;
       chi_2[j*nx+k][0]=0.0;  /* complex part is zero */ chi_2[j*nx+k][1]=0.0;
       chi_3[j*nx+k][0]=0.0;  /* complex part is zero */ chi_3[j*nx+k][1]=0.0;
     }else{
       chi_1[j*nx+k][0]=0.0;  /* complex part is zero */ chi_1[j*nx+k][1]=0.0;
       chi_2[j*nx+k][0]=0.0;  /* complex part is zero */ chi_2[j*nx+k][1]=0.0;
       chi_3[j*nx+k][0]=1.0;  /* complex part is zero */ chi_3[j*nx+k][1]=0.0;
      }
     }
    } //loop over k
 }//loop over j
 
 ofstream fout;

 // file out, initial condition
/*
fout.open("output/chi_000.txt",ios::out);
for (int j = 0; j < ny; j++ ){
  for (int k = 0; k < nx; k++ ){
		  fout << chi_1[j*nx+k][0] + 2.0 * chi_2[j*nx+k][0] + 3.0 * chi_3[j*nx+k][0] << " " ;
   }
  fout << std::endl;
 }
	fout.close(); 
*/
  
	//[2] FFT & THRESHOLDING & TIME MARCHING
	
	//FQ of each characteristic function 
	
	fftw_complex *FQ_chi_1;
	fftw_complex *FQ_chi_2;
	fftw_complex *FQ_chi_3;
	
	FQ_chi_1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny); 
	FQ_chi_2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny); 
	FQ_chi_3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny); 
	
	//Diffusion considering non-equal surface tension 
	//Will be determined by summation over j for  FQ_phi_i= (sigma_{ij} * FQ_chi_j) 
	
	fftw_complex *FQ_phi_1;
	fftw_complex *FQ_phi_2;
	fftw_complex *FQ_phi_3;
	
	FQ_phi_1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny); 
	FQ_phi_2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny); 
	FQ_phi_3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx * ny); 
	
	// phi_i will be constructed by IFFT of FQ_phi_i
	// and will be used as thresholding criteria
	fftw_complex *phi_1;
	fftw_complex *phi_2;
	fftw_complex *phi_3;
	
	phi_1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny);
	phi_2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny);
	phi_3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny);
	
	double dt =100.0;
	double D = 1;  
	
	for (int timeStep=1; timeStep < totalTimeStep; timeStep++){
		
		//[2]-[a] : FFT of characteristic functions 
		fftw_plan my_fft_plan; 
     
		my_fft_plan = fftw_plan_dft_2d (nx,ny, chi_1, FQ_chi_1, FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute ( my_fft_plan );
		
    my_fft_plan = fftw_plan_dft_2d (nx,ny, chi_2, FQ_chi_2, FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute ( my_fft_plan );
	
  	my_fft_plan = fftw_plan_dft_2d (nx,ny, chi_3, FQ_chi_3, FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute ( my_fft_plan );
		
    int j_star; int k_star;
		
    //[2]-[b] : Operate Gaussian Kernel
    for (int j = 0; j < ny; j++ ){
      for (int k = 0; k < nx; k++ ){
      
        if(j <ny/2)
          j_star=j;
        else
          j_star = j-ny;
			   
        if(k <nx/2)
          k_star=k;
        else
          k_star= k-ny;
           
       //real part & complex part
       FQ_chi_1[j*nx+k][0] *= exp(-4.0*dt*D*pow(PI,2)* ( pow(j_star,2)+pow(k_star,2) ) / (nx*ny) );
       FQ_chi_1[j*nx+k][1] *= exp(-4.0*dt*D*pow(PI,2)* ( pow(j_star,2)+pow(k_star,2) ) / (nx*ny) );
			   
       FQ_chi_2[j*nx+k][0] *= exp(-4.0*dt*D*pow(PI,2)* ( pow(j_star,2)+pow(k_star,2) ) / (nx*ny) );
       FQ_chi_2[j*nx+k][1] *= exp(-4.0*dt*D*pow(PI,2)* ( pow(j_star,2)+pow(k_star,2) ) / (nx*ny) );
			   
       FQ_chi_3[j*nx+k][0] *= exp(-4.0*dt*D*pow(PI,2)* ( pow(j_star,2)+pow(k_star,2) ) / (nx*ny) );
       FQ_chi_3[j*nx+k][1] *= exp(-4.0*dt*D*pow(PI,2)* ( pow(j_star,2)+pow(k_star,2) ) / (nx*ny) );
        
       //assemble FQ_phi_i, be carefure indexing in sigma
       FQ_phi_1[j*nx+k][0] = sigma[0][1] * FQ_chi_2[j*nx+k][0]+ sigma[0][2] * FQ_chi_3[j*nx+k][0];
       FQ_phi_1[j*nx+k][1] = sigma[0][1] * FQ_chi_2[j*nx+k][1]+ sigma[0][2] * FQ_chi_3[j*nx+k][1];
				    
       FQ_phi_2[j*nx+k][0] = sigma[1][0] * FQ_chi_1[j*nx+k][0]+ sigma[1][2] * FQ_chi_3[j*nx+k][0];
       FQ_phi_2[j*nx+k][1] = sigma[1][0] * FQ_chi_1[j*nx+k][1]+ sigma[1][2] * FQ_chi_3[j*nx+k][1];
			   
       FQ_phi_3[j*nx+k][0] = sigma[2][0] * FQ_chi_1[j*nx+k][0]+ sigma[2][1] * FQ_chi_2[j*nx+k][0];
       FQ_phi_3[j*nx+k][1] = sigma[2][0] * FQ_chi_1[j*nx+k][1]+ sigma[2][1] * FQ_chi_2[j*nx+k][1];
    }
   }
		

  //[2]-[c]: take inverse FFT
  fftw_plan my_ifft_plan;
	
  my_ifft_plan = fftw_plan_dft_2d (nx,ny,FQ_phi_1,phi_1,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute ( my_ifft_plan );
	
  my_ifft_plan = fftw_plan_dft_2d (nx,ny,FQ_phi_2,phi_2,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute ( my_ifft_plan );
	
  my_ifft_plan = fftw_plan_dft_2d (nx,ny,FQ_phi_3,phi_3,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute ( my_ifft_plan );
	
	
  //[2]-[d] Thresholding
  // Choose minimum phi
		 
  for (int j = 0; j < ny; j++ ){
    for (int k = 0; k < nx; k++ ){
 	     		
    phi_1[j*nx+k][0] *= 1.0/(nx*ny);
    phi_2[j*nx+k][0] *= 1.0/(nx*ny);
    phi_3[j*nx+k][0] *= 1.0/(nx*ny);
				
      if (phi_1[j*nx+k][0] <= phi_2[j*nx+k][0] ){
        if(phi_1[j*nx+k][0] <= phi_3[j*nx+k][0]){
          chi_1[j*nx+k][0]=1.0; chi_2[j*nx+k][0]=0.0; chi_3[j*nx+k][0]=0.0;
      }else{
        chi_1[j*nx+k][0]=0.0; chi_2[j*nx+k][0]=0.0; chi_3[j*nx+k][0]=1.0;
      }}
      else{
        if (phi_2[j*nx+k][0] <= phi_3[j*nx+k][0]){
          chi_1[j*nx+k][0]=0.0; chi_2[j*nx+k][0]=1.0; chi_3[j*nx+k][0]=0.0;
      }else{
          chi_1[j*nx+k][0]=0.0; chi_2[j*nx+k][0]=0.0; chi_3[j*nx+k][0]=1.0;
       }
      }
    }// end loop k
  }// end loop j
	
 
  //outputs the result
  if( (timeStep+1)%20==0){
    char filename[50];
    sprintf(filename, "output/chi_%03d.txt", timeStep);
	 	fout.open(filename,ios::out);

    for (int j = 0; j < ny; j++){
      for (int k = 0; k < nx; k++){
        fout << chi_1[j*nx+k][0] + 2.0 * chi_2[j*nx+k][0] + 3.0 * chi_3[j*nx+k][0] << " " ;
      }
      fout << std::endl;
    }

    fout.close();
  }
	
    cout << "Running with timestep " << timeStep << endl;
	
	}//end of time step, loop over "for (timeStep=0; timeStep < 1; timeStep++)"
	
  return;
}


int main ()
{
	MyGaussianKernel (200);

	std::cout <<"Program has ended sucessfully" <<std::endl;		
    return 0; 
}

