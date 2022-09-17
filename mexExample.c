#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

#define PI acos(-1.0)


// Example mex function 
// Take orientation of two grains 
// and provide energy using the linear interpolation

double surfacetension(double id1, double id2, double* covarianceData);
double linear_interpolate(double x, unsigned int dataNumber, double *x_k, double *y_k);

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
  double id1; double id2; 
  id1 = mxGetScalar(prhs[0]); 
  id2 = mxGetScalar(prhs[1]); 
  
  // Prepare variables to bring input matrix 
  
  double *a; 
  const mwSize *a_dims;
  int a_rows, a_cols; 
  
  //The structure of a is as follow 
  
  /* misorientation angle  -- energy
    0 0, 0.20
    0.5, 0.25229  
    1.0, 0.31
    1.5, 0.34
    2.0, 0.39
    2.5, 0.43
    3.0, 0.47
    3.5, 0.49
    4.0, 0.50
   */
  
  a = mxGetPr(prhs[2]); 
  a_dims = mxGetDimensions(prhs[2]);
  a_rows = (int) a_dims[0];
  a_cols = (int) a_dims[1];

  
  mexPrintf("Size of input is  %d rows and %d colums", a_rows,a_cols);
  
  for (int i = 0; i < a_rows; i++) {
    mexPrintf("angle: %03f  energy: %03f \n", a[i], a[a_rows+i]);
   }
  

   // Output array:
   // Two types of pointer is needed 
   mxArray *b_out;
   double *b; 
  
   // The outsize is defined here 
   b_out = plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
   b = mxGetPr(b_out);
   
   double misOrient= fabs(id1-id2); 
   double *x_k = &a[0]; 
   double *y_k = &a[a_rows]; 
   
   b[0] = linear_interpolate(misOrient,a_rows, x_k,y_k); 
  
  return; 
}



double linear_interpolate(double x, unsigned int dataNumber, double *x_k, double *y_k )
{
    //I assume the x_k are sorted in a increaseing order
    //Find n, such that x_k[n], x_k[n+1] that includes x    
    unsigned int n=WINT_MAX;
    for(int i=0; i<dataNumber; i++)
    {
        if(x_k[i] <= x && x_k[i+1]> x)
        n=i;
    }
   
    if(n==WINT_MAX){  // if not found, apply periodic data
       x = x - x_k[dataNumber-1];  
      // search again
      for(int i=0; i<dataNumber; i++)
     {
        if(x_k[i] <= x && x_k[i+1]> x)
        n=i;
     }
     
    }
  
    if(n==WINT_MAX){ //still not found, terminate program
     exit(1);
    }
    double slope=(y_k[n+1]-y_k[n])/(x_k[n+1]-x_k[n]);
    double value = y_k[n] + slope * (x - x_k[n]);
    return value;
}

