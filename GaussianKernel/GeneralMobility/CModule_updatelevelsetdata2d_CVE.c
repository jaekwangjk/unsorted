#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

#define PI acos(-1.0)

double surfacetension(double ori1, double ori2, unsigned int dataNumber, double *x_k, double *y_k);

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
 /*  
  *  updatelevelsetdata2d(presence,grains,ID,ori,alpha,beta,angBrandon,option);
  *
  *  CAUTION: MODIFIES MATLAB INPUT *IN PLACE*.
  */

  mxArray *indices, *grainlevsetvals, *grainconvvals1, *grainconvvals2, *locs;
  double *pindices, *plocs, *pgrainlevsetvals, *pgrainconvvals1, *pgrainconvvals2, *id, *ori;
  double sum, mink, st, m, a, b, alpha, beta, angBrandon;
  int N,i,j,k,ell,dims,n,nograins,gind,gind2,idk,idell, option;
  double temp1[100], temp2[100], phi[100], minphi[100];
    
  dims = mxGetM(prhs[0]); /* Number of pixels. */
  n = (int) sqrt(dims);   /* Dimension of grid. */
  N = mxGetM(prhs[1]);    /* Number of grains. */  

  // ID is just simply grain labels 
  // Here is just copying, using pointers
  id = (double *) mxGetData(prhs[2]);  /* List of grain IDs. */
  ori = (double *) mxGetData(prhs[3]); /* Grain orientations. */
  alpha = mxGetScalar(prhs[4]);
  beta = mxGetScalar(prhs[5]);
  angBrandon = mxGetScalar(prhs[6]);
  option = (int) mxGetScalar(prhs[7]);
  
  
  // Prepare variables to bring input matrix 
  double *energy; 
  const mwSize *energy_dims; 
  int eg_rows, eg_cols; 
  
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
  
  energy = mxGetPr(prhs[8]); // covariance energy data  
  energy_dims = mxGetDimensions(prhs[8]); 
  eg_rows = (int) energy_dims[0]; 
  eg_cols = (int) energy_dims[1]; 
  double *misorientation_k = &energy[0]; 
  double *energy_k = &energy[eg_rows]; 
  
  
  
  for (j=0;j<dims;j++){ /* Loop over pixels. */
    //presence가 넘어오는데 
    indices = mxGetCell(prhs[0],j); /* Grains near this pixel. */
    pindices = (double *) mxGetData(indices);
    locs = mxGetCell(prhs[0],3*dims+j); /* Location of pixel in grain's data. */
    plocs = (double *) mxGetData(locs);
    nograins = mxGetN(indices); /* Number of grains near this pixel. */
    
    for (k=0;k<nograins;k++){ /* Loop over grains. */
        gind = (int) pindices[k]; /* Index of grain in list of all grains. */
        i = (int) plocs[k]-1; /* Location of pixel within grain's data. */
      grainconvvals1 = mxGetCell(prhs[1],2*N+gind-1);
      pgrainconvvals1 = (double *) mxGetData(grainconvvals1);
      grainconvvals2 = mxGetCell(prhs[1],3*N+gind-1);
      pgrainconvvals2 = (double *) mxGetData(grainconvvals2);
      temp1[k] = pgrainconvvals1[i];
      temp2[k] = pgrainconvvals2[i];
    }

    /* These lines implement the redistribution step in Esedoglu-Otto algorithm: */
    /* Form the "phi" functions: */
    for (k=0;k<nograins;k++){
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      sum = 0.0;
      idk = (int) id[gind-1]; /* id of the k-th grain in the local list. */
      for (ell=0;ell<nograins;ell++){
        if (ell != k) {
          gind = (int) pindices[ell]; /* Index of grain in list of all grains. */
          idell = (int) id[gind-1]; /* id of the ell-th grain in the local list. */
          
          st = surfacetension(ori[idk-1],ori[idell-1], eg_rows, misorientation_k, energy_k); 
          //mexPrintf("orientation angle %d and %d has energy %04f \n", id1,id2, st);
          
          switch (option)
          {
              case 1:{
                  m = 1.0;
                  break;
              }
              case 2:{
                  m  = 1/st;
                  break;
              }
          }
          a = sqrt(PI)*sqrt(alpha)/(alpha-beta)*(st-beta/m);
          b = sqrt(PI)*sqrt(beta)/(alpha-beta)*(-st+alpha/m);
          sum = sum + a*temp1[ell] + b*temp2[ell];
        }
      }
      phi[k] = sum;
    }
    
    /* Minimization over the "phi" functions involved in forming level set functions: */
    for (k=0;k<nograins;k++){
      mink = 1e100;
      for (ell=0;ell<nograins;ell++){
        if (ell != k) {
          mink = min( mink , phi[ell] );
        }
      }
      minphi[k] = mink;
    }    

    /* Form the level set functions: */
    for (k=0;k<nograins;k++){
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      i = (int) plocs[k]-1; /* Location of pixel within grain's data. */
      grainlevsetvals = mxGetCell(prhs[1],N+gind-1);
      pgrainlevsetvals = (double *) mxGetData(grainlevsetvals);

      pgrainlevsetvals[i] = (minphi[k] - phi[k]); // if phi[k]=minphi[k] this is 0, other wise, this is negative...
      if (nograins==1) pgrainlevsetvals[i] = temp1[k]; // this is convolution value 야
    }
    
  } /* (for j). Loop over pixels ends. */
    
}


// Do linear interpolation using the existing data 
double surfacetension(double ori1, double ori2, unsigned int dataNumber, double *x_k, double *y_k)
{
  
  // orientation is random between 0 to 70 
  double misorient = fabs(ori1-ori2); 
  
  //I assume the x_k are sorted in a increaseing order
  //Find n, such that x_k[n], x_k[n+1] that includes x    
  unsigned int n=WINT_MAX;

  for(int i=0; i<dataNumber; i++)
  {
      if(x_k[i] <= misorient && x_k[i+1]> misorient)
      n=i;
  } 
 
  if(n==WINT_MAX){ //still not found, terminate program
   exit(1);
  }
 
 
  double slope=(y_k[n+1]-y_k[n])/(x_k[n+1]-x_k[n]);
  double st= y_k[n] + slope * (misorient - x_k[n]);
  //mexPrintf("orientation angle %f and %f has misorientation %01f  energy %03f \n", ori1,ori2,misorient, st);

  return st;
}

