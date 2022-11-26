#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

#define PI acos(-1.0)


double surfacetension(int id1, int id2);


void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
 /* 
  *  testCfunction(id1,id2,covarienceEnergy);
  *
  *  CAUTION: MODIFIES MATLAB INPUT *IN PLACE*.
  */

 /*
  mxArray *indices, *grainlevsetvals, *grainconvvals, *locs;
  double *pindices, *plocs, *pgrainlevsetvals, *pgrainconvvals, *ori, *id;
  double sum, mink, st, angBrandon;
  int N,i,j,k,ell,dims,n,nograins,gind,gind2,idk,idell;
  double temp[100], phi[100], minphi[100];
  
  dims = mxGetM(prhs[0]); /* Number of pixels. */
  //n = (int) sqrt(dims);   /* Dimension of grid. */
  //N = mxGetM(prhs[1]);    /* Number of grains. */  
  // ID is just simply grain labels 
  // Here is just copying, using pointers
  //id = (double *) mxGetData(prhs[2]);  /* List of grain IDs. */
  //ori = (double *) mxGetData(prhs[3]); /* Grain orientations. */
  //angBrandon = mxGetScalar(prhs[4]);

  double *covarianceData; 
  double id1; double id2; 
  id1 = mxGetScalar(prhs[0]); 
  id2 = mxGetScalar(prhs[1]); 

  covarianceData = 

}


double surfacetension(int id1, int id2)
{
  double st = 0.5; 
  return st;
}
