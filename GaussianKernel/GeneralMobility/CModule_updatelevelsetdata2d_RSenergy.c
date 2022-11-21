#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

#define PI acos(-1.0)

double surfacetension(double ori1, double ori2, double angBrandon);

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
          st = surfacetension(ori[idk-1],ori[idell-1],angBrandon);
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

double surfacetension(double ori1, double ori2,double angBrandon)
{
  double ang1, minang, st, angBrandonRad;
  angBrandonRad = angBrandon*PI/180;

  if (ori2>ori1) ang1 = 2*PI - ori2 + ori1;
  else ang1 = 2*PI - ori1 + ori2;
  minang = min( ang1, fabs(ori1-ori2) );
  
  /* Read-Shockley */
  
  if (minang>angBrandonRad) st = 1;
  else st = minang / angBrandonRad * ( 1-log(minang/angBrandonRad) );
  
  
  return st;
}