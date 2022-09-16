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
  *  updatelevelsetdata2doriginal(presence,grains,ID,ORI,angBrandon);
  *
  *  CAUTION: MODIFIES MATLAB INPUT *IN PLACE*.
  */

  mxArray *indices, *grainlevsetvals, *grainconvvals, *locs;
  double *pindices, *plocs, *pgrainlevsetvals, *pgrainconvvals, *ori, *id;
  double sum, mink, st, angBrandon;
  int N,i,j,k,ell,dims,n,nograins,gind,gind2,idk,idell;
  double temp[100], phi[100], minphi[100];
  
  dims = mxGetM(prhs[0]); /* Number of pixels. */
  n = (int) sqrt(dims);   /* Dimension of grid. */
  N = mxGetM(prhs[1]);    /* Number of grains. */  
  // ID is just simply grain labels 
  // Here is just copying, using pointers
  id = (double *) mxGetData(prhs[2]);  /* List of grain IDs. */
  ori = (double *) mxGetData(prhs[3]); /* Grain orientations. */
  angBrandon = mxGetScalar(prhs[4]);

  for (j=0;j<dims;j++){ /* Loop over pixels. */
    //presence가 넘어오는데 
    indices = mxGetCell(prhs[0],j); /* Grains near this pixel. */
    pindices = (double *) mxGetData(indices);
    locs = mxGetCell(prhs[0],2*dims+j); /* Location of pixel in grain's data. */
    plocs = (double *) mxGetData(locs);
    nograins = mxGetN(indices); /* Number of grains near this pixel. */
    
    for (k=0;k<nograins;k++){ /* Loop over (near grains) */
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      i = (int) plocs[k]-1; /* Location of pixel within grain's data. */
      grainconvvals = mxGetCell(prhs[1],2*N+gind-1);
      pgrainconvvals = (double *) mxGetData(grainconvvals);
      temp[k] = pgrainconvvals[i];
    }

    /* These lines implement the redistribution step in Esedoglu-Otto algorithm: */
    /* Form the "phi" functions: */
    
    // work again from here !!!! 
    
    for (k=0;k<nograins;k++){
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      sum = 0.0;
      idk = (int) id[gind-1]; /* id of the k-th grain in the local list. */
      for (ell=0;ell<nograins;ell++){
        if (ell != k) {
          gind = (int) pindices[ell]; /* Index of grain in list of all grains. */
          idell = (int) id[gind-1]; /* id of the ell-th grain in the local list. */
          //st = surfacetension(ori[idk-1],ori[idell-1],angBrandon);
          int id1 = id[idk-1]; 
          int id2 = id[idell-1]; 
          
          st = surfacetension(id1,id2); 
          
          sum = sum + st*temp[ell];
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
      minphi[k] = mink;  // [단위] convolution 값들의 합
    }    
     //min k를 찾았다. 
    /* Form the level set functions: */
    for (k=0;k<nograins;k++){
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      i = (int) plocs[k]-1; /* Location of pixel within grain's data. */
      grainlevsetvals = mxGetCell(prhs[1],N+gind-1);
      pgrainlevsetvals = (double *) mxGetData(grainlevsetvals);

      pgrainlevsetvals[i] = (minphi[k] - phi[k]); // if phi[k]=minphi[k] this is 0, other wise, this is negative...
      if (nograins==1) pgrainlevsetvals[i] = temp[k]; // this is convolution value 야
    }
    
  } /* (for j). Loop over pixels ends. */
    
}


double surfacetension(int id1, int id2)
{
  int id1Type, id2Type; 
  
  // id1%3 =0 - type A 
  // id2%3 =1 - type B
  // id1%3 =2 - type C
  
  if(id1%3==0)
    id1Type=0; 
  else if( id1%3==1)
    id1Type=1; 
  else
    id1Type=2; 
   
  if(id2%3==0)
    id2Type=0; 
  else if(id2%3==1)
    id2Type=1; 
  else
    id2Type=2; 
  
  
  double st; 
  double gamma1=1.0; 
  double gamma2=1.1;
  double gamma3=1.3; 
  
  
  if(id1Type==id2Type)
      st = gamma1; 
  else if (id1Type== 2 || id2Type==2)
      st = gamma2; 
  else 
      st = gamma3; 
  
  return st;
}