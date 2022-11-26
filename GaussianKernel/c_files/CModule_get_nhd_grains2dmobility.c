#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

double CVAL[100]; /* The max. number of neighboring grains allowed is 100. */
int INC[100];     /* The max. number of neighboring grains allowed is 100. */

int comparison(int *a, int *b);
void set_INC_and_CVAL(double *cval, int m);

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
  /*
    MATLAB INTERFACE:
      presence = get_nhd_grains2d(grains,dim1*dim2)
    where dim1 and dim2 are the dimensions of the grid.
    The output "presence" is an dim1*dim2-by-4 cell array:
      presence{index,1} = Grains in a nhd. of pixel at loc. index.
      presence{index,2} = Conv. vals. of those grains, in same order.
      presence{index,3} = Conv. vals. of those grains, in same order.
      presence{index,4} = Used to locate this pixel in each grains's 
                          data structure.
    Example:
      If k=presence{index,1}(i) and j=presence{index,4}(i), then
      grains{k,1}(j) = index. ??
  */
  mxArray *gind, *gloc, *gcval1, *gcval2, *cell1, *cell2, *cell3, *cell4;
  mwSize ndim, vecdim[2];
  int i,j,k,m,n,N,dims;
  int **ind, **loc;
  double **cval1, **cval2;
  int *nopresent; /* Number of grains present in a nhd. of each pixel. */
  double *pgind, *pgloc, *pgcval1, *pgcval2; /* C pointer versions of gind, gloc, gcval. */
  int index, location;
  double convolution1, convolution2, *p1, *p2, *p3, *p4;
  int *limit;
  
  N = mxGetM(prhs[0]); /* Number of grains. */
  dims = (int) mxGetScalar(prhs[1]); /* Number of pixels. */

  /* Local, C version of data structure: */
  ind = (int **) calloc(dims,sizeof(int *));
  loc = (int **) calloc(dims,sizeof(int *));
  cval1 = (double **) calloc(dims,sizeof(double *));
  cval2 = (double **) calloc(dims,sizeof(double *));
  nopresent = (int *) calloc(dims,sizeof(int));
  limit = (int *) calloc(dims,sizeof(int));
  
  /* Initial allocation: 5 neighbors. */
  for (i=0;i<dims;i++){
    limit[i] = 5;
    ind[i] = (int *) calloc(limit[i],sizeof(int));
    loc[i] = (int *) calloc(limit[i],sizeof(int));
    cval1[i] = (double *) calloc(limit[i],sizeof(double));
    cval2[i] = (double *) calloc(limit[i],sizeof(double));
  }
  
  /* Loop over grains: */
  for (k=0;k<N;k++){
    gind = mxGetCell(prhs[0],k);
    pgind = (double *) mxGetData(gind);
    gcval1 = mxGetCell(prhs[0],2*N+k);
    pgcval1 = (double *) mxGetData(gcval1);
    gcval2 = mxGetCell(prhs[0],3*N+k);
    pgcval2 = (double *) mxGetData(gcval2);
        
    m = mxGetM(gind); /* Number of pixels in buffered grain */
    
    for (i=0;i<m;i++) { /* Loop over pixels in buffered grain. */
      index = pgind[i]-1;  /* Grid index of current pixel. */
      convolution1 = pgcval1[i];
      convolution2 = pgcval2[i];
      j = nopresent[index]; /* Number of grains already at pixel. */            
      /* Reallocate memory, in chunks, if needed: */
      if (j>=limit[index]) {
        limit[index] = limit[index] + 5;
        ind[index] = realloc( ind[index], limit[index] * sizeof(int) );
        loc[index] = realloc( loc[index], limit[index] * sizeof(int) );
        cval1[index] = realloc( cval1[index], limit[index] * sizeof(double) );
        cval2[index] = realloc( cval2[index], limit[index] * sizeof(double) );
      } /* end if */
      ind[index][j] = k;              /* k-th grain is near pix. "index". */
      loc[index][j] = i;              /* Pix. "index" in i-th loc. */
      cval1[index][j] = convolution1;
      cval2[index][j] = convolution2;
      nopresent[index] = j+1;
    } /* end for i. */
  } /* end for k. */

  /* Copy data to MATLAB structure: */
  ndim = 2;
  vecdim[0] = dims;
  vecdim[1] = 4;
  plhs[0] = mxCreateCellArray(ndim,vecdim);
  for (index=0;index<dims;index++){
    m = nopresent[index];
    cell1 = mxCreateDoubleMatrix(1,m,mxREAL);
    cell2 = mxCreateDoubleMatrix(1,m,mxREAL);
    cell3 = mxCreateDoubleMatrix(1,m,mxREAL);
    cell4 = mxCreateDoubleMatrix(1,m,mxREAL);
    p1 = mxGetPr(cell1);
    p2 = mxGetPr(cell2);
    p3 = mxGetPr(cell3);
    p4 = mxGetPr(cell4);
    
    /* stopped here */
    
    /* Sort the data w.r.t. cval: */
    set_INC_and_CVAL(cval1[index],m);
    set_INC_and_CVAL(cval2[index],m);
    qsort(INC,m,sizeof(int),(int (*)(const void *, const void *)) comparison);
    /* Place data in MATLAB data structure: */
    for (i=0;i<m;i++){
      p1[i] = ind[index][INC[i]]+1;
      p2[i] = cval1[index][INC[i]];
      p3[i] = cval2[index][INC[i]];
      p4[i] = loc[index][INC[i]]+1;
    }
    mxSetCell(plhs[0],index,cell1);
    mxSetCell(plhs[0],index+dims,cell2);
    mxSetCell(plhs[0],index+2*dims,cell3);
    mxSetCell(plhs[0],index+3*dims,cell4);
  }      

  /* Free dynamically allocated memory: */    
  for (i=0;i<dims;i++){
    free(ind[i]);
    free(loc[i]);
    free(cval1[i]);
    free(cval2[i]);
  }
  free(ind);
  free(loc);
  free(cval1);
  free(cval2);
  free(nopresent);
  free(limit);
  
}

int comparison(int *a, int *b)
{
  /* Needed for sorting via qsort. */
  if (CVAL[*a] < CVAL[*b]) return 1;
  else if (CVAL[*a] == CVAL[*b]) return 0;
  return -1;
}

void set_INC_and_CVAL(double *cval, int m)
{
  /* Needed for sorting via qsort. */
  int i;
  
  for (i=0;i<m;i++){
    INC[i] = i;
    CVAL[i] = cval[i];
  }
}