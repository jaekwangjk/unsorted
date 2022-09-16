#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

int find_minimum(int[], int);
int find_maximum(int[], int);


struct grid2dsubset {
  int *x; /* List of element indices. */
  int *y; /* Ditto. */
  int *bx;/* List of bdry element indices. */
  int *by;/* Ditto. */
  int *p; /* Stores position of pixels in element list. */
  int m; /* Inner dimension of grid. */
  int n; /* Outer dimension of grid. */
  int N; /* Number of elements. */
  int bN; /* Number of boundary elements. */
  int maxN; /* Max number of elements allowed. */
};

struct intarray1d {
  int *D; /* Points to beginning of data. */
  int N;  /* Number of elements. */
};

void create_grid2dsubset(struct grid2dsubset *in, int m, int n, int Maxn, int *W, int *Wx, int *Wy, int *Wbx, int *Wby);
void destroy_grid2dsubset(struct grid2dsubset *in);
void add_element_grid2dsubset(struct grid2dsubset *in, int x, int y);
int isbdry(struct grid2dsubset *in, int x, int y);
void updatebdry_grid2dsubset(struct grid2dsubset *in);
void grow_grid2dsubset(struct grid2dsubset *in, int w, int *Wcx, int *Wcy);
void create_intarray1d(struct intarray1d *in, int N, int *W);
void destroy_intarray1d(struct intarray1d *u);

/*************************************************************************/
void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{

  /* **********************************************************************
   * MATLAB Interface:
   * **********************************************************************
   * [xo,yo] = pgrow(x,y,w,W);
   * **********************************************************************
   * Grows (dilates) the subset of grid by w pixels outwards.
   * Periodic boundary conditions used.
   * The vectors x and y should be COLUMN VECTORS of type int32!
   * **********************************************************************
   */
  
  int *x, *y, w, *W;
  int *Wx, *Wy, *Wbx, *Wby, *Wcx, *Wcy;
  double *xout, *yout;
  int i, j, ind, numpix, maxN, m, n;
  struct grid2dsubset S;
  
  /* Input variables: */
  x = (int *) mxGetData(prhs[0]); /* Array of x coords. */
  y = (int *) mxGetData(prhs[1]); /* Array of y coords. */
  
  // w= R = Minimum distance to maintain between grains in the same family.
  w = (int) mxGetScalar(prhs[2]); /* Width of growth. */
  W = (int *) mxGetData(prhs[3]); /* Workspace. */
  Wx = (int *) mxGetData(prhs[4]); /* Workspace for (*in).x */
  Wy = (int *) mxGetData(prhs[5]); /* Workspace for (*in).y */
  Wbx = (int *) mxGetData(prhs[6]); /* Workspace for (*in).bx */
  Wby = (int *) mxGetData(prhs[7]); /* Workspace for (*in).by */
  Wcx = (int *) mxGetData(prhs[8]); /* Workspace for (*in).cx */
  Wcy = (int *) mxGetData(prhs[9]); /* Workspace for (*in).cy */
  
  numpix = (int) mxGetM(prhs[0]);
  m = (int) mxGetM(prhs[3]); /* Dimension of grid. */
  n = (int) mxGetN(prhs[3]); /* Dimension of grid. */
  maxN = mxGetM(prhs[4]);  /* Number of pixels currently in the subset. */
  
  create_grid2dsubset(&S, m, n, maxN, W, Wx, Wy, Wbx, Wby);
  
  for (j=0;j<numpix;j++){ /* Fill the subset with points from list of coords. */
    add_element_grid2dsubset(&S,x[j]-1,y[j]-1);
  }
  
  /*printf("numpix = %d, %d\n",numpix,2*(maxx-minx+1+2*w)*(maxy-miny+1+2*w));*/
  updatebdry_grid2dsubset(&S); /* Boundary data structure made consistent. */
  
  grow_grid2dsubset(&S,w,Wcx,Wcy); /* Real work happens here. */
  
  /* Allocate memory and get pointer for output variables. */
  plhs[0] = mxCreateDoubleMatrix(S.N,1,mxREAL);	
  xout = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(S.N,1,mxREAL);
  yout = mxGetPr(plhs[1]);
  
  for (i=0;i<S.N;i++){ /* Fill output variables with coordinates. */
    *(xout+i) = S.x[i]+1;
    *(yout+i) = S.y[i]+1;
  }

  destroy_grid2dsubset(&S);
}


/*************************************************************************/
void create_grid2dsubset(struct grid2dsubset *in, int m, int n, int maxN, int *W, int *Wx, int *Wy, int *Wbx, int *Wby)
{
  int i,j,ind;
  
  /* Allocate memory. Set to empty subset. */
  (*in).x = Wx; /* Workspace array that is provided by caller is used. */
  (*in).y = Wy; /* Workspace array that is provided by caller is used. */
  (*in).bx = Wbx; /* Workspace array that is provided by caller is used. */
  (*in).by = Wby; /* Workspace array that is provided by caller is used. */
  (*in).m = m;
  (*in).n = n;
  (*in).p = W; /* Workspace array that is provided by caller is used. */
  (*in).N = 0;
  (*in).bN = 0;
  (*in).maxN = maxN;
      
}
/*************************************************************************/
void destroy_grid2dsubset(struct grid2dsubset *in)
{
  int i, n, ind, x, y;
      
  /* Clear the workspace array: */
  n = (*in).n;
  for (i=0;i<(*in).N;i++) {
    x = (*in).x[i];
    y = (*in).y[i];
    ind = x*n + y;
    (*in).p[ind] = -1;
    /*(*in).x[i] = 0;
    (*in).y[i] = 0;
    (*in).bx[i] = 0;
    (*in).by[i] = 0;*/
  }
}
/*************************************************************************/
void add_element_grid2dsubset(struct grid2dsubset *in, int x, int y)
{
  /* Adds element to grid subset. Note that boundary is NOT updated. */
  int m,n,N,ind;  
  N = (*in).N;
  
  m = (*in).m;
  n = (*in).n;
  
  if (N == 0) {/* If subset currently empty... */
    (*in).x[0] = x;
    (*in).y[0] = y;
    ind = x*n+y;
    (*in).p[ind] = 0;
    ((*in).N)++;
  }
  else { /* If subset already has elements... */
    ind = x*n+y;
    if ( (*in).p[ind] == -1 ) { /* If new loc. not already in subset... */
      (*in).p[ind] = N; /* Stores position of new element in the list. */
      (*in).x[N] = x; /* Add new x coord to list of x coords. */
      (*in).y[N] = y; /* Add new y coord to list of y coords. */
      ((*in).N)++;
    } /* end if */
  } /* end else. */
}
/*************************************************************************/
int isbdry(struct grid2dsubset *in, int x, int y)
{
  int neighmin, flag, N, S, E, W, m, n, ind;

  m = (*in).m;
  n = (*in).n;

  E = x+1;
  if (E == m) E = 0;
  W = x-1;
  if (W == -1) W = m-1;
  N = y+1;
  if (N == n) N = 0;
  S = y-1;
  if (S == -1) S = n-1;
  
  neighmin = 0;
  ind = E*n+y;
  neighmin = min( neighmin, (*in).p[ind] );
  ind = W*n+y;
  neighmin = min( neighmin, (*in).p[ind] );
  ind = x*n + N;
  neighmin = min( neighmin, (*in).p[ind] );
  ind = x*n + S;
  neighmin = min( neighmin, (*in).p[ind] );
    
  if (neighmin < 0) flag = 1;
  else flag = 0;

  return flag;
}
/*************************************************************************/
void updatebdry_grid2dsubset(struct grid2dsubset *in)
{
  int i,x,y,N,bN;
  
  N = (*in).N; /* Number of elements in grid subset. */
  bN = 0;  

  for (i=0;i<N;i++){
    x = (*in).x[i];
    y = (*in).y[i];
    if ( isbdry(in,x,y) == 1 ) {
      (*in).bx[bN] = x;
      (*in).by[bN] = y;
      bN++;
    }
  }
  
  (*in).bN = bN;
}
/*************************************************************************/
void grow_grid2dsubset(struct grid2dsubset *in, int w, int *Wcx, int *Wcy)
{
  int i,j,m,n,t,x,y,bx,by,pos,count,bN,ind,E,W,N,S,maxN;
  struct intarray1d candidate_x, candidate_y;

  m = (*in).m;
  n = (*in).n;
  maxN = (*in).maxN;
  
  create_intarray1d(&candidate_x,maxN,Wcx);
  create_intarray1d(&candidate_y,maxN,Wcy);
  
  
  
  for (t=0;t<w;t++){

    count = 0; /* Counts how many new points added to the subset
		              on current pass. */
    for (i=0;i<(*in).bN;i++){ /* Go through current boundary points. */
      bx = (*in).bx[i];
      by = (*in).by[i];

      /* Add neighbors of point (bx,by) if they aren't already in subset.
       * Four-neighborhood is used.
       */

      /* East neighbor: */
      E = bx+1;
      if (E==m) E=0;
      ind = E*n+by;
      if ( (*in).p[ind] < 0 ) {
        pos = (*in).N;
        (*in).x[pos] = E;
        (*in).y[pos] = by;
        (*in).N++;
        (*in).p[ind] = pos;
        candidate_x.D[count] = E;
        candidate_y.D[count] = by;
        count++;
      }
      
      /* West neighbor: */
      W = bx-1;
      if (W==-1) W=m-1;
      ind = W*n+by;
      if ( (*in).p[ind] < 0 ) {
        pos = (*in).N;
        (*in).x[pos] = W;
        (*in).y[pos] = by;
        (*in).N++;
        (*in).p[ind] = pos;
        candidate_x.D[count] = W;
        candidate_y.D[count] = by;
        count++;
      }
      
      /* North neighbor: */
      N = by+1;
      if (N==n) N=0;
      ind = bx*n+N;
      if ( (*in).p[ind] < 0 ) {
        pos = (*in).N;
        (*in).x[pos] = bx;
        (*in).y[pos] = N;
        (*in).N++;
        (*in).p[ind] = pos;
        candidate_x.D[count] = bx;
        candidate_y.D[count] = N;
        count++;
      }
      
      /* South neighbor: */
      S = by-1;
      if (S==-1) S=n-1;
      ind = bx*n+S;
      if ( (*in).p[ind] < 0 ) {
        pos = (*in).N;
        (*in).x[pos] = bx;
        (*in).y[pos] = S;
        (*in).N++;
        (*in).p[ind] = pos;
        candidate_x.D[count] = bx;
        candidate_y.D[count] = S;
        count++;
      }
      
    } /* i loop (through current bdry points) ends. */
    
    /* Check the newly added points; new bdry is a subset of those. */
    bN = 0;
    for (i=0;i<count;i++){
      x = candidate_x.D[i];
      y = candidate_y.D[i];
      if ( isbdry(in,x,y) == 1 ) {
        (*in).bx[bN] = x;
        (*in).by[bN] = y;
        bN++;
      }
    }
    (*in).bN = bN;
        
    
  } /* t loop ends. */
  /*destroy_intarray1d(&candidate_x);
  destroy_intarray1d(&candidate_y);*/
  
    
}
/*************************************************************************/
void create_intarray1d(struct intarray1d *in, int N, int *W)
{
  (*in).D = W;
  (*in).N = N;
}
/*************************************************************************/
void destroy_intarray1d(struct intarray1d *in)
{
  int i;
  /* Clear the workspace array: */
  for (i=0;i<(*in).N;i++) {
    (*in).D = 0;
  }
}
/*************************************************************************/