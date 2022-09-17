// Simple example of a MEX function. Save this as mexfunc.cpp and compile with
// mex; it can then be called as mexfunc(something) from MATLAB.

#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        mxArray *a_in, *b_out;
        const mwSize *dims;
        double *a, *b; 
        int rows, cols;

        // Input array:
        a_in = mxDuplicateArray(prhs[0]);

        // Get dimensions of input array:
        dims = mxGetDimensions(prhs[0]);
        rows = (int) dims[0];
        cols = (int) dims[1];

        // Output array:
        b_out = plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
        
        // Access the contents of the input and output arrays:
        a = mxGetPr(a_in);
        b = mxGetPr(b_out);

        // Compute some elementwise function of the input array:
        for (int i = 0; i < rows*cols; i++) {
                b[i] = 3*sin(a[i]) + 2*cos(a[i]);
    }
}
