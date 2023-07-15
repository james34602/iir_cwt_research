#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "Npole1zero.h"
#include <mex.h>
#include <float.h>
#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif
#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    size_t i, j, m, n;
    // Check for proper number of arguments
    if (nrhs != 4)
        mexErrMsgIdAndTxt("MATLAB:procNPole1Zero:invalidNumInputs", "2 input arguments are required.");
    else if (nlhs > 1)
        mexErrMsgIdAndTxt("MATLAB:procNPole1Zero:maxlhs", "Too many output arguments.");
    
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (!mxIsUint8(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || (MIN(m, n) != 1))
        mexErrMsgIdAndTxt("MATLAB:procNPole1Zero:incorrectType", "Requires a unsigned char vector as state variable.");
    char *yp = (char*)mxGetPr(prhs[0]);
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    if (!mxIsDouble(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[3]) || mxIsSparse(prhs[3]) || (MIN(m, n) != 1))
        mexErrMsgIdAndTxt("MATLAB:procNPole1Zero:incorrectType", "Array must be double vector");
    m = mxGetM(prhs[3]);
    n = mxGetN(prhs[3]);
    unsigned int numBands = max(m, n);
    double *b = mxGetPr(prhs[1]);
    double *a = mxGetPr(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    mxComplexDouble *y = mxGetComplexDoubles(plhs[0]);
    if(!mxIsComplex(prhs[3]))
    {
        double *x = mxGetPr(prhs[3]);
        procNpole1ZeroReal(yp, b, a, x, y);
    }
    else
    {
        mxComplexDouble *x = mxGetComplexDoubles(prhs[3]);
        procNpole1ZeroCplx(yp, b, a, x, y);
    }
    //y[0].imag += DBL_TRUE_MIN; // Fix MATLAB BUG that zero in zero out go slow
    return;
}