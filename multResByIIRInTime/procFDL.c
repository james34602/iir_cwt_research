#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "fdl.h"
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
    if (nrhs != 2)
        mexErrMsgIdAndTxt("MATLAB:procFDL:invalidNumInputs", "2 input arguments are required.");
    else if (nlhs > 1)
        mexErrMsgIdAndTxt("MATLAB:procFDL:maxlhs", "Too many output arguments.");
    
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (!mxIsUint8(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || (MIN(m, n) != 1))
        mexErrMsgIdAndTxt("MATLAB:procFDL:incorrectType", "Requires a unsigned char vector as state variable.");
    fractionalDelayLine *yp = (fractionalDelayLine*)mxGetPr(prhs[0]);
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    if (!mxIsDouble(prhs[1]) || mxIsSparse(prhs[1]) || (MIN(m, n) != 1))
        mexErrMsgIdAndTxt("MATLAB:procFDL:incorrectType", "Array must be double vector");
    unsigned int numBands = max(m, n);
    plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    mxComplexDouble *y = mxGetComplexDoubles(plhs[0]);
    if(!mxIsComplex(prhs[1]))
    {
        double *x = mxGetPr(prhs[1]);
        fractionalDelayLineProcessReal(yp, x, y);
    }
    else
    {
        mxComplexDouble *x = mxGetComplexDoubles(prhs[1]);
        fractionalDelayLineProcessCplx(yp, x, y);
    }
    y[0].imag += DBL_TRUE_MIN; // Fix MATLAB BUG that zero in zero out go slow
    return;
}