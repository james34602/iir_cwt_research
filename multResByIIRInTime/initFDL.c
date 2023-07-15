#include <stdio.h>
#include <math.h>
#include <float.h>
#include "fdl.h"
#include <mex.h>
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
        mexErrMsgIdAndTxt("MATLAB:initFDL:invalidNumInputs", "2 input arguments required.");
    else if (nlhs > 1)
        mexErrMsgIdAndTxt("MATLAB:initFDL:maxlhs", "Too many output arguments.");
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    size_t numElements = MAX(m, n);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || (MIN(m, n) < 1))
        mexErrMsgIdAndTxt("MATLAB:initCplxMovAvg:incorrectType", "Specific a array of delay in double");
    double *lag = (double*)mxGetPr(prhs[0]);
    // Create a matrix for the return argument
    plhs[0] = mxCreateNumericMatrix(1, sizeof(fractionalDelayLine), mxUINT8_CLASS, mxREAL);
    // Assign pointers to the various parameters
	char *yp = (char*)mxGetPr(plhs[0]);
    // Do the actual computations in a subroutine
	fractionalDelayLineInit(yp, numElements, lag, *((double*)mxGetPr(prhs[1])));
	fractionalDelayLine_setDelay(yp, lag);
    return;
}