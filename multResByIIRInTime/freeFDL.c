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
    if (nrhs != 1)
        mexErrMsgIdAndTxt("MATLAB:freeFDL:invalidNumInputs", "No input arguments required.");
    else if (nlhs > 0)
        mexErrMsgIdAndTxt("MATLAB:freeFDL:maxlhs", "No output argument is required");
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (!mxIsUint8(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || (MIN(m, n) != 1))
        mexErrMsgIdAndTxt("MATLAB:procCplxMovAvg:incorrectType", "Requires a unsigned char vector as state variable.");
    fractionalDelayLine *yp = (fractionalDelayLine*)mxGetPr(prhs[0]);
	fractionalDelayLine_free(yp);
    return;
}