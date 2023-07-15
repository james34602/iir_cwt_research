#include <stdio.h>
#include <math.h>
#include <float.h>
#include "Npole1zero.h"
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
        mexErrMsgIdAndTxt("MATLAB:initNPole1Zero:invalidNumInputs", "Please specific nBins and orders");
    else if (nlhs > 1)
        mexErrMsgIdAndTxt("MATLAB:initNPole1Zero:maxlhs", "Too many output arguments.");
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (!mxIsUint32(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || m > 1 || n > 1)
        mexErrMsgIdAndTxt("MATLAB:initCplxMovAvg:incorrectType", "Please specific nBins");
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    if (!mxIsUint32(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || m > 1 || n > 1)
        mexErrMsgIdAndTxt("MATLAB:initCplxMovAvg:incorrectType", "Please specific orders");
    unsigned int *nBands = (unsigned int*)mxGetPr(prhs[0]);
    unsigned int *orders = (unsigned int*)mxGetPr(prhs[1]);
    size_t dataStructSize = estNpole1Zero(*nBands, *orders);
    // Create a matrix for the return argument
    plhs[0] = mxCreateNumericMatrix(1, dataStructSize, mxUINT8_CLASS, mxREAL);
    // Assign pointers to the various parameters
	char *yp = (char*)mxGetPr(plhs[0]);
    // Do the actual computations in a subroutine
	initNpole1Zero(yp, *nBands, *orders);
    return;
}