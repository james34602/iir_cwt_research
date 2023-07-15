#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <float.h>
#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif
#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
typedef struct
{
	double real, imag;
} cplx;
void procNpole1ZeroCplx(cplx *x, double *X_inst_freqs, cplx *y, unsigned int len, char interpolation)
{
	unsigned int i;
	if (interpolation)
	{
		for (i = 0; i < len; i++)
		{
			unsigned int range1 = round(X_inst_freqs[i]);
			if (range1 >= (len - 1))
			{
				y[len - 1].real = y[len - 1].real + x[i].real;
				y[len - 1].imag = y[len - 1].imag + x[i].imag;
				continue;
			}
			double beta = X_inst_freqs[i] - (double)range1;
			if (beta > 0.0)
			{
				y[range1 + 1].real = y[range1 + 1].real + x[i].real * beta;
				y[range1 + 1].imag = y[range1 + 1].imag + x[i].imag * beta;
				y[range1].real = y[range1].real + x[i].real * (1 - beta);
				y[range1].imag = y[range1].imag + x[i].imag * (1 - beta);
			}
			else if (beta < 0.0)
			{
				y[range1 - 1].real = y[range1 - 1].real + x[i].real * -beta;
				y[range1 - 1].imag = y[range1 - 1].imag + x[i].imag * -beta;
				y[range1].real = y[range1].real + x[i].real * (1 + beta);
				y[range1].imag = y[range1].imag + x[i].imag * (1 + beta);
			}
			else
			{
				y[range1].real = y[range1].real + x[i].real;
				y[range1].imag = y[range1].imag + x[i].imag;
			}
		}
	}
	else
	{
		for (i = 0; i < len; i++)
		{
			unsigned int range1 = round(X_inst_freqs[i]);
			if (range1 >= (len - 1))
			{
				y[len - 1].real = y[len - 1].real + x[i].real;
				y[len - 1].imag = y[len - 1].imag + x[i].imag;
				continue;
			}
			y[range1].real = y[range1].real + x[i].real;
			y[range1].imag = y[range1].imag + x[i].imag;
		}
	}
}
void procNpole1ZeroReal(double *x, double *X_inst_freqs, double *y, unsigned int len, char interpolation)
{
	unsigned int i;
	if (interpolation)
	{
		for (i = 0; i < len; i++)
		{
			unsigned int range1 = round(X_inst_freqs[i]);
			if (range1 >= (len - 1))
			{
				y[len - 1] = y[len - 1] + x[i];
				continue;
			}
			double beta = X_inst_freqs[i] - (double)range1;
			if (beta > 0.0)
			{
				y[range1 + 1] = y[range1 + 1] + x[i] * beta;
				y[range1] = y[range1] + x[i] * (1 - beta);
			}
			else if (beta < 0.0)
			{
				y[range1 - 1] = y[range1 - 1] + x[i] * -beta;
				y[range1] = y[range1] + x[i] * (1 + beta);
			}
			else
				y[range1] = y[range1] + x[i];
		}
	}
	else
	{
		for (i = 0; i < len; i++)
		{
			unsigned int range1 = round(X_inst_freqs[i]);
			if (range1 >= (len - 1))
			{
				y[len - 1] = y[len - 1] + x[i];
				continue;
			}
			y[range1] = y[range1] + x[i];
		}
	}
}
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	size_t i, j, m, n;
	// Check for proper number of arguments
#ifndef NOCHECK
	if (nrhs != 3)
		mexErrMsgIdAndTxt("MATLAB:cplxReassignment:invalidNumInputs", "3 input arguments are required.");
	else if (nlhs > 1)
		mexErrMsgIdAndTxt("MATLAB:cplxReassignment:maxlhs", "Too many output arguments.");
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	if (!mxIsDouble(prhs[0]) || (MIN(m, n) != 1))
		mexErrMsgIdAndTxt("MATLAB:cplxReassignment:incorrectType", "First argument is X");
#endif
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
#ifndef NOCHECK
	if (!mxIsDouble(prhs[0]) || (MIN(m, n) != 1))
		mexErrMsgIdAndTxt("MATLAB:cplxReassignment:incorrectType", "Second argument is X_inst_freqs");
#endif
	unsigned int len = MAX(m, n);
	double *X_inst_freqs = mxGetPr(prhs[1]);
	double *outputMode = mxGetPr(prhs[2]);
	if (!mxIsComplex(prhs[0]))
	{
		plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
		double *X = mxGetPr(prhs[0]);
		double *y = mxGetPr(plhs[0]);
		procNpole1ZeroReal(X, X_inst_freqs, y, len, *outputMode == 1);
	}
	else
	{
		plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
		mxComplexDouble *X = mxGetComplexDoubles(prhs[0]);
		mxComplexDouble *y = mxGetComplexDoubles(plhs[0]);
		procNpole1ZeroCplx(X, X_inst_freqs, y, len, *outputMode == 1);
	}
	return;
}