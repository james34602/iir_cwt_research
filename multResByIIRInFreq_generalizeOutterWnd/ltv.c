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
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	size_t i, j, m, n;
	// Check for proper number of arguments
#ifndef NOCHECK
	if (nrhs != 6)
		mexErrMsgIdAndTxt("MATLAB:ltv:invalidNumInputs", "6 input arguments are required.");
	else if (nlhs > 1)
		mexErrMsgIdAndTxt("MATLAB:ltv:maxlhs", "Too many output arguments.");
#endif
	// First input is complex spectrum
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	unsigned int len = (unsigned int)MAX(m, n);
#ifndef NOCHECK
	if (!mxIsDouble(prhs[0]) || (MIN(m, n) != 2))
		mexErrMsgIdAndTxt("MATLAB:ltv:incorrectType", "First argument is X");
	if (!mxIsDouble(prhs[0]) || (MIN(m, n) != 2))
		mexErrMsgIdAndTxt("MATLAB:ltv:incorrectType", "Second argument is tmp");
#endif
	// Other input are IIR coefficients
	double *b = mxGetPr(prhs[2]);
	double *a = mxGetPr(prhs[3]);
	double *c1 = mxGetPr(prhs[4]);
    double *c2 = mxGetPr(prhs[5]);
    mxComplexDouble o[2], st;
    mxComplexDouble z1, z1_;
    mxComplexDouble z2, z2_;
	if (mxIsEmpty(prhs[2]))
	{
		// First order
		if (mxIsComplex(prhs[0]))
		{
			plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
			mxComplexDouble *y = mxGetComplexDoubles(plhs[0]);
			mxComplexDouble *y2 = y + m;
			mxComplexDouble *X = mxGetComplexDoubles(prhs[0]);
			mxComplexDouble *X2 = X + m;
			mxComplexDouble *tmp = mxGetComplexDoubles(prhs[1]);
			mxComplexDouble *tmp2 = tmp + m;
            z1.real = z1_.real = z1.imag = z1_.imag = 0;
			for (i = 0; i < m; i++)
			{
				z1.real = z1.real + a[i] * (X[i].real - z1.real);
				z1.imag = z1.imag + a[i] * (X[i].imag - z1.imag);

				z1_.real = z1_.real + a[i] * (X2[i].real - z1_.real);
				z1_.imag = z1_.imag + a[i] * (X2[i].imag - z1_.imag);

				tmp[i].real = z1.real;
				tmp[i].imag = z1.imag;
				tmp2[i].real = z1_.real;
				tmp2[i].imag = z1_.imag;
			}
			z1.real = z1_.real = z1.imag = z1_.imag = 0;
			for (i = m; i-- > 0; )
			{
				z1.real = z1.real + a[i] * (tmp[i].real - z1.real);
				z1.imag = z1.imag + a[i] * (tmp[i].imag - z1.imag);

				z1_.real = z1_.real + a[i] * (tmp2[i].real - z1_.real);
				z1_.imag = z1_.imag + a[i] * (tmp2[i].imag - z1_.imag);

				y[i].real = z1.real;
				y[i].imag = z1.imag;
				y2[i].real = z1_.real;
				y2[i].imag = z1_.imag;
			}
		}
		else
		{
			plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
			double *y = mxGetPr(plhs[0]);
			double *y2 = y + m;
			double *X = mxGetPr(prhs[0]);
			double *X2 = X + m;
			double *tmp = mxGetPr(prhs[1]);
			double *tmp2 = tmp + m;
			z1.real = z1_.real = 0;
			for (i = 0; i < m; i++)
			{
				z1.real = z1.real + a[i] * (X[i] - z1.real);
				z1_.real = z1_.real + a[i] * (X2[i] - z1_.real);
				tmp[i] = z1.real;
				tmp2[i] = z1_.real;
			}
			z1.real = z1_.real = 0;
			for (i = m; i-- > 0; )
			{
				z1.real = z1.real + a[i] * (tmp[i] - z1.real);
				z1_.real = z1_.real + a[i] * (tmp2[i] - z1_.real);
				y[i] = z1.real;
				y2[i] = z1_.real;
			}
		}
	}
	else
	{
		// Second order
		if (mxIsComplex(prhs[0]))
		{
			plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
			mxComplexDouble *y = mxGetComplexDoubles(plhs[0]);
			mxComplexDouble *y2 = y + m;
			mxComplexDouble *X = mxGetComplexDoubles(prhs[0]);
			mxComplexDouble *X2 = X + m;
			mxComplexDouble *tmp = mxGetComplexDoubles(prhs[1]);
			mxComplexDouble *tmp2 = tmp + m;
			z1.real = z1_.real = z1.imag = z1_.imag = 0;
			z2.real = z2_.real = z2.imag = z2_.imag = 0;
			for (i = 0; i < m; i++)
			{
				o[0].real = X[i].real - z1.real - z2.real;
				o[0].imag = X[i].imag - z1.imag - z2.imag;
				st.real = c2[i] * z1.real;
				st.imag = c2[i] * z1.imag;
				tmp[i].real = b[i] * o[0].real + 2 * st.real + z2.real;
				tmp[i].imag = b[i] * o[0].imag + 2 * st.imag + z2.imag;
				z2.real = z2.real + st.real;
				z2.imag = z2.imag + st.imag;
				z1.real = z1.real + c1[i] * o[0].real;
				z1.imag = z1.imag + c1[i] * o[0].imag;

				o[1].real = X2[i].real - z1_.real - z2_.real;
				o[1].imag = X2[i].imag - z1_.imag - z2_.imag;
				st.real = c2[i] * z1_.real;
				st.imag = c2[i] * z1_.imag;
				tmp2[i].real = b[i] * o[1].real + 2 * st.real + z2_.real;
				tmp2[i].imag = b[i] * o[1].imag + 2 * st.imag + z2_.imag;
				z2_.real = z2_.real + st.real;
				z2_.imag = z2_.imag + st.imag;
				z1_.real = z1_.real + c1[i] * o[1].real;
				z1_.imag = z1_.imag + c1[i] * o[1].imag;
			}
			z1.real = z1_.real = z1.imag = z1_.imag = 0;
			z2.real = z2_.real = z2.imag = z2_.imag = 0;
			for (i = m; i-- > 0; )
			{
				o[0].real = tmp[i].real - z1.real - z2.real;
				o[0].imag = tmp[i].imag - z1.imag - z2.imag;
				st.real = c2[i] * z1.real;
				st.imag = c2[i] * z1.imag;
				y[i].real = b[i] * o[0].real + 2 * st.real + z2.real;
				y[i].imag = b[i] * o[0].imag + 2 * st.imag + z2.imag;
				z2.real = z2.real + st.real;
				z2.imag = z2.imag + st.imag;
				z1.real = z1.real + c1[i] * o[0].real;
				z1.imag = z1.imag + c1[i] * o[0].imag;

				o[1].real = tmp2[i].real - z1_.real - z2_.real;
				o[1].imag = tmp2[i].imag - z1_.imag - z2_.imag;
				st.real = c2[i] * z1_.real;
				st.imag = c2[i] * z1_.imag;
				y2[i].real = b[i] * o[1].real + 2 * st.real + z2_.real;
				y2[i].imag = b[i] * o[1].imag + 2 * st.imag + z2_.imag;
				z2_.real = z2_.real + st.real;
				z2_.imag = z2_.imag + st.imag;
				z1_.real = z1_.real + c1[i] * o[1].real;
				z1_.imag = z1_.imag + c1[i] * o[1].imag;
			}
		}
		else
		{
			plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
			double *y = mxGetPr(plhs[0]);
			double *y2 = y + m;
			double *X = mxGetPr(prhs[0]);
			double *X2 = X + m;
			double *tmp = mxGetPr(prhs[1]);
			double *tmp2 = tmp + m;
			z1.real = z1_.real = 0;
			z2.real = z2_.real = 0;
			for (i = 0; i < m; i++)
			{
				o[0].real = X[i] - z1.real - z2.real;
				st.real = c2[i] * z1.real;
				tmp[i] = b[i] * o[0].real + 2 * st.real + z2.real;
				z2.real = z2.real + st.real;
				z1.real = z1.real + c1[i] * o[0].real;

				o[1].real = X2[i] - z1_.real - z2_.real;
				st.real = c2[i] * z1_.real;
				tmp2[i] = b[i] * o[1].real + 2 * st.real + z2_.real;
				z2_.real = z2_.real + st.real;
				z1_.real = z1_.real + c1[i] * o[1].real;
			}
			z1.real = z1_.real = 0;
			z2.real = z2_.real = 0;
			for (i = m; i-- > 0; )
			{
				o[0].real = tmp[i] - z1.real - z2.real;
				st.real = c2[i] * z1.real;
				y[i] = b[i] * o[0].real + 2 * st.real + z2.real;
				z2.real = z2.real + st.real;
				z1.real = z1.real + c1[i] * o[0].real;

				o[1].real = tmp2[i] - z1_.real - z2_.real;
				st.real = c2[i] * z1_.real;
				y2[i] = b[i] * o[1].real + 2 * st.real + z2_.real;
				z2_.real = z2_.real + st.real;
				z1_.real = z1_.real + c1[i] * o[1].real;
			}
		}
	}
	return;
}