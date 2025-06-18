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
#define M_PI	3.14159265358979323846
double hz2mel(double hz)
{
	return 2595 * log10(1.0 + (hz / 700));
}
double mel2hz(double mel)
{
	return 700 * (pow(10, (mel / 2595.0)) - 1.0);
}
#include <float.h>
double wrightOmegaq(double X)
{
	double EXP1 = exp(1);
	double EXP2 = exp(2);
	double LN2 = log(2);
	double OMEGA = 0.5671432904097838;
	double ONE_THIRD = 1.0 / 3.0;
	double W;
	// Special values
	if (X > pow(2.0, 59.0))
		W = X; // W self-saturates: X > 2^59 (abs(Y) > 2^54 too)
	else if (X == 0.0)
		W = OMEGA; // Omega constant
	else if (X == 1.0)
		W = 1;
	else if (X == 1.0 + EXP1)
		W = EXP1;
	else
	{
		if (X < log(DBL_EPSILON * DBL_MIN) - LN2)
			W = 0.0; // Z -> -Inf
		else
		{
			// W used in order retain datatype
			if (X <= -2.0)
			{
				// Region 3: series about -Inf
				double x = exp(X);
				W = x * (1.0 - x * (1.0 - x * (36.0 - x * (64.0 - 125.0 * x)) / 24.0));
				// Series is exact, X < -exp(2)
				if (X < -EXP2)
					return W;
			}
			else if (X > M_PI + 1)
			{
				// Region 7: log series about Z = Inf
				double x = log(X);
				double lzi = x / X;
				W = X - x + lzi * (1.0 + lzi * (0.5 * x - 1.0 + lzi * ((ONE_THIRD * x - 1.5) * x + 1)));
			}
			else
			{
				// Region 4: series about Z = 1
				double x = X - 1.0;
				W = 1.0 + x * (1.0 / 2.0 + x * (1.0 / 16.0 - x * (1.0 / 192.0 + x * (1.0 / 3072.0 - (13.0 / 61440.0) * x))));
			}
			// Residual
			double r = X - (W + log(W));
			if (fabs(r) > DBL_EPSILON)
			{
				// FSC-type iteration, N = 3, (Fritsch, Shafer, & Crowley, 1973)
				double w1 = 1.0 + W;
				double w2 = w1 + 2 * ONE_THIRD * r;
				W = W * (1.0 + r * (w1 * w2 - 0.5 * r) / (w1 * (w1 * w2 - r)));
				// Test residual
				r = X - (W + log(W));
				// Second iterative improvement via FSC method, if needed
				if (fabs(r) > DBL_EPSILON)
				{
					w1 = 1.0 + W;
					w2 = w1 + 2 * ONE_THIRD * r;
					W = W * (1.0 + r * (w1 * w2 - 0.5 * r) / (w1 * (w1 * w2 - r)));
				}
			}
		}
	}
	return W;
}
double linLgTransform(double x, double a)
{
	return x * a + log(x) * (1.0 - a);
}
double linExpTransform(double x, double a)
{
	return -(wrightOmegaq(-log(-(a - 1.0) / a) - x / (a - 1)) * (a - 1.0)) / a;
}
double lgTransform(double x, double linLogRatio)
{
	//return log(x);
	//return hz2mel(x);
	return linLgTransform(x, linLogRatio);
}
double expTransform(double x, double linLogRatio)
{
	//return exp(x);
	//return mel2hz(x);
	return linExpTransform(x, linLogRatio);
}
typedef struct
{
	double real, imag;
} cplx;
static double gridMap(int halfLen, int magindex, double min_freq, double max_freq, int samplerate, double linLogRatio)
{
	// Log scale
	double low2 = lgTransform(min_freq, linLogRatio);
	double high2 = lgTransform(max_freq, linLogRatio);
	double lgVect = low2 + (high2 - low2) * magindex / (halfLen - 1);
	double freq = expTransform(lgVect, linLogRatio);
	return (freq * halfLen / (samplerate / 2));
}
static void gridTransformReal(double *output, double *input, double *map, int halfLen, double lowerBoundHz, double upperBoundHz, double samplerate, double linLogRatio)
{
	for (int k = 0; k < halfLen; k++)
	{	/* Average the pixels in the range it comes from */
		double this = gridMap(halfLen, k, lowerBoundHz, upperBoundHz, samplerate, linLogRatio);
		double next = gridMap(halfLen, k + 1, lowerBoundHz, upperBoundHz, samplerate, linLogRatio);
		map[k] = (this + next) * 0.5;

		/* Range check: can happen if --max-freq > samplerate / 2 */
		if (this > halfLen)
		{
			output[k] = 0.0;
			return;
		}

		if (next > this + 1)
		{	/* The output indices are more sparse than the input indices
			** so average the range of input indices that map to this output,
			** making sure not to exceed the input array (0..halfLen inclusive)
			*/
			/* Take a proportional part of the first sample */
			double count = 1.0 - (this - floor(this));
			double sumRe = input[(int)this] * count;

			while ((this += 1.0) < next && (int)this <= halfLen)
			{
				sumRe += input[(int)this];
				count += 1.0;
			}
			/* and part of the last one */
			if ((int)next <= halfLen)
			{
				sumRe += input[(int)next] * (next - floor(next));
				count += next - floor(next);
			}

			output[k] = sumRe / count;
		}
		else
			/* The output indices are more densely packed than the input indices
			** so interpolate between input values to generate more output values.
			*/
			/* Take a weighted average of the nearest values */
			output[k] = input[(int)this] * (1.0 - (this - floor(this))) + input[(int)this + 1] * (this - floor(this));
	}
}
static void gridTransformComplex(mxComplexDouble *output, mxComplexDouble *input, double *map, int halfLen, double lowerBoundHz, double upperBoundHz, double samplerate, double linLogRatio)
{
	for (int k = 0; k < halfLen; k++)
	{	/* Average the pixels in the range it comes from */
		double this = gridMap(halfLen, k, lowerBoundHz, upperBoundHz, samplerate, linLogRatio);
		double next = gridMap(halfLen, k + 1, lowerBoundHz, upperBoundHz, samplerate, linLogRatio);
		map[k] = (this + next) * 0.5;

		/* Range check: can happen if --max-freq > samplerate / 2 */
		if (this > halfLen)
		{
			output[k].real = 0.0;
			output[k].imag = 0.0;
			return;
		}

		if (next > this + 1)
		{	/* The output indices are more sparse than the input indices
			** so average the range of input indices that map to this output,
			** making sure not to exceed the input array (0..halfLen inclusive)
			*/
			/* Take a proportional part of the first sample */
			double count = 1.0 - (this - floor(this));
			double sumRe = input[(int)this].real * count;
			double sumIm = input[(int)this].imag * count;

			while ((this += 1.0) < next && (int)this <= halfLen)
			{
				sumRe += input[(int)this].real;
				sumIm += input[(int)this].imag;
				count += 1.0;
			}
			/* and part of the last one */
			if ((int)next <= halfLen)
			{
				sumRe += input[(int)next].real * (next - floor(next));
				sumIm += input[(int)next].imag * (next - floor(next));
				count += next - floor(next);
			}

			output[k].real = sumRe / count;
			output[k].imag = sumIm / count;
		}
		else
		{
			output[k].real = input[(int)this].real * (1.0 - (this - floor(this))) + input[(int)this + 1].real * (this - floor(this));
			output[k].imag = input[(int)this].imag * (1.0 - (this - floor(this))) + input[(int)this + 1].imag * (this - floor(this));
		}
	}
}
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	size_t i, j, m, n;
	// Check for proper number of arguments
#ifndef NOCHECK
	if (nrhs != 5)
		mexErrMsgIdAndTxt("MATLAB:ltv:invalidNumInputs", "5 input arguments are required.");
	else if (nlhs > 2)
		mexErrMsgIdAndTxt("MATLAB:ltv:maxlhs", "Too many output arguments.");
#endif
	// First input is complex spectrum
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	unsigned int len = (unsigned int)MAX(m, n);
#ifndef NOCHECK
	if (!mxIsDouble(prhs[0]) || (MIN(m, n) != 1))
		mexErrMsgIdAndTxt("MATLAB:ltv:incorrectType", "First argument is X");
#endif
	// Other input are parameters
	double a = *(mxGetPr(prhs[1]));
	if (a < DBL_EPSILON)
		a = DBL_EPSILON;
	if (a > (1.0 - DBL_EPSILON))
		a = (1.0 - DBL_EPSILON);
	double lb = *(mxGetPr(prhs[2]));
	double ub = *(mxGetPr(prhs[3]));
	double fs = *(mxGetPr(prhs[4]));
	plhs[1] = mxCreateDoubleMatrix(len, 1, mxREAL);
	double *map = mxGetPr(plhs[1]);
	if (mxIsComplex(prhs[0]))
	{
		plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
		mxComplexDouble *y = mxGetComplexDoubles(plhs[0]);
		mxComplexDouble *X = mxGetComplexDoubles(prhs[0]);
		gridTransformComplex(y, X, map, len, lb, ub, fs, a);
	}
	else
	{
		plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
		double *y = mxGetPr(plhs[0]);
		double *X = mxGetPr(prhs[0]);
		gridTransformReal(y, X, map, len, lb, ub, fs, a);
	}
	return;
}