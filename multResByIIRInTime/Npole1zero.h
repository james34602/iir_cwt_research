#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct
{
	double real, imag;
} cplx;
size_t estNpole1Zero(unsigned int nBands, unsigned int order)
{
	return sizeof(unsigned int) * 2 + sizeof(cplx) * nBands * order;
}
void initNpole1Zero(char *me, unsigned int nBands, unsigned int order)
{
	*(unsigned int*)me = nBands;
	*(((unsigned int*)me) + 1) = order;
	cplx *Z = (cplx*)(me + sizeof(unsigned int) * 2);
	memset(Z, 0, nBands * order * sizeof(cplx));
}
void procNpole1ZeroCplx(char *me, double *b, double *aPacked, cplx *x, cplx *y)
{
	unsigned int i, j;
	unsigned int nBands = *(unsigned int*)me;
	unsigned int order = *(unsigned int*)(me + sizeof(unsigned int));
	cplx *ZPacked = (cplx*)(me + sizeof(unsigned int) * 2);
	for (i = 0; i < nBands; i++)
	{
		cplx *Z = &ZPacked[i * order];
		double *a = &aPacked[i * (order + 1)];
		double YiRe = x[i].real + Z[0].real;
		double YiIm = x[i].imag + Z[0].imag;
		for (j = 1; j < order; j++) // Update conditions
		{
			Z[j - 1].real = Z[j].real - a[j] * YiRe;
			Z[j - 1].imag = Z[j].imag - a[j] * YiIm;
		}
		Z[order - 1].real = -a[order] * YiRe;
		Z[order - 1].imag = -a[order] * YiIm;
		y[i].real = b[i] * YiRe;
		y[i].imag = b[i] * YiIm;
	}
}
void procNpole1ZeroReal(char *me, double *b, double *aPacked, double *x, cplx *y)
{
	unsigned int i, j;
	unsigned int nBands = *(unsigned int*)me;
	unsigned int order = *(unsigned int*)(me + sizeof(unsigned int));
	cplx *ZPacked = (cplx*)(me + sizeof(unsigned int) * 2);
	for (i = 0; i < nBands; i++)
	{
		cplx *Z = &ZPacked[i * order];
		double *a = &aPacked[i * (order + 1)];
		double YiRe = x[i] + Z[0].real;
		double YiIm = Z[0].imag;
		for (j = 1; j < order; j++) // Update conditions
		{
			Z[j - 1].real = Z[j].real - a[j] * YiRe;
			Z[j - 1].imag = Z[j].imag - a[j] * YiIm;
		}
		Z[order - 1].real = -a[order] * YiRe;
		Z[order - 1].imag = -a[order] * YiIm;
		y[i].real = b[i] * YiRe;
		y[i].imag = b[i] * YiIm;
	}
}