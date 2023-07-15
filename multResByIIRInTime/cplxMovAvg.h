#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct
{
	double real, imag;
} cplx;
size_t estCplxMovAvg(unsigned int *wndSize, unsigned int nBands)
{
	unsigned int i;
	unsigned int totalWndSize = 0;
	for (i = 0; i < nBands; i++)
		totalWndSize += wndSize[i];
	return sizeof(unsigned int) * 2 + sizeof(unsigned int) * nBands * 3 + sizeof(double) * nBands + sizeof(cplx) * nBands * 3 + sizeof(cplx) * totalWndSize;
}
void initCplxMovAvg(char *me, unsigned int *wndSize2, unsigned int nBands)
{
	unsigned int i;
	*(unsigned int*)me = nBands;
	unsigned int *wndSize = (unsigned int*)(me + sizeof(unsigned int) * 2);
	unsigned int *writeIdx = wndSize + nBands;
	unsigned int *Pos = writeIdx + nBands;
	double *g = (double*)(Pos + nBands);
	cplx *firstInput = (cplx*)(g + nBands);
	cplx *dl = firstInput + nBands;
	cplx *acc = dl + nBands;
	unsigned int totalWndSize = 0;
	for (i = 0; i < nBands; i++)
	{
		wndSize[i] = wndSize2[i];
		Pos[i] = totalWndSize;
		totalWndSize += wndSize[i];
		printf("%llu %llu\n", i ,totalWndSize);
		firstInput[i].real = 0.0;
		firstInput[i].imag = 0.0;
		dl[i].real = 0.0;
		dl[i].imag = 0.0;
		acc[i].real = 0.0;
		acc[i].imag = 0.0;
		writeIdx[i] = 0;
		g[i] = 1.0 / wndSize[i];
	}
	*(((unsigned int*)me) + 1) = totalWndSize;
	cplx *buf = acc + nBands;
	printf("%llu\n", totalWndSize);
	memset(buf, 0, totalWndSize * sizeof(cplx));
}
void procCplxMovAvgCplx(char *me, cplx *x, cplx *y)
{
	unsigned int i;
	unsigned int nBands = *(unsigned int*)me;
	unsigned int totalWndSize = *(unsigned int*)(me + sizeof(unsigned int));
	unsigned int *wndSize = (unsigned int*)(me + sizeof(unsigned int) * 2);
	unsigned int *writeIdx = wndSize + nBands;
	unsigned int *Pos = writeIdx + nBands;
	double *g = (double*)(Pos + nBands);
	cplx *firstInput = (cplx*)(g + nBands);
	cplx *dl = firstInput + nBands;
	cplx *acc = dl + nBands;
	cplx *buf = acc + nBands;
	for (i = 0; i < nBands; i++)
	{
        if (wndSize[i] <= 1)
        {
            y[i].real = x[i].real;
            y[i].imag = x[i].imag;
            continue;
        }
		unsigned int pos = Pos[i];
		unsigned int readIdx = writeIdx[i] + 1;
		if (readIdx >= wndSize[i])
			readIdx = 0;
		dl[i].real = firstInput[i].real;
		dl[i].imag = firstInput[i].imag;
		firstInput[i].real = buf[pos + readIdx].real;
		firstInput[i].imag = buf[pos + readIdx].imag;
		buf[pos + writeIdx[i]].real = x[i].real;
		buf[pos + writeIdx[i]].imag = x[i].imag;
		writeIdx[i]++;
		if (writeIdx[i] >= wndSize[i])
			writeIdx[i] = 0;
		acc[i].real = acc[i].real - dl[i].real + x[i].real;
		acc[i].imag = acc[i].imag - dl[i].imag + x[i].imag;
		y[i].real = acc[i].real * g[i];
		y[i].imag = acc[i].imag * g[i];
	}
}
void procCplxMovAvgReal(char *me, double *x, cplx *y)
{
	unsigned int i;
	unsigned int nBands = *(unsigned int*)me;
	unsigned int totalWndSize = *(unsigned int*)(me + sizeof(unsigned int));
	unsigned int *wndSize = (unsigned int*)(me + sizeof(unsigned int) * 2);
	unsigned int *writeIdx = wndSize + nBands;
	unsigned int *Pos = writeIdx + nBands;
	double *g = (double*)(Pos + nBands);
	cplx *firstInput = (cplx*)(g + nBands);
	cplx *dl = firstInput + nBands;
	cplx *acc = dl + nBands;
	cplx *buf = acc + nBands;
	for (i = 0; i < nBands; i++)
	{
        if (wndSize[i] <= 1)
        {
            y[i].real = x[i];
            y[i].imag = 0;
            continue;
        }
		unsigned int pos = Pos[i];
		unsigned int readIdx = writeIdx[i] + 1;
		if (readIdx >= wndSize[i])
			readIdx = 0;
		dl[i].real = firstInput[i].real;
		dl[i].imag = firstInput[i].imag;
		firstInput[i].real = buf[pos + readIdx].real;
		firstInput[i].imag = buf[pos + readIdx].imag;
		buf[pos + writeIdx[i]].real = x[i];
		buf[pos + writeIdx[i]].imag = 0;
		writeIdx[i]++;
		if (writeIdx[i] >= wndSize[i])
			writeIdx[i] = 0;
		acc[i].real = acc[i].real - dl[i].real + x[i];
		acc[i].imag = acc[i].imag - dl[i].imag;
		y[i].real = acc[i].real * g[i];
		y[i].imag = acc[i].imag * g[i];
	}
}