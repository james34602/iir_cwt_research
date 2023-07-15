#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
typedef struct
{
	double real, imag;
} cplx;
typedef struct
{
	cplx *inputs;
	long nBands, *delayLineLength, *offset, *inPoint, *outPoint;
	double *alpha;
} fractionalDelayLine;
void fractionalDelayLine_free(fractionalDelayLine *p)
{
	free(p->inputs);
	free(p->delayLineLength);
	free(p->offset);
	free(p->inPoint);
	free(p->outPoint);
	free(p->alpha);
}
void fractionalDelayLineInit(fractionalDelayLine *p, unsigned int nBands, double *lag, double initValue)
{
	p->nBands = nBands;
	p->delayLineLength = (long*)malloc(p->nBands * sizeof(long));
	p->offset = (long*)malloc(p->nBands * sizeof(long));
	p->inPoint = (long*)malloc(p->nBands * sizeof(long));
	p->outPoint = (long*)malloc(p->nBands * sizeof(long));
	p->alpha = (double*)malloc(p->nBands * sizeof(double));
	unsigned int totalLen = 0;
	for (unsigned int i = 0; i < p->nBands; i++)
	{
		if (lag[i] < 0.0)
			p->delayLineLength[i] = (long)ceilf(0.0) + 1;
		else
			p->delayLineLength[i] = (long)ceilf(lag[i]) + 1;
		p->offset[i] = totalLen;
		totalLen += p->delayLineLength[i];
		p->inPoint[i] = 0;
		p->outPoint[i] = p->delayLineLength[i] >> 1;
	}
	p->inputs = (cplx*)malloc(totalLen * sizeof(cplx));
	for (unsigned int i = 0; i < totalLen; i++)
	{
        p->inputs[i].real = initValue;
        p->inputs[i].imag = 0.0;
    }
}
void fractionalDelayLine_setDelay(fractionalDelayLine *p, double *lag)
{
	for (unsigned int i = 0; i < p->nBands; i++)
	{
		double outputPointer;
		if (lag[i] < 0.0)
			outputPointer = p->inPoint[i]; /* read chases write, + 1 for interp. */
		else
			outputPointer = p->inPoint[i] - lag[i]; /* read chases write, + 1 for interp. */
		while (outputPointer < 0.0)
			outputPointer += (double)p->delayLineLength[i];           /* modulo maximum length */
		while (outputPointer >= (double)p->delayLineLength[i])
			outputPointer -= p->delayLineLength[i];           /* modulo maximum length */
		p->outPoint[i] = (long)outputPointer;            /* integer part */
		p->alpha[i] = outputPointer - (double)p->outPoint[i]; /* fractional part */
	}
}
void fractionalDelayLineProcessCplx(fractionalDelayLine *p, cplx *x, cplx *y)
{
	for (unsigned int i = 0; i < p->nBands; i++)
	{
		cplx *buf = &p->inputs[p->offset[i]];
		buf[p->inPoint[i]].real = x[i].real;
		buf[p->inPoint[i]++].imag = x[i].imag;
		if (p->inPoint[i] == p->delayLineLength[i])
			p->inPoint[i] -= p->delayLineLength[i];
		double omAlpha = 1.0 - p->alpha[i];
		y[i].real = buf[p->outPoint[i]].real * omAlpha;
		y[i].imag = buf[p->outPoint[i]++].imag * omAlpha;
		if (p->outPoint[i] < p->delayLineLength[i])
		{
			y[i].real += buf[p->outPoint[i]].real * p->alpha[i];
			y[i].imag += buf[p->outPoint[i]].imag * p->alpha[i];
		}
		else
		{
			y[i].real += buf[0].real * p->alpha[i];
			y[i].imag += buf[0].imag * p->alpha[i];
			p->outPoint[i] -= p->delayLineLength[i];
		}
	}
}
void fractionalDelayLineProcessReal(fractionalDelayLine *p, double *x, cplx *y)
{
	for (unsigned int i = 0; i < p->nBands; i++)
	{
		cplx *buf = &p->inputs[p->offset[i]];
		buf[p->inPoint[i]].real = x[i];
		buf[p->inPoint[i]++].imag = 0.0;
		if (p->inPoint[i] == p->delayLineLength[i])
			p->inPoint[i] -= p->delayLineLength[i];
		double omAlpha = 1.0 - p->alpha[i];
		y[i].real = buf[p->outPoint[i]].real * omAlpha;
		y[i].imag = buf[p->outPoint[i]++].imag * omAlpha;
		if (p->outPoint[i] < p->delayLineLength[i])
		{
			y[i].real += buf[p->outPoint[i]].real * p->alpha[i];
			y[i].imag += buf[p->outPoint[i]].imag * p->alpha[i];
		}
		else
		{
			y[i].real += buf[0].real * p->alpha[i];
			y[i].imag += buf[0].imag * p->alpha[i];
			p->outPoint[i] -= p->delayLineLength[i];
		}
	}
}