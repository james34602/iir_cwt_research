#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <mex.h>
#include <float.h>
#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif
#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
#define PI 3.14159265358979323846264338327950288419716939937510582097494459072381640628620899862803482534211706798
#define sqrt2PI 2.50662827463100050241576528
#define IPI4 0.75112554446
typedef struct
{
    float fb;

    float ifb, fb2;
    int width;
    float four_wavelen;
    float* mother;
} Morlet;
typedef struct
{
    float* scales;
    int fs;
    float fourwavl;
    int nscales;
} Scales;
inline int find2power(long long n)
{
    long long m, m2;

    m = 0;
    m2 = 1 << m; /* 2 to the power of m */
    while (m2 - n < 0) {
        m++;
        m2 <<= 1; /* m2 = m2*2 */
    }
    return(m);
}
void InitMorlet(Morlet* mwt, float bandwidth)
{
    mwt->four_wavelen = 0.9876f;
    mwt->fb = bandwidth;
    mwt->fb2 = 2.0f * mwt->fb * mwt->fb;
    mwt->ifb = 1.0f / mwt->fb;
    mwt->mother = NULL;
}

void generate1(Morlet* mwt, int size) {
    //Frequency domain, because we only need size. Default scale is always 2;
    mwt->width = size;

    float tmp1;
    float toradians = (2 * PI) / (float)size;
    float norm = sqrt(2 * PI) * IPI4;

    mwt->mother = (float*)malloc(sizeof(float) * mwt->width);

    //calculate array
    for (int w = 0; w < mwt->width; w++) {
        tmp1 = (2.0f * ((float)w * toradians) * mwt->fb - 2.0f * PI * mwt->fb);
        tmp1 = -(tmp1 * tmp1) / 2;
        mwt->mother[w] = (norm * exp(tmp1));
    }
}
int getSupport(float fb, float scale)
{
    return (int)(fb * scale * 3.0f);
}
static inline void calculate_logscale_array(float* scales, float base, float four_wavl, int fs, float f0, float f1, int fn)
{
    //If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
    float nf0 = f0;
    float nf1 = f1;
    float s0 = (fs / nf1);
    float s1 = (fs / nf0);

    //Cannot pass the nyquist frequency
    assert(("Max frequency cannot be higher than the Nyquist frequency (fs/2)", f1 <= fs / 2));

    float power0 = log(s0) / log(base);
    float power1 = log(s1) / log(base);
    float dpower = power1 - power0;

    for (int i = 0; i < fn; i++) {
        float power = power0 + (dpower / (fn - 1)) * i;
        scales[i] = pow(base, power);
    }
}
static inline void calculate_linfreq_array(float* scales, float four_wavl, int fs, float f0, float f1, int fn) {

    float nf0 = f0;
    float nf1 = f1;
    //If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;

    //Cannot pass the nyquist frequency
    assert(("Max frequency cannot be higher than the Nyquist frequency (fs/2)", f1 <= fs / 2));
    float df = nf1 - nf0;

    for (int i = 0; i < fn; i++) {
        scales[fn - i - 1] = (((float)fs) / (nf0 + (df / fn) * (float)i));
    }
}
void InitScales(Scales* scl, Morlet* wav, unsigned char st, int afs, float af0, float af1, int afn)
{
    scl->fs = afs;
    scl->scales = (float*)malloc(afn * sizeof(float));
    scl->fourwavl = wav->four_wavelen;
    scl->nscales = afn;
    if (!st)
        calculate_logscale_array(scl->scales, 2.0f, wav->four_wavelen, afs, af0, af1, afn);
    else
        calculate_linfreq_array(scl->scales, wav->four_wavelen, afs, af0, af1, afn);
}
void getScales(Scales* scl, float* pfreqs, int pnf) {
    for (int i = 0; i < pnf; i++) {
        pfreqs[i] = scl->scales[i];
    }
}
void getFrequencies(Scales* scl, float* pfreqs, int pnf) {
    for (int i = 0; i < pnf; i++) {
        pfreqs[i] = ((float)scl->fs) / scl->scales[i];
    }
}
unsigned int fftbased(int isize, float scale)
{
    float endpointf = fmin(isize / 2.0, ((isize * 2.0) / scale));
    int endpoint = ((int)endpointf);
    //printf("%d %d\n", 0, endpoint);
    for (int q1 = 0; q1 < endpoint; q1++) {
        //output[q1][0] = input[q1][0] * mother[(int)tmp]; // 1 multiplication
        //output[q1][1] = input[q1][1] * mother[(int)tmp]; // 1 multiplication, ? * (1 - 2 * imaginary);
    }
	return endpoint;
}
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	// Check for proper number of arguments
#ifndef NOCHECK
	if (nrhs != 7)
		mexErrMsgIdAndTxt("MATLAB:fcwtComplexity:invalidNumInputs", "7 input arguments are required.");
	else if (nlhs > 1)
		mexErrMsgIdAndTxt("MATLAB:fcwtComplexity:maxlhs", "Too many output arguments.");
#endif
	// Other input are IIR coefficients
	long long n = (long long)(*mxGetPr(prhs[0]));
	double fs = (double)(*mxGetPr(prhs[1]));
	double f0 = (double)(*mxGetPr(prhs[2]));
    double f1 = (double)(*mxGetPr(prhs[3]));
    long long fn = (long long)(*mxGetPr(prhs[4]));
	unsigned char linLogLinFreq = (unsigned char)(*mxGetPr(prhs[5]));
    if (f1 > fs/2)
		f1 = fs/2;
    //Initialize a Morlet wavelet having sigma=1.0;
    Morlet morl;
    InitMorlet(&morl, 1.0);
    //dist      - FCWT_LOGSCALES | FCWT_LINSCALES for logarithmic or linear distribution of scales across frequency range
    //fs        - sample frequency
    //f0        - beginning of frequency range
    //f1        - end of frequency range
    //fn        - number of wavelets to generate across frequency range
    Scales scs;
    InitScales(&scs, &morl, linLogLinFreq, fs, f0, f1, fn);
    //Find nearest power of 2
    const unsigned long long nt = find2power(n);
    const unsigned long long newsize = 1 << nt;
    //Generate mother wavelet function
    generate1(&morl, newsize);
	plhs[0] = mxCreateNumericMatrix(scs.nscales, 1, mxDOUBLE_CLASS, mxREAL);
	double *Y = (double*)mxGetPr(plhs[0]);
    for (unsigned int i = 0; i < scs.nscales; i++) {
        //FFT-base convolution in the frequency domain
        Y[i] = fftbased(newsize, scs.scales[i]);
    }
    free(morl.mother);
    free(scs.scales);
	return;
}