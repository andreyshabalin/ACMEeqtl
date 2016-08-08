#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdbool.h>

# Add use of R_CheckUserInterrupt

const int outrows = 8;
const bool isdebug = false;

// C analog of sum(x)
inline double sum(double *x, int n) {
	double sum = 0;
	for( int i=0; i<n; i++)
		sum += x[i];
	return(sum);
}

// C analog of sum(x^2)
inline double sumsq(double *x, int n) {
	double sum = 0;
	for( int i=0; i<n; i++)
		sum += x[i] * x[i];
	return(sum);
}

// C analog of sum(x*y)
inline double sumxy(double *x, double *y, int n) {
	double sum = 0;
	for( int i=0; i<n; i++)
		sum += x[i] * y[i];
	return(sum);
}

static R_CallMethodDef callMethods[] = {
	{"vecmatProductC", (DL_FUNC) &vecmatProductC, 2},
	{"orthogonolizeVectorC", (DL_FUNC) &orthogonolizeVectorC, 2},
	{NULL, NULL, 0}
};

void R_init_ACMEeqtl(DllInfo *info)	{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
