#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdbool.h>

# Add use of R_CheckUserInterrupt

const int outrows = 8;
const bool isdebug = false;


static R_CallMethodDef callMethods[] = {
	{"vecmatProductC", (DL_FUNC) &vecmatProductC, 2},
	{"orthogonolizeVectorC", (DL_FUNC) &orthogonolizeVectorC, 2},
	{NULL, NULL, 0}
};

void R_init_ACMEeqtl(DllInfo *info)	{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
