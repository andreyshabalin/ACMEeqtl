#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdbool.h>
#include <R_ext/Rdynload.h>

// Add use of R_CheckUserInterrupt

// Length of output vector
const int outrows = 8;
// Show debug messages
const bool isdebug = false;

// Product of vector 'xx' and matrix 'cc'
// xx - n-vector
// cc - n*ncvrt matrix
// out - ncvrt vector
void vecmatProduct(double *xx, double *cc, double *out, int n, int ncvrt) {
	for( int j=0; j<ncvrt; j++) {
		double rez = 0;		
		for( int i=0; i<n; i++) {
			rez += xx[i] * cc[i + j*n];
		}
		out[j] = rez;
	}
}

// R wrapper for vecmatProduct, for testing
SEXP vecmatProductC(SEXP x, SEXP cvrtqr){
	x = PROTECT(coerceVector(x, REALSXP));
	cvrtqr = PROTECT(coerceVector(cvrtqr, REALSXP));
	
	int n = length(x);
	int ncvrt = length(cvrtqr) / length(x);
	
	SEXP rez = PROTECT(allocVector(REALSXP, ncvrt)); 
	
	vecmatProduct(REAL(x), REAL(cvrtqr), REAL(rez), n, ncvrt);
	
	UNPROTECT(3);
	return(rez);
}


// Orthogonolize 'xx' with respect to 'cc'
// xx - n-vector
// cc - n*ncvrt matrix
// temp - ncvrt vector for temporary data (inner product)
// xxout - output 
void orthogonolizeVector(double *xx, double *cc, double *temp, double *xxout, int n, int ncvrt) {
	// temp = xx %*% cc
	for( int j=0; j<ncvrt; j++) {
		double val = 0;
		for( int i=0; i<n; i++)
			val += xx[i] * cc[i + j*n];
		temp[j] = val;
	}
	// xxout = temp - cc %*% temp
	for( int i=0; i<n; i++) {					
		double value = xx[i];
		for( int j=0; j<ncvrt; j++)
			value -= temp[j] * cc[i + j*n];
		xxout[i] = value;
	}
}

// Orthogonalize vector 'x' with respect to orthonormalized vectors in 'cvrtqr'
SEXP orthogonolizeVectorC(SEXP x, SEXP cvrtqr){
	x = PROTECT(coerceVector(x, REALSXP));
	cvrtqr = PROTECT(coerceVector(cvrtqr, REALSXP));
	
	int n = length(x);
	int ncvrt = length(cvrtqr) / length(x);
	
	SEXP tmp = PROTECT(allocVector(REALSXP, ncvrt)); 
	SEXP rez = PROTECT(allocVector(REALSXP, n)); 
	
	orthogonolizeVector(REAL(x), REAL(cvrtqr), REAL(tmp), REAL(rez), n, ncvrt); 

	UNPROTECT(4);
	return(rez);
}

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


void internalEstimate(double *xx, double *yy, double *cc, 
                      int ncvrt, int n, 
                      double *fit, double *ydiff, double *ycvrt,  
                      double* out){
	double beta_cur = 0;
	double beta_pre = 0;
	double sse_pre = 1.797693e+308;
	double sse_cur = 0;
	int iterations_counter = 0;
	int overshoot_counter = 0;
	double SST = -1;
	
	if(isdebug)Rprintf("n = %d\n", n);
	if(isdebug)Rprintf("ncvrt = %d\n", ncvrt);
	
	double xmin = -2.220446e-16; // Not zero
	double xmax =  2.220446e-16; // Not zero
	for( int i=0; i<n; i++) {
		if(xx[i] < xmin)
			xmin = xx[i];
		if(xx[i] > xmax)
			xmax = xx[i];
	}
	if(isdebug) Rprintf("xmin = %f\n", xmin);
	if(isdebug) Rprintf("xmax = %f\n", xmax);
	
	// 1 + beta * xx > 0
	// beta * xx > -1
	
	// xx > 0: beta > -1 / xx
	// xx < 0: beta < -1 / xx
	
	double beta_min = -1 / xmax;
	double beta_max = -1 / xmin;
	
	if(isdebug) Rprintf("beta_min = %f\n", beta_min);
	if(isdebug) Rprintf("beta_max = %f\n", beta_max);
	
	double beta_dif = 0;
	for( iterations_counter=0; iterations_counter<300; iterations_counter++) {
		
		if(isdebug) Rprintf("iterations_counter = %d\n", iterations_counter);
		
		beta_cur = beta_cur + beta_dif;
		beta_dif = 0;
		
		if(beta_cur <= beta_min)
			beta_cur = (beta_pre + beta_min)/2;
		
		if(beta_cur >= beta_max)
			beta_cur = (beta_pre + beta_max)/2;
		
		if(isdebug) Rprintf("beta_cur = %f\n", beta_cur);
		
		// Calculate fit
		for( int i=0; i<n; i++) {
			fit[i] = 1 + beta_cur * xx[i];
			ydiff[i] = yy[i] - log1p(beta_cur * xx[i]); // ydiff[i] = yy[i] - log(fit[i]);
		}
		
		orthogonolizeVector(ydiff, cc, ycvrt, ydiff, n, ncvrt);
		
		if(isdebug) Rprintf("ydiff tail: %f\n", ydiff[n-1]);
		// void orthogonolizeVector(double *xx, double *cc, double *temp, double *xxout, int n, int ncvrt) {
		
		sse_cur = sumsq(ydiff, n);
		if(isdebug) Rprintf("sse_cur = %f\n", sse_cur);
		
		if(iterations_counter==0)
			SST = sse_cur;
		
		if( fabs(sse_cur / sse_pre - 1) < (1e-12 * n) ) {
			if(isdebug) Rprintf("sse_cur / sse_pre = %f / %f = 1\n", sse_cur, sse_pre);
			break;
		}
		
		if (sse_cur > sse_pre) {
			beta_cur = (beta_cur + beta_pre) / 2;
			overshoot_counter++;
			if(isdebug) Rprintf("iter = %d, overshoot_counter = %d,  %f, %f, %e\n", iterations_counter, overshoot_counter, beta_cur, beta_pre, beta_cur-beta_pre);
			continue;
		}
		
		
		// temp:
		// fit = 1 + beta*x
		// ydiff - length n
		// ycvrt - length ncvrt
		
		// xfit = x / fit
		double *xfit = fit;
		for( int i=0; i<n; i++)
			xfit[i] = xx[i] / fit[i];
		double ff = sumsq(xfit, n);			// FF <- crossprod(x / fit)
		double xfyadj = sumxy(xfit, ydiff, n);	// xfyadj = crossprod(x / fit, ydiff_adj)
		
		if(isdebug) Rprintf("xfit tail: %f\n", xfit[54-1]);
		if(isdebug) Rprintf("ff: %f\n", ff);
		if(isdebug) Rprintf("xfyadj: %f\n", xfyadj);
		
		// ycvrt = crossprod(cvrtqr, x / fit)
		vecmatProduct(xfit, cc, ycvrt, n, ncvrt);
		if(isdebug) Rprintf("ycvrt tail: %f\n", ycvrt[ncvrt-1]);
		
		// XtX_inv_sub <- 1 / (FF - crossprod(CF))
		double XtX_inv_sub = 1 / (ff - sumsq(ycvrt, ncvrt));
		if(isdebug) Rprintf("XtX_inv_sub: %f\n", XtX_inv_sub);
		
		beta_dif = XtX_inv_sub * xfyadj;	// beta_dif == XtX_inv_sub * crossprod(x / fit, ydiff_adj)
		if(isdebug) Rprintf("beta_dif: %f\n", beta_dif);
		
		
		if(isdebug) Rprintf("iter = %d, sse %f -> %f, %e\n", iterations_counter, sse_cur, sse_pre, sse_pre-sse_cur);
		
		sse_pre = sse_cur;
		beta_pre = beta_cur;
	}
	
	if(isdebug) Rprintf("iter = %d\n", iterations_counter);
	if(isdebug) Rprintf("beta_pre = %f\n", beta_pre);
	if(isdebug) Rprintf("beta_cur = %f\n", beta_cur);
	//if(isdebug) Rprintf("beta_dif = %f\n", beta_dif);
	
	
	
	// sumlogfit += log(1 + beta_cur * xx[i]);
	double sumlogfit = 0;
	for( int i=0; i<n; i++)
		sumlogfit += log1p(beta_cur * xx[i]);
	double beta0 = exp(sum(yy,n)/n - sumlogfit/n);		// beta0 = exp(mean(y) - mean(log(1 + beta_cur * x)))
	double beta1 = beta_cur * beta0;
	
	if(isdebug) Rprintf("beta0 = %f\n", beta0);
	if(isdebug) Rprintf("beta1 = %f\n", beta1);
	
	//double *fit,  *ydiff,  *ycvrt;
	double H = 0;
	{	
		// x1Px1
		double *xfit = fit; // xfit = x / fit
		for( int i=0; i<n; i++)
			xfit[i] = xx[i] / (1 + beta_cur * xx[i]);
		orthogonolizeVector(xfit, cc, ycvrt, ydiff, n, ncvrt);
		
		H = sumsq(ydiff, n);
		if(isdebug) Rprintf("x1Px1 = %f\n", H);
		
		// (x12)Py1
		for( int i=0; i<n; i++) {
			ydiff[i] = yy[i] - log1p(beta_cur * xx[i]);	// y1 = y - log( fit )
			xfit[i] = xfit[i] * xfit[i];
		}
		orthogonolizeVector(xfit, cc, ycvrt, xfit, n, ncvrt);
		
		H += sumxy(xfit, ydiff, n);
		if(isdebug) Rprintf("H = %f\n", H);
	}
	
	double SE_eta = sqrt((sse_cur/ (n - ncvrt - 1)) / H);
	
	out[0] = beta0;
	out[1] = beta1;
	out[2] = iterations_counter;
	out[3] = sse_cur;
	out[4] = SST;
	out[5] = (SST - sse_cur)/(sse_cur/(n-ncvrt-1));
	out[6] = beta_cur;
	out[7] = SE_eta;
	
	return;
}

SEXP effectSizeSingleC(SEXP x, SEXP y, SEXP cvrtqr){
	x = PROTECT(coerceVector(x, REALSXP));
	y = PROTECT(coerceVector(y, REALSXP));
	cvrtqr = PROTECT(coerceVector(cvrtqr, REALSXP));
	double *xx = REAL(x);
	double *yy = REAL(y);
	double *cc = REAL(cvrtqr);
	int n = length(x);
	int ncvrt = length(cvrtqr) / length(x);
	//SEXP cvrtdim;
	//PROTECT(cvrtdim=getAttrib(cvrtqr,R_DimSymbol));
	
	double *fit   = (double*)malloc(n*sizeof(double));
	double *ydiff = (double*)malloc(n*sizeof(double));
	double *ycvrt = (double*)malloc(ncvrt*sizeof(double));
	
	SEXP rez = PROTECT(allocVector(REALSXP, outrows)); 
	
	internalEstimate(xx, yy, cc, ncvrt, n, fit, ydiff, ycvrt, REAL(rez));
	
	UNPROTECT(4);
	free(fit);
	free(ydiff);
	free(ycvrt);
	return rez;
}


SEXP effectSizeManyC(SEXP genemat, SEXP snpsmat, SEXP snpsmat_fr,SEXP snpsmat_to, SEXP cvrt_qr, SEXP outmatrix){
	
	genemat = PROTECT(coerceVector(genemat, REALSXP));
	snpsmat = PROTECT(coerceVector(snpsmat, REALSXP));
	snpsmat_fr = PROTECT(coerceVector(snpsmat_fr, REALSXP));
	snpsmat_to = PROTECT(coerceVector(snpsmat_to, REALSXP));
	cvrt_qr = PROTECT(coerceVector(cvrt_qr, REALSXP));
	outmatrix = PROTECT(coerceVector(outmatrix, REALSXP));
	
	double *s_fr = REAL(snpsmat_fr);
	double *s_to = REAL(snpsmat_to);
	double *out = REAL(outmatrix);
	double *xgene = REAL(genemat);
	double *xsnps = REAL(snpsmat);
	double *cc = REAL(cvrt_qr);
	
	int n, ng, ncvrt;
	
	SEXP genedim = getAttrib(genemat,R_DimSymbol);
	if(genedim == R_NilValue) {
		n = length(genemat);
		ng = 1;
	} else {
		genedim = PROTECT(coerceVector(genedim,INTSXP));
		n  = INTEGER(genedim)[0];
		ng = INTEGER(genedim)[1];
		UNPROTECT(1);
	}
	
	//ns = length(snpsmat) / n;
	ncvrt = length(cvrt_qr) / n;
	
	double *fit   = (double*)malloc(n*sizeof(double));
	double *ydiff = (double*)malloc(n*sizeof(double));
	double *ycvrt = (double*)malloc(ncvrt*sizeof(double));
	
	int counter = 0;
	
	for( int i=0; i<ng; i++) {
		for( int j=s_fr[i]; j<=s_to[i]; j++) {
			out[(outrows+2)*counter] = i;
			out[(outrows+2)*counter+1] = j;
			internalEstimate(xsnps + n*j, xgene + n*i, cc, ncvrt, n, fit, ydiff, ycvrt, out + (outrows+2)*counter + 2);
			counter++;
			//Rprintf("i,j  = %d,%d,%d,%d\n", i,j,s_fr[i],s_to[i]);
		}
	}
	
	free(fit);
	free(ydiff);
	free(ycvrt);	

//	Rprintf("n  = %d\n", n);
//	Rprintf("ng = %d\n", ng);
//	Rprintf("ns = %d\n", ns);
//	Rprintf("ncvrt = %d\n", ncvrt);
//	Rprintf("pairs = %d\n", counter);

	UNPROTECT(6);
	return genedim;
}



static R_CallMethodDef callMethods[] = {
	{"vecmatProductC", (DL_FUNC) &vecmatProductC, 2},
	{"orthogonolizeVectorC", (DL_FUNC) &orthogonolizeVectorC, 2},
	{"effectSizeSingleC", (DL_FUNC) &effectSizeSingleC, 3},
	{"effectSizeManyC", (DL_FUNC) &effectSizeManyC, 6},
	{NULL, NULL, 0}
};

void R_init_ACMEeqtl(DllInfo *info)	{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
