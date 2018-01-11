### R caller of C code
# x - genotype
# y - phenotype (expression)
# cvrt - matrix of covariates (a column per covariate)

.outputNames = c("beta0", "beta1", "nits", "SSE",
                "SST", "Ftest", "eta", "SE_eta");

effectSizeEstimationC = function(x, y, cvrt){
    stopifnot( length(x) == length(y) );
    stopifnot( (length(cvrt) %% length(x)) == 0 );
    stopifnot( NROW(cvrt) == length(x) );
    stopifnot( NCOL(cvrt)+3 <= length(x) );
    
    if(length(cvrt) > 0){
        cvrt_qr = cbind(1, cvrt);
        cvrt_qr = qr.Q(qr(cvrt_qr));
    } else {
        cvrt_qr = rep(1/sqrt(length(x)), length(x));
    }
    rez = .Call("effectSizeSingleC", x, y, cvrt_qr, PACKAGE = "ACMEeqtl");
    names(rez) = .outputNames;
    return(rez);
}

effectSizeEstimationR = function(x, y, cvrt){
    beta_cur = 0;
    beta_pre = 0;
    beta_dif = 0;
    SSR_pre = .Machine$double.xmax;
    nits = 0;
    nover = 0;
    
    # Prepare C_std
    C_new = cbind(1, cvrt);
    C_new = qr.Q(qr(C_new));
    
    xmin = min(x, -.Machine$double.eps); # Not zero
    xmax = max(x,  .Machine$double.eps); # Not zero
    beta_min = -1 / xmax;
    beta_max = -1 / xmin;
    
    
    for( iterations_counter in 0:300 ){ # iterations_counter = 0

        beta_cur = beta_cur + beta_dif;
        beta_dif = 0;
        
        if(beta_cur <= beta_min){
            beta_cur = (beta_pre + beta_min)/2;
            # message("Distance to min: ", beta_cur - beta_min);
        }
        
        if(beta_cur >= beta_max)
            beta_cur = (beta_pre + beta_max)/2;
        
        fit = 1 + beta_cur * x;

        # Residuals before covariates
        ydiff = y - log(fit);
        # Residuals with covariates
        ydiff_adj = ydiff - C_new %*% crossprod(C_new, ydiff);
        
        # plot(y - C_new %*% crossprod(C_new, y), fit)
        
        # Sum of squared errors
        SSR_cur = crossprod(ydiff_adj);
        # cat("nits",nits,"sse =",SSR_cur,"\n");
        
        # The SST = SSR at first interation
        if(iterations_counter == 0)
            SST = SSR_cur;
        
        # Stop if achieved stable position
        if(abs(SSR_cur / SSR_pre - 1) < 1e-12 * length(x))
            break;
        
        # Did we worsen the fit
        if(SSR_cur > SSR_pre){
            beta_cur = (beta_cur + beta_pre) / 2;
            # cat("Overshoot: SSR_cur - SSR_pre =", SSR_cur - SSR_pre, "\n");
            nover = nover + 1;
            next;
        }
        
        # NLLS step
        xoverfit = x / fit;
        FF = crossprod(xoverfit);
        CF = crossprod(C_new, xoverfit);
        XtX_inv_sub = 1 / (FF - crossprod(CF));
        
        beta_dif = as.vector(XtX_inv_sub * crossprod(xoverfit, ydiff_adj));
        # 
        # beta_pre = beta_cur;
        # beta_cur = beta_cur + as.vector(beta_dif);
        
        
        # change2 = crossprod(beta_dif);
        #if( change2 < 1e-8 )
        #  break;
        #cat(nits, "nits", change2, "change", beta_cur, "beta",
        # round(SSR_cur, 2), "SSR", SSR_cur / SSR_pre - 1, "SSR_diff\n")
        SSR_pre = SSR_cur;
        beta_pre = beta_cur;
    }
    beta0 = exp(mean(y) - mean(log(1 + beta_cur * x)));
    beta1 = beta_cur * beta0;
    
    # SE calculation
    P = diag(length(y)) - tcrossprod(C_new);
    n = NROW(x);
    xoverfit = x / fit;
    H = crossprod(xoverfit^2, crossprod(P, ydiff)) + 
        crossprod(xoverfit, crossprod(P, xoverfit));
    SE_H = sqrt(H^(-1) * SSR_cur / (n - ncol(C_new) - 1)); 
    # c("beta0", "beta1", "nits", "SSE", "SST", "F","eta","SE");
    Ftest = (SST - SSR_cur)/(SSR_cur/(n-ncol(C_new)-1));
    # cat("H = ",H,"\n");
    # cat("Hyx2 = ",crossprod((x/fit)^2, crossprod(P, ydiff)),"\n");
    # cat("Hxx = ",crossprod(x/fit, crossprod(P, x/fit)),"\n");
    rez = c(
        beta0=beta0,
        beta1=beta1,
        nits=iterations_counter,
        SSE=SSR_cur,
        SST = SST,
        Ftest = Ftest,
        eta = beta_cur,
        SE_eta = SE_H);
    stopifnot( names(rez) == .outputNames );
    return(rez)
}
