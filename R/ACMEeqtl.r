### R caller of C code
# x - genotype
# y - phenotype (expression)
# cvrt - matrix of covariates (a column per covariate)

.outputNames = c("beta0", "beta1", "nits", "SSE", "SST", "F","eta","SE_eta");

effectSizeEstimationC = function(x, y, cvrt) {
	stopifnot( length(x) == length(y) );
	stopifnot( (length(cvrt) %% length(x)) == 0 );
	stopifnot( NROW(cvrt) == length(x) );
	stopifnot( NCOL(cvrt)+3 <= length(x) );
	
	if(length(cvrt) > 0) {
		cvrt_qr <- cbind(1, cvrt);
		cvrt_qr = qr.Q(qr(cvrt_qr));
	} else {
		cvrt_qr = rep(1/sqrt(length(x)), length(x));
	}
	rez = .Call("effectSizeSingleC", x, y, cvrt_qr, PACKAGE = "ACMEeqtl");
	names(rez) = .outputNames;
	return(rez);
}

effectSizeEstimationR = function(x, y, cvrt) {
	beta_init = 0
	beta_cur = beta_init;
	beta_pre = -1;
	SSR_pre <- 1.797693e+308
	nits <- 0;
	nover <- 0
	
	# Prepare C_std
	C_new <- cbind(1, cvrt)
	C_new = qr.Q(qr(C_new))
	
	#break_all <- FALSE
	repeat {
		
		#loop_count <- 0    
		
		repeat {
			
			#loop_count <- loop_count + 1
			fit = 1 + beta_cur * x
			if(any(fit <= 1e-10)) {
				#cat("#------LINE89---------beta_cur = ", round(beta_cur, 2), "minfit = ", min(fit), "\n")
				# cat(nits,'beta_cur',beta_cur,'beta_pre',beta_pre,'\n');
				beta_cur = (beta_cur + beta_pre)/2;
				# cat(nits,'beta_cur',beta_cur,'beta_pre',beta_pre,'\n');
			} else {
				break;
			}
			#if (loop_count > 20)
			#break_all <- TRUE
		}
		
		#if (break_all)
		#break
		
		
		ydiff <- y - log(fit)
		ydiff_adj <- ydiff - C_new %*% crossprod(C_new, ydiff)
		SSR_cur <- crossprod(ydiff_adj);
		# cat('nits',nits,'sse =',SSR_cur,'\n');
		if (abs(SSR_cur / SSR_pre - 1) < 1e-15 * length(x))
			break
		if (SSR_cur > SSR_pre) {
			beta_cur <- (beta_cur + beta_pre) / 2
			cat("caught overshoot: SSR_cur - SSR_pre =", SSR_cur - SSR_pre, "\n")
			nover <- nover + 1
			next
		}
		if(nits==0)
			SST = SSR_cur;
		nits <- nits + 1;
		
		FF <- crossprod(x / fit)
		CF <- crossprod(C_new, x / fit)
		XtX_inv_sub <- 1 / (FF - crossprod(CF))
		
		beta_dif <- XtX_inv_sub * crossprod(x / fit, ydiff_adj)
		
		beta_pre = beta_cur;
		beta_cur = beta_cur + beta_dif;
		
		
		change2 = crossprod(beta_dif);
		#if( change2 < 1e-8 )
		#  break;
		#cat(nits, "nits", change2, "change", beta_cur, "beta", round(SSR_cur, 2), "SSR", SSR_cur / SSR_pre - 1, "SSR_diff\n")
		SSR_pre <- SSR_cur
		
	}
	beta0 <- exp(mean(y) - mean(log(1 + beta_cur * x)));
	beta1 <- beta_cur * beta0;
	
	# SE calculation
	P = diag(length(y)) - tcrossprod(C_new);
	n = NROW(x);
	H <- crossprod((x/fit)^2, crossprod(P, ydiff)) + 
		crossprod(x/fit, crossprod(P, x/fit))
	SE_H <- sqrt(H^(-1) * SSR_cur / (n - ncol(C_new) - 1)) 
	# c("beta0", "beta1", "nits", "SSE", "SST", "F","eta","SE");
	F = (SST - SSR_cur)/(SSR_cur/(n-ncol(C_new)-1));
	# cat('H = ',H,'\n');
	# cat('Hyx2 = ',crossprod((x/fit)^2, crossprod(P, ydiff)),'\n');
	# cat('Hxx = ',crossprod(x/fit, crossprod(P, x/fit)),'\n');
	rez = c(beta0=beta0, beta1=beta1, nits=nits, SSE=SSR_cur, SST = SST, F = F, eta = beta_cur, SE_eta = SE_H);
	stopifnot( names(rez) == .outputNames );
	return(rez)
}

