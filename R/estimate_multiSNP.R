# Takes a SNP matrix, expression response, and covariates and 
# estimates the corresponding multivariate ACME model via MLE
#
# Arguments:
#  SNPs - a matrix of genotypes; columns are samples
#  expr - a vector of *ACME-transformed* expression values
#  cvrt - a matrix of covariates; columns are samples
#  method - the optim method to use; if NULL, default will be used.
#           Only "BFGS" currently implemented. See ?optim for details.
#
# Returns a list object with components:
#  beta0hat - estimated beta0
#  betahat - estimated SNP coefficients
#  gamhat - estimated cvrt coefficients
#  sigma2hat - estimated noise variance
#  df - degrees of freedom
#-------------------------------------------------------------------------------



estimate_multiSNP <- function (SNPs, expr, cvrt, method = "BFGS") {
  
  # Prepping data
  expr <- as.vector(expr)
  ncov <- nrow(cvrt); nsnp <- nrow(SNPs); n <- length(expr)
  X <- rbind(rep(1, n), SNPs)
  
  # RSS function
  loglik <- function(BetaGam) {
    Beta <- BetaGam[1:(nsnp + 1)]
    Gam <- BetaGam[(nsnp + 2):(nsnp + 1 + ncov)]
    Resid <- expr - log(crossprod(X, Beta)) - crossprod(cvrt, Gam)
    return(sum(Resid^2))
  }
  
  # Gradient function
  Dloglik <- function(BetaGam) {
    Beta <- BetaGam[1:(nsnp + 1)]
    Gam <- BetaGam[(nsnp + 2):(nsnp + 1 + ncov)]
    Resid <- expr - log(crossprod(X, Beta)) - crossprod(cvrt, Gam)
    DBeta <- -2 * crossprod(t(X) / as.vector(crossprod(X, Beta)), Resid)
    DGam <- -2 * cvrt %*% Resid
    return(c(DBeta, DGam))
  }
  
  # Finding start value for beta0
  beta0start <- exp(mean(expr))
  parStart <- c(beta0start, rep(0, nsnp + ncov))
  
  if (method == "BFGS") {
    estimates <- optim(parStart, loglik, gr = Dloglik, method = "BFGS")$par
  } else {
    estimates <- optim(parStart, loglik)$par
  }
  
  # Calculating RSS and sigmahat
  RSSfinal <- loglik(estimates)
  DF <- n - nsnp - ncov - 1
  
  return(list("beta0hat" = estimates[1],
              "betahat" = estimates[2:(nsnp + 1)],
              "gamhat" = estimates[(nsnp + 2):(nsnp + ncov)],
              "sigmahat" = RSSfinal / DF,
              "DF" = DF))
  
}

#-------------------------------------------------------------------------------
# Sample data and usage

if (FALSE) {
  
  n <- 1000; ncov <- 5; p <- 0.2
  beta0 <- 1000; beta1 <- 100; beta2 <- 200
  cvrt <- matrix(rnorm(ncov * n), nrow = ncov)
  SNPs <- matrix(sample(c(0:2), n * 2, replace = TRUE), nrow = 2)
  X <- rbind(rep(1, n), SNPs); Beta = c(beta0, beta1, beta2)
  trueGam <- rnorm(ncov)
  expr <- crossprod(X, Beta) * exp(crossprod(cvrt, trueGam) + rnorm(n)) - 1
  log_expr <- log(1 + expr)
  
  estimate_multiSNP(SNPs, log_expr, cvrt)
  
}

#-------------------------------------------------------------------------------

msACME = function (
		genefm = 'gene',
		snpsfm = 'snps',
		glocfm = 'gene_loc',
		slocfm = 'snps_loc',
		cvrtfm = 'cvrt',
		acmefm = 'ACME',
		workdir = NULL,
		genecap = Inf,
		verbose = TRUE) {
	
	if(!is.null(workdir))
		olddir = setwd(workdir);
	
	### Orthonormalize covariates
	{
		if(verbose)	message('Loading and orthonormalizing covariates');
		cvrt = fm.load(cvrtfm);
		cvrt_qr = qr.Q(qr(cbind(1, cvrt)));
	} # cvrt, cvrt_qr
	
	### Gene/SNP locations	
	{
		if(verbose)	message('Loading gene/SNP locations');
		gene_loc = fm.load(glocfm);
		snps_loc = fm.load(slocfm);
		stopifnot( !is.unsorted(snps_loc) );
	} # gene_loc, snps_loc

	### Get matrix sizes and gene names
	{
		if(verbose)	message('Checking gene/SNP filematrices');
		gfm = fm.open(genefm, readonly = TRUE);
		sfm = fm.open(snpsfm, readonly = TRUE);
		
		gnames = colnames(gfm);
		snames = colnames(sfm);
		
		stopifnot( ncol(gfm) == nrow(gene_loc) );
		stopifnot( ncol(sfm) == nrow(snps_loc) );
		stopifnot( nrow(gfm) == nrow(cvrt) );
		stopifnot( nrow(sfm) == nrow(cvrt) );
		
		n = nrow(cvrt);
		p = ncol(cvrt);
		
		# Keeping gfm/sfm open for analysis loop
		# close(gfm);
		# close(sfm);
	}
	
	### Create output matrix
	{
		if(verbose)	message('Creating output filematrix')
		fm = fm.create(paste0(acmefm, '_multiSNP'), nrow = 4, ncol = 0);
		rownames(fm) = c("geneid","snp_id","eta", "forward_adjR2")
		close(fm)
	}
	
	### Doing multi-SNP estimation
	
	afm = fm.open(acmefm, readonly = TRUE);
	cur_start = 1;
	count = 1;
	
	while (count <= genecap && cur_start <= ncol(gfm)) {
		
		cur_gene = as.numeric(afm[1, cur_start]);
		expr = gfm[ , cur_gene];
		
		# Scanning afm for cur_gene
		if (verbose) message(paste0('Scanning for gene', count));
		cur_end = cur_start + 10 - 1;
		scan_count = 1;
		while(sum(as.vector(afm[1, cur_start:cur_end]) == cur_gene) == 
			  10 * scan_count) {
			cur_end = cur_end + 10;
			scan_count = scan_count + 1;
		}
		
		# Getting cur_gene locations and associated snps
		glocs = c(cur_start:cur_end)[afm[1, cur_start:cur_end] == cur_gene];
		afm_submat = afm[ , glocs];
		gsnps = afm_submat[2, ];
		
		# Initializing inner loop
		last_snp = 0;
		snp_list = integer(0);
		final_AR2 = numeric(0);
		remainingSNPs = 1:ncol(afm_submat);
		SST = afm_submat[7, 1];
		SSEs = afm_submat[6, ];
		adjR2 = 1 - (SSEs / SST) / ((n - p - length(snp_list) - 2) / (n - p - 1));
		old_AR2 = -Inf
		new_AR2 = max(adjR2);
		
		# Forward selection loop
		while (new_AR2 >= old_AR2 && length(snp_list) < n - p - 2) {
			
			next_snp = remainingSNPs[which.max(adjR2)];
			snp_list = c(snp_list, next_snp);
			final_AR2 = c(final_AR2, new_AR2);
			remainingSNPs = setdiff(remainingSNPs, next_snp);
			
			if (verbose) message(paste0('----Gene ', count, ', added SNP ', next_snp,
										' with AdjR2 ', round(max(adjR2), 4)));
			
			# Computing SSEs
			real_snp_list = afm_submat[2, snp_list];
			SNPmat = rbind(t(sfm[ , real_snp_list]), 0);
			SSEs = numeric(length(remainingSNPs));
			for (rsi in seq_along(remainingSNPs)) {
				
				real_snps = afm_submat[2, c(snp_list, remainingSNPs[rsi])];
				SNPmat[nrow(SNPmat) , ] = 
					as.vector(sfm[ , afm_submat[2, remainingSNPs[rsi]]]);
				suppressWarnings(
					{msacme_out = estimate_multiSNP(SNPmat / 1000, as.vector(expr), t(cvrt_qr))}
				);
				
				SSEs[rsi] = msacme_out$DF * msacme_out$sigmahat;
				
			}
			
			adjR2 = 1 - (SSEs / SST) / ((n - p - length(snp_list) - 2) / (n - p - 1));
			next_snp = remainingSNPs[which.max(adjR2)];
			old_AR2 = new_AR2;
			new_AR2 = max(adjR2);
			
		}
			
		# Computing final estimates and writing the results
		SNPmat = subset(SNPmat, c(rep(T, nrow(SNPmat) - 1), F))
		suppressWarnings(
			{final_msacme_out = estimate_multiSNP(SNPmat / 1000, as.vector(expr), t(cvrt))}
		);
		real_snps = afm_submat[2, snp_list];
		outmatrix = rbind(rep(cur_gene, length(real_snps)),
						  real_snps,
						  final_msacme_out$betahat / final_msacme_out$beta0hat,
						  final_AR2);
		fm = fm.open(paste0(acmefm, '_multiSNP'));
		fm$appendColumns(outmatrix);
		close(fm);
		
		# Re-setting loop
		cur_start = tail(glocs, 1) + 1;
		count = count + 1;
		
	}
	
}
