# Make a toy eQTL data set for a single tissue
# dependencies: filematrix
# 
# Arguments:
#   n - number of samples
#   ngene - number of genes (default = 5000)
#   nsnp - number of snps (default = 500000)
#   ncvrt - number of covariates (default = 10)
#   p - reference allele probability (default = .2)
#   saveDir - directory in which to save data files (default = cwd)
#             (function will make a subdirectory named "filematrices")
#   returnData - whether or not to return data in a list (default = FALSE)
#                if TRUE, data is returned in a list with components:
#                  'snps'    - toy SNP data
#                  'genes'   - toy gene data
#                  'cvrts'   - toy covariate data
#                  'snploc'  - dummy SNP locations
#                  'geneloc' - dummy gene locations
#   fmat - if TRUE, will save data to filematrices. Otherwise, saves to 
#          text files. (default = FALSE)
#   verbose - whether or not to print SNP filematrix writing progress
#             (default = FALSE)
#-------------------------------------------------------------------------------

create_artificial_data = function(
		nsample, 
		ngene = 5000, 
		nsnp = 500000, 
		ncvrt = 10, 
		minMAF = 0.2, 
		saveDir = getwd(), 
		returnData = FALSE, 
		savefmat = FALSE,
		savetxt = FALSE,
		verbose = TRUE){ 
	
	stopifnot( any(returnData, savefmat, savetxt) )
	stopifnot( nsample > ncvrt )

	# Making SNP data
	{
		if(verbose)
			message('Generating SNP matrix');
		# tic = proc.time();
		MAF = runif(n = nsample, min = minMAF, max = 0.5);
		snps = rbinom(n = nsnp * nsample, size = 2, prob = MAF);
		dim(snps) = c(nsnp, nsample);
		rm(MAF);
		# toc = proc.time();
	  	# show(toc-tic);
	} # snps
		
	# Making covariates and covariate effect
	{
		if(verbose)
			message('Generating covariates');
		cvrts = rnorm(nsample * ncvrt);
		dim(cvrts) = c(ncvrt, nsample);
		cvrtEffect = as.vector(crossprod(cvrts, rnorm(ncvrt)));
	} # cvrts, cvrtEffect
	
	# Generating genes
	{
		if(verbose)
			message('Generating gene expression matrix');
		# Making eQTL network
		
		causalSnps = as.integer(seq(1,nsnp, length.out = ngene))
		
		# Making effect sizes
		etas = (exp(rnorm(ngene)) - 1) / 2;
		
		# Making gene data
		# constant + snp effect + covariateEffect + noise
		genes = log(100) +
		        log(1 + etas * snps[causalSnps, , drop = FALSE]) +
		  		  rep(cvrtEffect, each = ngene) +
		   	  rnorm(ngene * nsample);
		tbleta = data.frame("gene" = seq_len(ngene), "snp" = causalSnps, eta = etas)
		rm(etas, causalSnps)
	} # genes, tbleta
	
	# Naming gene/snp rows and making gene/snp loc dummy dataframes
	{
		rownames(genes) = sprintf('Gene_%04d',1:ngene);
		rownames(snps)  = sprintf('SNP_%06d',1:nsnp);
		rownames(cvrts) = sprintf('Cvrt_%02d',1:ncvrt);
		colnames(genes) = sprintf('Sample_%04d',1:nsample);
		colnames(snps)  = sprintf('Sample_%04d',1:nsample);
		colnames(cvrts) = sprintf('Sample_%04d',1:nsample);
	}
	
	# Generate location tables
	{
		genepos = as.integer(seq(500, 250e6, length.out = ngene))
		geneloc = data.frame("geneid" = rownames(genes),
		                     "chr" = "chr1",
		                     "start" = genepos,
									"end" = genepos + 100);
		snpsloc <- data.frame("SNP" = rownames(snps),
		                      "chr" = "chr1",
		                      "pos" = as.integer(seq(500, 250e6, length.out = nsnp)));
	} # genepos, geneloc, snpsloc
	
	if(savetxt){
		# Making directory for txtfiles and going there
		if(!dir.exists(file.path(saveDir, "txtfiles")))
			dir.create(file.path(saveDir, "txtfiles"), recursive = TRUE)

		if(verbose) message('Saving gene/SNP location files (txt)');
		write.table(geneloc, file = file.path(saveDir, "txtfiles", "gene_loc.txt"),
		          quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t");
		write.table(snpsloc, file = file.path(saveDir, "txtfiles", "snps_loc.txt"),
		          quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t");
		
		if(verbose) message('Saving expression matrix (txt)');
		write.table(genes, file = file.path(saveDir, "txtfiles", "gene.txt"),
		          quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t");
		if(verbose) message('Saving covariates (txt)');
		write.table(cvrts, file = file.path(saveDir, "txtfiles", "cvrt.txt"),
		          quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t");
		if(verbose) message('Saving true effect info (txt)');
		write.table(tbleta, file = file.path(saveDir, "txtfiles", "etas.txt"),
		          quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t");
		if(verbose) message('Saving SNP matrix (txt)');
		write.table(snps,  file = file.path(saveDir, "txtfiles", "snps.txt"),
		          quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t");
		
	}
	if(savefmat){
		# genepos
		# snpspos = snpsloc$pos;
		
		# Making directory for filematrices and going there
		if(!dir.exists(file.path(saveDir, "filematrices")))
			dir.create(file.path(saveDir, "filematrices"), recursive = TRUE)
		# setwd(file.path(saveDir, "filematrices"))
		
		if(verbose) message('Saving gene/SNP location files (filematrix)');
		fm = fm.create.from.matrix(file.path(saveDir, "filematrices", 'gene_loc'), genepos)
		close(fm);
		fm = fm.create.from.matrix(file.path(saveDir, "filematrices", 'snps_loc'), snpsloc$pos)
		close(fm);

		# Save genes 
		if(verbose) message('Saving expression matrix (filematrix)');
		fm = fm.create.from.matrix(file.path(saveDir, "filematrices", 'gene'), t(genes))
		close(fm)
    
		# Save cvrt 
		if(verbose) message('Saving covariates (filematrix)');
		fm = fm.create.from.matrix(file.path(saveDir, "filematrices", 'cvrt'), t(cvrts))
		close(fm)
		
		if(verbose) message('Saving true effect info (filematrix)');
		fm = fm.create.from.matrix(file.path(saveDir, "filematrices", 'etas'), as.matrix(tbleta))
		close(fm)
	
		# Save SNPs in a filematrix
		if(verbose) message('Saving SNP matrix (filematrix)');
		fm = fm.create(file.path(saveDir, "filematrices", 'snps'), 
							nrow = ncol(snps), ncol = nrow(snps), type = 'integer', size = 1);

		step1 = 100000;
		mm = nrow(snps);
		nsteps = ceiling(mm/step1);
		for( part in 1:nsteps ) { # part = 1
			if (verbose)
				message('Filling SNP filematrix, ', part, 'of', nsteps);
			fr = (part-1)*step1 + 1;
			to = min(part*step1, mm);
			
			fm[,fr:to] = t(snps[fr:to,]);
			gc();
		}
		rm(part, step1, mm, nsteps, fr, to);
		colnames(fm) = rownames(snps);
		rownames(fm) = colnames(snps);
		close(fm);
		
	}
  
	if(returnData) {
		return(list(
			"snps" = snps, 
			"gene" = genes,
			"cvrt" = cvrts,
			"gene_loc" = genepos, 
			"snps_loc" = snpsloc$pos,
			"etas" = tbleta));
	} else {
		return(invisible(NULL));
	}
}
