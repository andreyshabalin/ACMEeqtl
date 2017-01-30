process_gene_block = function(param, genefm, snpsfm, cvrt_qr, acmefm) {
	
	library(ACMEeqtl)
	library(filematrix)
	
	fm = fm.open(genefm, readonly = TRUE);
	genemat = fm[, param$geneset];
	close(fm);

	fm = fm.open(snpsfm, readonly = TRUE);
	snpsmat = fm[,min(param$snps_l):max(param$snps_r)];
	if(fm$size==2)
		snpsmat = snpsmat / 1000;
	close(fm)
	
	firstSNP = min(param$snps_l);
	
	snpsmat_fr = param$snps_l - firstSNP;
	snpsmat_to = param$snps_r - firstSNP;
	
	noutput = sum(snpsmat_to-snpsmat_fr+1L);
	outmatrix = double(10*noutput);
	dim(outmatrix) = c(10,noutput);

	{
		tic = proc.time();
		rez = .Call("effectSizeManyC", genemat, snpsmat, snpsmat_fr, snpsmat_to, cvrt_qr, outmatrix, PACKAGE = "ACMEeqtl");
		toc = proc.time();
		# show(toc-tic);
		message('Finished in ', round((toc-tic)[3],3), ' seconds');
	}
	
	outmatrix[1,] = outmatrix[1,] + param$geneset[1];
	outmatrix[2,] = outmatrix[2,] + firstSNP;
	
	fm = fm.open(acmefm);
	fm$writeCols(start = param$offset+1, outmatrix);
	close(fm);
	
	# cat(format(Sys.time(), "%a %b %d %X %Y"), as.character(param$offset), 
	# 	 sum(snpsmat_to-snpsmat_fr+1), '+\n', file = paste0(dir,'/logs/log',Sys.getpid(),'.txt'), 
	# 	 append = TRUE);
	
	return(c(param$offset, param$offset + ncol(outmatrix)));
}

multithreadACME = function(
		genefm = 'gene',
		snpsfm = 'snps',
		glocfm = 'gene_loc',
		slocfm = 'snps_loc',
		cvrtfm = 'cvrt',
		acmefm = 'ACME',
		cisdist = 1e6,
		threads = -1){	

	if(threads <= 0)
		threads = detectCores();
	# dir.create(paste0(dir,'/logs'), showWarnings = FALSE, recursive = TRUE);
	
	### Orthonormalize covariates
	{
		message('Loading and orthonormalizing covariates');
		cvrt = fm.load(cvrtfm);
		cvrt_qr = qr.Q(qr(cbind(1, cvrt)));
	} # cvrt, cvrt_qr

	### Gene/SNP locations	
	{
		message('Loading gene/SNP locations');
		gene_loc = fm.load(glocfm);
		snps_loc = fm.load(slocfm);
	} # gene_loc, snps_loc

	### Get matrix sizes
	{
		message('Checking gene/SNP filematrices');
		gfm = fm.open(genefm, readonly = TRUE);
		sfm = fm.open(snpsfm, readonly = TRUE);
		
		stopifnot( ncol(gfm) == nrow(gene_loc) );
		stopifnot( ncol(sfm) == nrow(snps_loc) );
		stopifnot( nrow(gfm) == nrow(cvrt) );
		stopifnot( nrow(sfm) == nrow(cvrt) );
		
		close(gfm);
		close(sfm);
	}
		
	### Count number of SNPs for each gene
	{
		ind1 = findInterval(gene_loc, snps_loc + cisdist + 1) + 1L;
		ind2 = findInterval(gene_loc, snps_loc - cisdist);
		# total_pairs = sum(ind2-ind1+1);
		cum_pairs = c(0,cumsum(ind2-ind1+1));
		total_pairs = tail(cum_pairs, 1);
		message('Total ', total_pairs, ' local gene-SNP pairs')
	} # ind1, ind2, total_pairs, cum_pairs
	
	# genes_without_snps = ind2<ind1;
	
	### blocks of genes
	
	# from (last snp of the current) - (first of the next snp range)
	# gene_block_starts = seq_len( length(gene_loc)+1);

	# process genes by 100
	# gene_block_starts = c(seq(1, nrow(gene_loc), 100), nrow(gene_loc)+1L)
	
	# Split genes into groups with non-overlapping nearby SNPs
	gene_block_starts = which( c(ind1,1e100) - c(-1e9,ind2) > 0 );
	gb_npairs = diff(cum_pairs[gene_block_starts]);
	# gb_npairs = cum_pairs[gene_block_starts[-1]] - cum_pairs[gene_block_starts[-length(gene_block_starts)]];
	
	maxpairs = 3e5;
	overpairs = which(gb_npairs > 5e5);
	if( length(overpairs) > 0 ) {
		extralist = vector('list', length(overpairs)); 
		for( j in seq_along(overpairs) ) { # j = 1
			i = overpairs[j];
			# cat(j, i,'\n');
			set = gene_block_starts[i]:(gene_block_starts[i+1]);
			cumnum = cum_pairs[set] - cum_pairs[set[1]];
			# cumnum %/% maxpairs
			extralist[[j]] = set[ diff(cumnum %/% maxpairs)>0L ];
		}
		gene_block_starts = sort( c(gene_block_starts, unlist(extralist)));
		gb_npairs = diff(cum_pairs[gene_block_starts]);
		# gb_npairs = cum_pairs[gene_block_starts[-1]] - cum_pairs[gene_block_starts[-length(gene_block_starts)]];
	}

	# max(gb_npairs)
	# sum(gb_npairs>0)
	gb_offsets = cum_pairs[gene_block_starts]
	
	message('Task split into ',length(gene_block_starts),' parts')
	
	# ord = sort.list(gb_npairs, decreasing = TRUE);
	# ord = ord[gb_npairs[ord]>0];

	# ord = which(gb_npairs>0);
	# ord = ord[ order(pmin(gb_npairs[ord],maxpairs*0.8), decreasing = TRUE) ];

	# Create output matrix
	{
		message('Creating output filematrix')
		fm = fm.create(acmefm, nrow = 10, ncol = total_pairs);
		rownames(fm) = c("geneid","snp_id","beta0", "beta1", "nits", "SSE", "SST", "F","eta","SE")
		close(fm)
	}
	
	paramlist = vector('list', length(gb_npairs));
	for( i in seq_along(gb_npairs)) { # i=1
		# i = ord[j];
		geneset = gene_block_starts[i]:(gene_block_starts[i+1]-1L);
		paramlist[[i]] = list(
			geneset = geneset,
			snps_l = ind1[geneset],
			snps_r = ind2[geneset],
			offset = gb_offsets[i]
		);
	}
	
	# param = paramlist[[1]]
	# param = tail(paramlist,1)[[1]]
	
	threads = min(threads, length(gb_npairs));
	
	if( threads > 1 ) {
		message('Starting ACME eQTL analysis in ', threads, ' parallel jobs');
		tic = proc.time();
		cl = makeCluster(threads);
		# clusterExport(cl, varlist = c("genefm", "snpsfm", "cvrt_qr","acmefm"));
		# nmslist = clusterSplit(cl, nms);
		z = clusterApplyLB(cl, paramlist, process_gene_block, 
								 genefm = genefm,
								 snpsfm = snpsfm,
								 cvrt_qr = cvrt_qr,
								 acmefm = acmefm);
		stopCluster(cl)
		toc = proc.time();
		message('Finished in ', round((toc-tic)[3],3), ' seconds');
		# show(toc-tic);
	} else {
		for( i in seq_along(paramlist) ) { # i = 1
			message('Processing part ', i, ' of ', length(paramlist));
			process_gene_block(param = paramlist[[i]], genefm, snpsfm, cvrt_qr, acmefm);
		}
	}
	return(invisible(NULL));
}

