estimate_multiSNP = function (SNPs, expr, cvrt, method = "BFGS"){
    
    # Prepping data
    expr = as.vector(expr)
    ncov = nrow(cvrt); 
    nsnp = nrow(SNPs); 
    n = length(expr)
    X = rbind(rep(1, n), SNPs)
    
    if(max(SNPs) > 2)
        SNPs = SNPs / 1000;
    if(max(SNPs) > 2 || min(SNPs) < 0)
        stop("SNP data ill-formatted: ",
            "condition max(SNPs) > 2 || min(SNPs) < 0 is TRUE");
    
    # RSS function
    loglik = function(BetaGam){
        Beta = BetaGam[1:(nsnp + 1)];
        Gam = BetaGam[(nsnp + 2):(nsnp + 1 + ncov)];
        Resid = expr - log(crossprod(X, Beta)) - crossprod(cvrt, Gam);
        return(sum(Resid^2));
    }
    
    # Gradient function
    Dloglik = function(BetaGam){
        Beta = BetaGam[1:(nsnp + 1)];
        Gam = BetaGam[(nsnp + 2):(nsnp + 1 + ncov)];
        Resid = expr - log(crossprod(X, Beta)) - crossprod(cvrt, Gam);
        DBeta = -2 * crossprod(t(X) / as.vector(crossprod(X, Beta)), Resid);
        DGam = -2 * cvrt %*% Resid;
        return(c(DBeta, DGam));
    }
    
    # Finding start value for beta0
    beta0start = exp(mean(expr));
    parStart = c(beta0start, rep(0, nsnp + ncov));
    
    if(method == "BFGS"){
        estimates = optim(parStart, loglik, gr = Dloglik, method = "BFGS")$par;
    } else {
        estimates = optim(parStart, loglik)$par;
    }
    
    # Calculating RSS and sigmahat
    RSSfinal = loglik(estimates);
    DF = n - nsnp - ncov - 1;
    
    return(list("beta0hat" = estimates[1],
                "betahat" = estimates[2:(nsnp + 1)],
                "gamhat" = estimates[(nsnp + 2):(nsnp + ncov)],
                "sigmahat" = RSSfinal / DF,
                "DF" = DF));
    
}

#-------------------------------------------------------------------------------

multisnpACME = function (
        genefm = "gene",
        snpsfm = "snps",
        glocfm = "gene_loc",
        slocfm = "snps_loc",
        cvrtfm = "cvrt",
        acmefm = "ACME",
        workdir = ".",
        genecap = Inf,
        verbose = TRUE){

    ### Orthonormalize covariates
    {
        if(verbose)
            message("Loading and orthonormalizing covariates");
        cvrt = fm.load(file.path(workdir, cvrtfm));
        cvrt_qr = qr.Q(qr(cbind(1, cvrt)));
    } # cvrt, cvrt_qr
    
    ### Gene/SNP locations    
    {
        if(verbose)
            message("Loading gene/SNP locations");
        gene_loc = fm.load(file.path(workdir, glocfm));
        snps_loc = fm.load(file.path(workdir, slocfm));
        stopifnot( !is.unsorted(snps_loc) );
    } # gene_loc, snps_loc

    ### Get matrix sizes and gene names
    {
        if(verbose)
            message("Checking gene/SNP filematrices");
        sfm = fm.open(file.path(workdir, snpsfm), readonly = TRUE);
        gfm = fm.open(file.path(workdir, genefm), readonly = TRUE);
        
        snames = colnames(sfm);
        gnames = colnames(gfm);
        
        stopifnot( ncol(gfm) == nrow(gene_loc) );
        stopifnot( ncol(sfm) == nrow(snps_loc) );
        stopifnot( nrow(gfm) == nrow(cvrt) );
        stopifnot( nrow(sfm) == nrow(cvrt) );
        
        n = nrow(cvrt);
        p = ncol(cvrt);
    } # gfm, sfm, n, p
    
    ### Create output matrix
    {
        if(verbose)
            message("Creating output filematrix");
        fmms = fm.create(
                filenamebase = file.path(workdir, paste0(acmefm, "_multiSNP")),
                nrow = 5,
                ncol = 0);
        rownames(fmms) = c( "geneid", "snp_id", 
                            "beta0", "betas", "forward_adjR2");
    } # fm
    
    ### Doing multi-SNP estimation
    
    afm = fm.open(file.path(workdir, acmefm), readonly = TRUE);
    cur_start = 1;
    count = 1;
    
    while(count <= genecap && cur_start <= ncol(gfm)){
        
        cur_gene = as.numeric(afm[1, cur_start]);
        expr = gfm[, cur_gene];
        
        # Scanning afm for cur_gene
        if(verbose)
            message("Scanning for gene ", count);
        
        # Find all local genes
        tempstart = cur_start;
        snpcount = 0;
        
        repeat{
            gids = as.vector(afm[1, tempstart:(tempstart+9)]);
            matches = sum(cur_gene == gids);
            snpcount = snpcount + matches;
            if(matches != 10){
                cur_end = cur_start + snpcount - 1;
                break;
            }
            tempstart = tempstart + 10;
        }
        
        # Getting cur_gene locations and associated snps
        afm_submat = afm[, cur_start:cur_end];
        gsnps = afm_submat[2, ];
        
        # Initializing inner loop
        last_snp = 0;
        snp_list = integer(0);
        final_AR2 = numeric(0);
        remainingSNPs = 1:ncol(afm_submat);
        SST = afm_submat[7, 1];
        SSEs = afm_submat[6, ];
        adjR2 = 1 - (SSEs / SST) / 
                    ((n - p - length(snp_list) - 2) / (n - p - 1));
        old_AR2 = -Inf
        new_AR2 = max(adjR2);
        
        # Forward selection loop
        while(
                (new_AR2 > old_AR2) && 
                (length(snp_list) < n - p - 2) && 
                (length(remainingSNPs) > 1)){
            
            next_snp = remainingSNPs[which.max(adjR2)];
            snp_list = c(snp_list, next_snp);
            final_AR2 = c(final_AR2, new_AR2);
            remainingSNPs = setdiff(remainingSNPs, next_snp);
            
            if(verbose) 
                message("----Gene ", count, 
                        ", added SNP ", next_snp,
                        " with AdjR2 ", round(max(adjR2), 4));
            
            # Computing SSEs
            real_snp_list = afm_submat[2, snp_list];
            SNPmat = rbind(t(sfm[ , real_snp_list]), 0);
            SSEs = numeric(length(remainingSNPs));
            for(rsi in seq_along(remainingSNPs)){
                
                real_snps = afm_submat[2, c(snp_list, remainingSNPs[rsi])];
                SNPmat[nrow(SNPmat) , ] = 
                    as.vector(sfm[ , afm_submat[2, remainingSNPs[rsi]]]);
                suppressWarnings(
                    {msacme_out = estimate_multiSNP(SNPmat, as.vector(expr), t(cvrt_qr))}
                );
                
                SSEs[rsi] = msacme_out$DF * msacme_out$sigmahat;
            }
            
            adjR2 = 1 - (SSEs / SST) / ((n - p - length(snp_list) - 2) / (n - p - 1));
            next_snp = remainingSNPs[which.max(adjR2)];
            old_AR2 = new_AR2;
            new_AR2 = max(adjR2);
        }
        
        if(length(remainingSNPs) == 1){ 
            # Decide whether or not to add the last SNP
            if(new_AR2 <= old_AR2){
                SNPmat = subset(SNPmat, c(rep(T, nrow(SNPmat) - 1), FALSE));
            } else {
                snp_list = c(snp_list, next_snp);
                final_AR2 = c(final_AR2, new_AR2);
            }
        } else {
            SNPmat = subset(SNPmat, c(rep(T, nrow(SNPmat) - 1), FALSE))
        }
            
        # Computing final estimates and writing the results
        suppressWarnings(
            {final_msacme_out = estimate_multiSNP(SNPmat, as.vector(expr), t(cvrt))}
        );
        real_snps = afm_submat[2, snp_list];
        outmatrix = rbind(
                        rep(cur_gene, length(real_snps)),
                        real_snps,
                        rep(final_msacme_out$beta0hat, length(real_snps)),
                        final_msacme_out$betahat,
                        final_AR2);
        fmms$appendColumns(outmatrix);

        
        # Re-setting loop
        cur_start = cur_end + 1;
        count = count + 1;
    }
    close(fmms);    
}
