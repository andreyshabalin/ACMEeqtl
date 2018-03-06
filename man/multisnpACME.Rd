\name{multisnpACME}
\alias{multisnpACME}
\title{
    Estimation of Multi-SNP ACME Model for Full-Tissue Genome and All Local SNPs
}
\description{
    This function estimates gene-wise multi-SNP ACME models. 
    It requires the output of \code{\link{multithreadACME}} 
    to know all the local SNPs for each gene. 
    It then performs forward step-wise selection of the local SNPs, 
    based on the adjusted-R-squared at each step.
    The arguments closely mirror those of \code{\link{multithreadACME}}
    their values must correspond to a set of output files from that function 
    (as well as the input files which originally produced the output).
    It saves the data in filematrix format, 
    similar to the output of \code{\link{multithreadACME}}.
    
    Note that each multi-SNP model will contain at least one SNP, 
    even if that initial SNP was not significant under the single-SNP models.
    This initial SNP will be the one with the highest adjusted-R-squared 
    value among the single-SNP models. However, after the initial SNP, 
    further SNPs are added only if the combined model's adjusted-R-squared is 
    greater than that from the previous combined model.
}
\usage{
multisnpACME(
    genefm = "gene",
    snpsfm = "snps",
    glocfm = "gene_loc",
    slocfm = "snps_loc",
    cvrtfm = "cvrt",
    acmefm = "ACME",
    workdir = ".",
    genecap = Inf,
    verbose = TRUE)
}
\arguments{
    \item{genefm}{
        Name of the filematrix with gene expression data.
        One column per gene and one row per sample.
    }
    \item{snpsfm}{
        Name of the filematrix with SNP data.
        One column per SNP and one row per sample.
    }
    \item{glocfm}{
        Name of the filematix with gene location information. 
        Must contain two columns, 
        first with gene start location and second with the gene end.
        The locations must be stored as numbers, 
        the locations for different chromosomes must differ greatly.
        We suggest encoding 
        (location = 1e9 * chromosome + position_on_chromosome).
        The rows must match the columns of the \code{genefm} filematrix.
    }
    \item{slocfm}{
        Name of the filematrix with SNP locations.
        Must have one column and rows matching columns of 
        \code{snpsfm} filematrix.
        See the instructions for \code{glocfm} above.
    }
    \item{cvrtfm}{
        Name of the filematirx with covariates.
        Must not include constant (it is added automatically).
        One column per covariate and one row per sample.
    }
    \item{acmefm}{
        Name of the filematrix to in which the ACME estimates are stored.
        A new file matrix with the name \code{paste0(acmefm, "_multiSNP")} 
        will be created to store the multi-SNP ACME estimates.
        If the filematrix exists, it will be overwritten.
    }
    \item{workdir}{
        Directory where the input filematrices are located.
    }
    \item{genecap}{
        Number of genes to estimate multi-SNP model for.
    }
    \item{verbose}{
        Set to \code{TRUE} to indicate progress.
    }
}
\value{
    The function creates a filematrix named 
    \code{paste0(acmefm, "_multiSNP")} with 4 rows
    and a column for a SNP when it is included in a mult-SNP model. 
    If the SNP is included in more than one multi-SNP model,
    it will appear multiple times in the matrix 
    (but with different beta estimates, corresponding to the paritular models). 
    The rows contain gene-SNP ids, 
    step-wise adjusted-R-squared statistics, 
    and beta estimates:
    \item{geneid}{
        The gene id - the column number for the gene 
        in the \code{genefm} filematrix.
    }
    \item{snp_id}{
        The SNP id - the column number for the SNP 
        in the \code{snpsfm} filematrix.
    }
    \item{beta0}{
        The beta0 estimate in the full model.
    }
    \item{beta}{
        The beta estimate for the SNP in the full model 
        (after all chosen SNPs have been added).
    }
    \item{forward_adjR2}{
        The step-wise adjusted-R-squared, 
        computed for the full model when the SNP was added.
    }
}
\references{
    The manuscript is available at:
    \url{http://onlinelibrary.wiley.com/doi/10.1111/biom.12810/full}
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}, John Palowitch
}
\note{
  The rows of \code{genefm}, \code{snpsfm}, and \code{cvrtfm} 
  filematrices must match.
  The SNPs must have increasing locations.
}
\seealso{
    For package overview and code examples see the package vignette via:\cr
    \code{browseVignettes("ACMEeqtl")}\cr
    or\cr
    \code{RShowDoc("doc/ACMEeqtl.html", "html", "ACMEeqtl")}
}
\examples{
# First we generate a eQTL dataset in filematrix format 
tempdirectory = tempdir()
z = create_artificial_data(
    nsample = 50,
    ngene = 11,
    nsnp = 51,
    ncvrt = 1,
    minMAF = 0.2,
    saveDir = tempdirectory,
    returnData = FALSE,
    savefmat = TRUE,
    savetxt = FALSE,
    verbose = FALSE)

# Then we run multithreadACME to obtain single-SNP estimates.
# In this example, we use 2 CPU cores (threads) 
# for testing of all gene-SNP pairs within 100,000 bp.
multithreadACME(
    genefm = "gene",
    snpsfm = "snps",
    glocfm = "gene_loc",
    slocfm = "snps_loc",
    cvrtfm = "cvrt",
    acmefm = "ACME",
    cisdist = 10e+06,
    threads = 1, # Use more for faster run
    workdir = file.path(tempdirectory, "filematrices"),
    verbose = FALSE)

# Now the filematrix `ACME` holds estimations for all local gene-SNP pairs.

fm = fm.open(file.path(tempdirectory, "filematrices", "ACME"))
TenResults = fm[,1:10]
rownames(TenResults) = rownames(fm)
close(fm)

show(t(TenResults))

# Now we can estimate multi-SNP ACME models for each gene:
multisnpACME(
    genefm = "gene",
    snpsfm = "snps",
    glocfm = "gene_loc",
    slocfm = "snps_loc",
    cvrtfm = "cvrt",
    acmefm = "ACME",
    workdir = file.path(tempdirectory, "filematrices"),
    genecap = Inf,
    verbose = TRUE)
}