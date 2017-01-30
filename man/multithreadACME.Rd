\name{multithreadACME}
\alias{multithreadACME}
\title{
	Fast Parallelized Estimation of ACME Model for All Local Gene-SNP Pairs
}
\description{
This function estimates the ACME model (see the vignette for model details)
for all gene-SNP pairs within pre-defined distance (\code{cisdist}).
The input data must be stored in filematrices 
(see \link[filematrix]{filematrix} package)
and the results are also saved in a filematrix.
This allows the function to perform estimation
using multiple CPU cores in parallel
without having to duplicate the data across all jobs.
}
\usage{
multithreadACME(
    genefm = "gene",
    snpsfm = "snps",
    glocfm = "gene_loc",
    slocfm = "snps_loc",
    cvrtfm = "cvrt",
    acmefm = "ACME",
    cisdist = 1e+06, 
    threads = -1,
    workdir = NULL)
}
\arguments{
  \item{genefm}{
    Name of the filematrix with gene expression data. One column per gene and one row per sample.
}
  \item{snpsfm}{
    Name of the filematrix with SNP data. One column per SNP and one row per sample.
}
  \item{glocfm}{
    Name of the filematix with gene location information. 
    Must contain two columns, first with gene start location and second with the gene end.
    The locations must be stored as numbers, the locations for different chromosomes must differ greatly.
    We suggest encoding (location = 1e9 * chromosome + position_on_chromosome).
    The rows must match the columns of the \code{genefm} filematrix.
}
  \item{slocfm}{
    Name of the filematrix with SNP locations.
    Must have one column and rows matching columns of \code{snpsfm} filematrix.
    See the instructions for \code{glocfm} above.
}
  \item{cvrtfm}{
    Name of the filematirx with covariates.
    Must not include constant (it is added automatically).
    One column per covariate and one row per sample.
}
  \item{acmefm}{
    Name of the filematrix to store the estimates. 
    The filemarix will be created.
    If the filematrix exists, it will be overwritten.
}
  \item{cisdist}{
    The maximum allowed distance between genes and SNPs. 
    Gene-SNP pairs further than \code{cisdist} apart will not be tested.
}
  \item{threads}{
    The number of local jobs (CPU cores) used for calculation.
    If negative, \code{threads} is set to the number of cores of the host machine. 
}
  \item{workdir}{
    Directory where the input filematrices are located.
}
}
\details{
  
}
\value{
  The function creates a filematrix named \code{acmefm} with 10 rows
  and a column for each tested gene-SNP pair.
  The rows contain gene-SNP ids and the estimates by \code{\link{effectSizeEstimation}}:
    \item{geneid}{The gene id - the column number for the gene in the \code{genefm} filematrix.}
    \item{snp_id}{The SNP id - the column number for the SNP in the \code{snpsfm} filematrix.}
	\item{beta0}{The constant parameter in the non-linear model (within the logarithm).}
	\item{beta1}{The effect size parameter in the non-linear model (within the logarithm).}
	\item{nits}{Number of iterations till convergence of the estimation algorithm.}
	\item{SSE}{Sum of squared residuals of the fitted model.}
	\item{SST}{Sum of squared residuals of the model with zero effect.}
	\item{F}{The F test for the significance of the genotype effect.}
	\item{eta}{The effect size parameter for simplified model (beta1/beta0).}
	\item{SE_eta}{Standard error of the eta estimate.}
}
\references{
	The manuscript is available at: \url{http://arxiv.org/abs/1605.08799}
}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\note{
  The rows of \code{genefm}, \code{snpsfm}, and \code{cvrtfm} filematrices must match.
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
tempdirectory = tempdir();
z = create_artificial_data(
    nsample = 100,
    ngene = 500,
    nsnp = 5000,
    ncvrt = 1,
    minMAF = 0.2,
    saveDir = tempdirectory,
    returnData = FALSE,
    savefmat = TRUE,
    savetxt = FALSE,
    verbose = FALSE)

# In this example, we use 2 CPU cores (threads) 
# for testing of all gene-SNP pairs within 100,000 bp.
multithreadACME(
    genefm = "gene",
    snpsfm = "snps",
    glocfm = "gene_loc",
    slocfm = "snps_loc",
    cvrtfm = "cvrt",
    acmefm = "ACME",
    cisdist = 100e+03, 
    threads = 2,
    workdir = paste0(tempdirectory,"/filematrices"),
    verbose = FALSE)

Now the filematrix `ACME` holds estimations for all local gene-SNP pairs.

fm = fm.open(paste0(tempdirectory,"/filematrices/ACME"))
TenResults = fm[,1:10];
rownames(TenResults) = rownames(fm);
close(fm);

show(t(TenResults))
}