\name{msACME}
\alias{msACME}
\title{
	Estimation of Multi-SNP ACME Model for Full-Tissue Genome and All Local SNPs
}
\description{
This function performs stepwise selection of a multi-SNP ACME model
for all local SNPs, as given by the output of \code{\link[ACMEeqtl]{estimate_multiSNP}}. 
The arguments closely mirror those of \code{\link[ACMEeqtl]{estimate_multiSNP}}, and must
correspond to a certain output set from that function.
The input data must be stored in filematrices 
(see \link[filematrix]{filematrix} package)
and the results are also saved in a filematrix.
This allows the function to perform estimation
using multiple CPU cores in parallel
without having to duplicate the data across all jobs.
}
\usage{
msACME(
	genefm = 'gene',
	snpsfm = 'snps',
	glocfm = 'gene_loc',
	slocfm = 'snps_loc',
	cvrtfm = 'cvrt',
	acmefm = 'ACME',
	workdir = paste0(tempdirectory, "/filematrices"),
	genecap = Inf,
	verbose = TRUE)
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
    Name of the filematrix to in which the ACME estimates are stored.
    A new file matrix with the name \code{paste0(acmefm, 'multiSNP_')} will be created to store the multi-SNP ACME estimates.
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
  The function creates a filematrix named \code{paste0(acmefm, 'multiSNP_')} with 4 rows
  and a column for SNP when it is included in a mult-SNP model. If the SNP is included in more than one
  multi-SNP model, it will appear that many times. in the matrix. 
  The rows contain gene-SNP ids, step-wise adjusted R^2 statistics, and beta estimates:
    \item{geneid}{The gene id - the column number for the gene in the \code{genefm} filematrix.}
    \item{snp_id}{The SNP id - the column number for the SNP in the \code{snpsfm} filematrix.}
	\item{eta}{The beta estimate in the full model (after all SNPs have been added).}
	\item{forward_adjR2}{The step-wise adjusted R^2, computed for the full model when the SNP was added.}
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
    cisdist = 100e+03, 
    threads = 2,
    workdir = paste0(tempdirectory,"/filematrices"),
    verbose = FALSE)

# Now the filematrix `ACME` holds estimations for all local gene-SNP pairs.

fm = fm.open(paste0(tempdirectory,"/filematrices/ACME"))
TenResults = fm[,1:10];
rownames(TenResults) = rownames(fm);
close(fm);

show(t(TenResults))

Now we can estimate multi-SNP ACME models for each gene:
msACME(
	genefm = 'gene',
	snpsfm = 'snps',
	glocfm = 'gene_loc',
	slocfm = 'snps_loc',
	cvrtfm = 'cvrt',
	acmefm = 'ACME',
	workdir = paste0(tempdirectory, "/filematrices"),
	genecap = Inf,
	verbose = TRUE)
}