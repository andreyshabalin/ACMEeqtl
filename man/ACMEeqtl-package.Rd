\name{ACMEeqtl-package}
\alias{ACMEeqtl-package}
\docType{package}
\title{
    Estimation of Interpretable eQTL Effect Sizes Using a Log of Linear Model
}
\description{
We use a non-linear model, termed ACME, 
that reflects a parsimonious biological model for 
allelic contributions of cis-acting eQTLs.
With non-linear least-squares algorithm we 
estimate maximum likelihood paramters. The ACME model
provides interpretable effect size estimates and
p-values with well controlled Type-I error.	
Includes both R and (faster) C implementations.
}
\details{
	\tabular{ll}{
		Package: \tab ACMEeqtl \cr
		Type: \tab Package \cr
		License: \tab LGPL-3 \cr
	}
}
\author{
	Andrey A Shabalin \email{andrey.shabalin@gmail.com}, John Palowitch
}
\references{
	The manuscript is available at: \url{http://arxiv.org/abs/1605.08799}
}
\seealso{
    For package overview and code examples see the package vignette via:\cr
    \code{browseVignettes("ACMEeqtl")}\cr
    or\cr
	\code{RShowDoc("doc/ACMEeqtl.html", "html", "ACMEeqtl")}
	
	For fast testing of all local gene-SNP pairs (local eQTL) see \code{\link{multithreadACME}}.
}
