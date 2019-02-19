# ACMEeqtl: Estimation of Interpretable eQTL Effect Sizes Using a Log of Linear Model

We use a non-linear model, termed ACME,
that reflects a parsimonious biological model for
allelic contributions of cis-acting eQTLs.
With non-linear least-squares algorithm we
estimate maximum likelihood parameters. The ACME model
provides interpretable effect size estimates and
p-values with well controlled Type-I error.
Includes both R and (much faster) C implementations.
For more details see
[Palowitch *et al.* (2017)
](http://onlinelibrary.wiley.com/doi/10.1111/biom.12810/abstract).

## Installation

### Install CRAN Version

To install the
[CRAN version](https://CRAN.R-project.org/package=ACMEeqtl)
of `ACMEeqtl`, run

```
install.packages("ACMEeqtl")
```

### Install GitHub Version

To install `ACMEeqtl` directly from GitHub, run

```
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("andreyshabalin/ACMEeqtl")
```
