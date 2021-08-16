\name{TwinMeta_simulate}
\alias{TwinMeta_simulate}
\title{
	Generate data for test running the package.
}
\description{
	\code{TwinMeta_simulate} function creates artificial data set for test running
	the main package function \code{\link{TwinMeta_testAll}}.
}
\usage{
    TwinMeta_simulate(Nm, Nd, Ns, Ngene, Nsnps, Ncvrt, ACEparam = NULL, MAF = NULL)
}
\arguments{
	\item{Nm}{Number of MZ twin pairs.}
	\item{Nd}{Number of DZ twin pairs.}
	\item{Ns}{Number of singleton samples.}
	\item{Ngene}{Number of genes in gene expression data.}
	\item{Nsnps}{Number of SNPs in genotype data.}
	\item{Ncvrt}{Number of covariates.}
	\item{ACEparam}{
	    ACE model parameters. \cr
	    Must be either length 3 vector with \eqn{(\sigma^2_A, \sigma^2_C, \sigma^2_E)} or an \code{Nsnps} by 3 matrix of such rows for each SNP.
	}
	\item{MAF}{
	    Minor allele frequency for the generated SNPs. \cr
	    Must be either a single value or a vector of length \code{Nsnps}. \cr
	    The values must be between 0 and 1.
	}
}
\value{
	Returns a list with generated gene expression (\code{gene}), genotype (\code{snps}), covariates (\code{cvrt}),
	and data frame \code{twininfo} with information about the twin pairs in the data.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}

	Maintainer: Andrey A Shabalin <andrey.shabalin@gmail.com>
}
\references{
	The package website: \url{https://github.com/andreyshabalin/TwinMeta}
}
\examples{

library(TwinMeta)

Nm = 1
Nd = 2
Ns = 3

Ngene = 10
Nsnps = 10
Ncvrt = 1

sim = TwinMeta_simulate(Nm, Nd, Ns, Ngene, Nsnps, Ncvrt)

# Show gene expression data
head(sim$gene)

# Show genotypes
head(sim$snps)

# Show covariates
head(sim$snps)

# Show twininfo data frame
head(sim$twininfo)

}
