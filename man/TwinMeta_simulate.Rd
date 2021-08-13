\name{TwinMeta_simulate}
\alias{TwinMeta_simulate}
\title{
	Generate data for test running the package.
}
\description{
	\code{TwinMeta_simulate} function creates artificial data for test running
	the function \code{\link{TwinMeta_testAll}}.
}
\usage{
    TwinMeta_simulate(Nm, Nd, Ns, Ngene, Nsnps, Ncvrt)
}
\arguments{
	\item{Nm}{Number of MZ twin pairs.}
	\item{Nd}{Number of DZ twin pairs.}
	\item{Ns}{Number of singleton samples.}
	\item{Ngene}{Number of genes in gene expression data.}
	\item{Nsnps}{Number of SNPs in genotype data.}
	\item{Ncvrt}{Number of covariates.}
}
\value{
	Returns a list with \code{
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}

	Maintainer: Andrey A Shabalin <andrey.shabalin@gmail.com>
}
\references{
	The package website: \url{https://github.com/andreyshabalin/TwinMeta}
}
\seealso{
	The code below is the sample code for CorrMeta analysis NOT using gene/SNP locations.	
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

show(sim)
}
