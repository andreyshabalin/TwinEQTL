\name{TwinMeta_testAll}
\alias{TwinMeta_testAll}
\title{
	Main function for fast eQTL analysis of correlated samples in TwinMeta package.
}
\description{
	\code{TwinMeta_testAll} function tests for association 
	between every gene and every SNP in the input data.
}
\usage{
    TwinMeta_testAll(gene, snps, cvrt, twininfo, pvthreshold)
}
\arguments{
    \item{gene}{
        A matrix with gene expression data, one column per sample and one row per gene. \cr
        The matrix must have both row (gene) and column (sample) names.
    }
    \item{snps}{
        A matrix with genotype data, one column per sample and one row per SNP. \cr
        The matrix must have both row (SNP) and column (sample) names.
    }
    \item{cvrt}{
        A matrix with covariate data, one column per sample and one row per covariate. \cr
        The matrix must have both row (covariate) and column (sample) names. \cr
        Can be \code{NULL} to indicate no covariates (constant still included).
    }
    \item{twininfo}{
        A data frame with the information about all twin pairs, one twin pair per row. Like this: \cr
        \tabular{rrr}{
            \strong{first}    \tab \strong{second}   \tab \strong{type}   \cr
            Sample1  \tab Sample2  \tab MZ      \cr
            Sample8  \tab Sample9  \tab DZ      \cr
            \ldots   \tab \ldots   \tab \ldots \cr
        } 
        
        Sample ID in the first and second columns must match the column names of \code{gene} and \code{snps}.\cr
        The third column can contain only \code{"MZ"} and \code{"DZ"} values.
    }
	\item{pvthreshold}{
		\code{numeric}. Significance threshold. Only associations significant at this level will be returned.
	}
}
\details{
	Note that the columns of \code{gene}, \code{snps}, and \code{cvrt} (if not \code{NULL}) must match.
}
\value{
	Returns data frame with the significant associations (z-scores and p-values).
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}

	Maintainer: Andrey A Shabalin <andrey.shabalin@gmail.com>
}
\references{
	The package website: \url{https://github.com/andreyshabalin/TwinMeta}
}
\seealso{
	The code below is the sample code for TwinMeta analysis on artificial data.	
}
\examples{
library(TwinMeta)

# Number of MZ twin pairs
Nm = 1000

# Number of DZ twin pairs
Nd = 2000

# Number of singleton samples
Ns = 3000

# Number of genes
Ngene = 1000

# Number of SNPs
Nsnps = 1000

# Number of covariates
Ncvrt = 10

# Gerenate artificial data
sim = TwinMeta_simulate( Nm, Nd, Ns, Ngene, Nsnps, Ncvrt)

# Pick a p-value threshold
pvthreshold = 1000 / (Ngene * Nsnps)

# Run eQTL analysis on the data with twins
eqtls = TwinMeta_testAll(
    gene = sim$gene,
    snps = sim$snps,
    cvrt = sim$cvrt,
    twininfo = sim$twininfo,
    pvthreshold = pvthreshold)

# Display the results
head(eqtls)
}
