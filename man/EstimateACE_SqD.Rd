\name{EstimateACE_SqD}
\alias{EstimateACE_SqD}
\title{
	Estimate ACE model for each gene
}
\description{
    Estimate ACE model using SqD method.
}
\usage{
    EstimateACE_SqD(gene, cvrt, twininfo)
}
\arguments{
    \item{gene}{
        A matrix with gene expression data, one column per sample and one row per gene. \cr
    }
    \item{cvrt}{
        A matrix with covariate data, one column per sample and one row per covariate. \cr
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
}
\details{
	Note that the columns of \code{gene}, and \code{cvrt} (if not \code{NULL}) must match.
}
\value{
	Returns a matrix with \code{nrow(gene)} rowa and 3 columds.\cr
	Each row contains estimates for \eqn{(\sigma^2_A, \sigma^2_C, \sigma^2_E)} parameters for the corresponding gene.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}

	Maintainer: Andrey A Shabalin <andrey.shabalin@gmail.com>
}
\references{
	The package website: \url{https://github.com/andreyshabalin/TwinMeta}
}
\seealso{
	The code below estimates ACE model parameters for 1000 genes and 
	inspects of the distribution of the estimates.
}
\examples{
library(TwinMeta)

## Model parameters

# Number of MZ twin pairs
Nm = 1000

# Number of DZ twin pairs
Nd = 2000

# Number of singleton samples
Ns = 3000

# Number of genes
Ngene = 1000

# Number of SNPs (minimal, as we don't need genotypes)
Nsnps = 1

# Number of covariates
Ncvrt = 1

# ACE model parameters, same for all genes
ACEparam = c(3,4,5);

# Simulate the data
sim = TwinMeta_simulate(Nm, Nd, Ns, Ngene, Nsnps, Ncvrt, ACEparam = ACEparam);

# Perform ACE model estimation
ace = EstimateACE_SqD(gene = sim$gene, cvrt = sim$cvrt, twininfo = sim$twininfo);

# Scale estimates to the same sum as the true value.
ace = ace / rowSums(ace) * sum(ACEparam);

# Look at the distribution of A estimates
hist(ace[,1], 100, col = 'lightblue', main = 'Estimates of A');
abline(v = ACEparam[1], lwd = 3, col = 'red');

# Look at the distribution of C estimates
hist(ace[,2], 100, col = 'lightblue', main = 'Estimates of C');
abline(v = ACEparam[2], lwd = 3, col = 'red');

# Look at the distribution of E estimates
hist(ace[,3], 100, col = 'lightblue', main = 'Estimates of E');
abline(v = ACEparam[3], lwd = 3, col = 'red');
}
