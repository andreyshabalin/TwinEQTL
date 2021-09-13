# TwinEQTL: Ultra Fast and Powerful Association Analysis for eQTL and GWAS in Twin Studies

## Introduction

TwinEQTL is a computationally efficient alternative to a linear mixed-effects model (LMM)
for twin genome-wide association study (GWAS) or expression quantitative trait loci (eQTL) analyses.
Instead of analyzing all twin samples together with LMM,
TwinEQTL first randomly splits twin samples into two independent
groups on which multiple linear regression analysis is performed separately,
followed by an appropriate meta-analysis to combine the two non-independent test results.
Our approaches provide a huge leap in terms of computing performance for GWAS data with twin pairs. 

## Key features
1. Similar to meta-analysis, only summarized SNP level test statistics are necessary
2. Fast alternative to linear mixed effect model with no inflation of type I error and negligible power loss
3. Fast standard GWAS analysis for twin or correlated subjects
4. Fast expression quantitative trait loci (eQTL) analysis for twin or correlated subjects
5. Implemented as an easy-to-use R package similar to MatrixEQTL

## Install GitHub Version
To install `TwinEQTL` directly from GitHub, run

```r
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("andreyshabalin/TwinEQTL@main")
```

The package includes reference manual, sample data and a Vignette.

## Basic Usage

```r
library(TwinEQTL)

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
sim = TwinEQTL_simulate(Nm, Nd, Ns, Ngene, Nsnps, Ncvrt)

# Pick a p-value threshold
pvthreshold = 1000 / (Ngene * Nsnps)

# Run eQTL analysis on the data with twins
eqtls = TwinEQTL_testAll(
    gene = sim$gene,
    snps = sim$snps,
    cvrt = sim$cvrt,
    twininfo = sim$twininfo,
    pvthreshold = pvthreshold)

# Display the results
head(eqtls)
```

## Citation
K Xia, AA Shabalin, ..., F Zou. Ultra Fast and Powerful Association
Analysis for eQTL and GWAS in Twin Studies. (Submitted)

## Contact
Kai Xia: kxia@med.unc.edu

Andrey A Shabalin: andrey.shabalin@gmail.com

Fei Zou: fzou@bios.unc.edu

