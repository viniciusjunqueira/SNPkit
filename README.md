# SNPkit

<!-- badges: start -->
[![R-CMD-check](https://github.com/viniciusjunqueira/SNPkit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/viniciusjunqueira/SNPkit/actions/workflows/R-CMD-check.yaml) [![CRAN status](https://www.r-pkg.org/badges/version/SNPkit)](https://CRAN.R-project.org/package=SNPkit) [![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/SNPkit)](https://CRAN.R-project.org/package=SNPkit) [![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://viniciusjunqueira.github.io/SNPkit/) [![GitHub issues](https://img.shields.io/github/issues/viniciusjunqueira/SNPkit)](https://github.com/viniciusjunqueira/SNPkit/issues)
<!-- badges: end -->

`SNPkit` is an R package designed for manipulation, organization, and analysis of genotypic data, with a strong focus on integration with tools such as **FImpute** and **PLINK**.

It provides robust S4-based data structures for storing genotypes and marker maps, along with functions to combine different genotype panels, summarize data, and prepare files for imputation and selection pipelines.

Key capabilities:

-   Import Illumina `FinalReport.txt` files (any panel density) and merge multiple genotype panels into a single object.
-   Quality control on SNPs and samples (call rate, MAF, HWE, monomorphic, duplicated positions, chromosome filters).
-   Prepare and run **FImpute** imputation and export to **PLINK**.
-   PCA (`runPCA()`) and anticlustering (`runAnticlusteringPCA()`) utilities for exploring structure and building balanced groups (e.g. batch design).

------------------------------------------------------------------------

## 📦 Installation

`SNPkit` depends on [`snpStats`](https://bioconductor.org/packages/snpStats/), which is distributed through **Bioconductor**. Install it first:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("snpStats")
```

Then install the stable release from [CRAN](https://CRAN.R-project.org/package=SNPkit):

``` r
install.packages("SNPkit")
```

Or install the development version (latest features) from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("viniciusjunqueira/SNPkit")
```

### Optional: faster PCA

`runPCA()` and `runAnticlusteringPCA()` can use [`RSpectra`](https://CRAN.R-project.org/package=RSpectra) for a much faster, low-memory truncated PCA on wide genotype data. It is optional — install it to enable the fast path:

``` r
install.packages("RSpectra")
```

------------------------------------------------------------------------

## 📖 Documentation

The full package website with detailed function reference and vignettes is available at:

-   [SNPkit site](https://viniciusjunqueira.github.io/SNPkit/)

Key pages:

-   [Reference index](https://viniciusjunqueira.github.io/SNPkit/reference/index.html)
-   [Vignettes and Tutorials](https://viniciusjunqueira.github.io/SNPkit/articles/)

------------------------------------------------------------------------

## 📄 License

SNPkit is licensed under the [GPL-3](https://www.gnu.org/licenses/gpl-3.0.html) license.
