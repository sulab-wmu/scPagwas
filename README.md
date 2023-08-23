# About scPagwas

**scPagwas** employs a polygenic regression model to prioritize a set of
trait-relevant genes and uncover trait-relevant cell subpopulations by
incorporating pathway activity transformed single-cell RNA sequencing
(scRNA-seq) data with genome-wide association studies (GWAS) summary
data.

<img src="./man/figures/Figure1.png" width="100%" style="display: block; margin: auto;" />

Please cite this article in press as: Ma et al.,Polygenic regression
uncovers trait-relevant cellular contexts through pathway activation
transformation of single-cell RNA sequencing data,Cell Genomics
(2023),<https://doi.org/10.1016/j.xgen.2023.100383>

Code for reproducing the analysis from the paper is available
[here](https://github.com/dengchunyu/scPagwas_reproduce), or
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8137370.svg)](https://doi.org/10.5281/zenodo.8137370)

For further usage on the scPagwas package, you can visit the
[website](https://dengchunyu.github.io/about/). A vignette for using
also can be accessed using browseVignettes(“scPagwas”)

Some important data can be download from
[here](https://drive.google.com/drive/folders/1z7uQtzjnieJhYLhLmgmIz2jk6a6hes4l?usp=drive_link)

## Installation

You can install the released version of scPagwas from
[github](https://github.com/sulab-wmu/scPagwas) with:

``` r
#install some dependence packages
install.packages("Seurat")
install.packages("ggpubr")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("IRanges")

devtools::install_github("sulab-wmu/scPagwas")
```

In many cases, installing packages using `devtools::install_github` may
fail.

``` r
library(devtools)
install_git("https://github.com/sulab-wmu/scPagwas.git", ref = "main")
```

Or, download the package file from
[here](https://api.github.com/repos/sulab-wmu/scPagwas/tarball/HEAD)
Then install it locally.

``` r
devtools::install_local("sulab-wmu-scPagwas-****.tar.gz")
```

## Usage

quick-start example:

``` r
library(scPagwas)
system.time(
 #1.start to run the wrapper functions for example.
 Pagwas_data<-scPagwas_main(Pagwas = NULL,
                     gwas_data =system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas"), # The GWAS Summary statistics files 
                     Single_data =system.file("extdata", "scRNAexample.rds", package = "scPagwas"),# scRNA-seq data in seruat format with "RNA" assays and normalized.
                     output.prefix="test", # the prefix name for output files
                     output.dirs="scPagwastest_output",# the directory file's name for output
                     block_annotation = block_annotation,# gene position in chromosome is provided by package.
                     assay="RNA", # the assays for scRNA-seq data to use.
                     Pathway_list=Genes_by_pathway_kegg,# pathway list is provided by package, including gene symbols.
                     n.cores=1,
                     iters_singlecell = 10,
                     chrom_ld = chrom_ld,# The LD data is provided by package.
                     singlecell=T, # Whether to run the singlecell process.
                     celltype=T# Whether to run the celltype process.
)
)
```
