
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SimBu

<!-- badges: start -->

[![R-CMD-check](https://github.com/omnideconv/simulator/workflows/R-CMD-check/badge.svg)](https://github.com/omnideconv/simulator/actions)
[![docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://alex-d13.github.io/index.html)
<!-- badges: end -->

The goal of SimBu is to simulate pseudo-bulk RNAseq datasets with
variable cell-type fractions baed on public or private single-cell
RNAseq datasets.

## Installation

To install the developmental version of the package, run:

``` r
install.packages("devtools")
devtools::install_github("omnideconv/SimBu") 
```

To install from Bioconductor:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SimBu")
```

## Usage

Create a dataset-object with local data and simulate a pseudo-bulk
dataset

``` r
library(SimBu)
# use local data to build dataset
dataset <- dataset(annotation = annotation_dataframe, count_matrix = expression_matrix, name = "test_dataset")
simulation <- simulate_bulk(data = dataset, scenario = "random", scaling_factor = "NONE")
```

For more detailed usage instructions, check out [Getting
started](http://omnideconv.org/SimBu/articles/simulator_documentation.html).

You can find more information on the simulation setup and mRNA bias in
our [SimBu preprint](https://doi.org/10.1101/2022.05.06.490889).
