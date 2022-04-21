
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SimBu

<!-- badges: start -->

[![R-CMD-check](https://github.com/omnideconv/simulator/workflows/R-CMD-check/badge.svg)](https://github.com/omnideconv/simulator/actions)
<!-- badges: end -->

The goal of SimBu is to simulate pseudo-bulk RNAseq datasets from public
or private single-cell RNAseq datasets.

## Installation

To install this package, use the following command:

``` r
install.packages("devtools")
devtools::install_github("omnideconv/SimBu") 
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
