
<!-- README.md is generated from README.Rmd. Please edit that file -->
calcWCD
=======

calcWCD calculates weighted copublication distance between a gene/protein and a disease/term in PubMed literature.

Installation
------------

You can install the released version of calcWCD from this repo with:

``` r
devtools::install_github("ed-lau/calcWCD")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
require(calcWCD)
#> Loading required package: calcWCD
wcd(pmids = test_pmids, annotation = test_annotation)
#> Loading required package: dplyr
#> Warning: package 'dplyr' was built under R version 3.5.1
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Joining, by = "GeneID"
#> # A tibble: 695 x 8
#>    GeneID Term_count Total_count WCD_n WCD_d   WCD     Z      P
#>     <int>      <dbl>       <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
#>  1   5972     241.        508.    1.15  2.56 0.448 -2.72 0.0032
#>  2   1636     169.        498.    1.30  2.57 0.507 -2.33 0.0098
#>  3    183     154.        445.    1.34  2.62 0.513 -2.29 0.011 
#>  4   8991      82.3       153.    1.62  3.08 0.524 -2.22 0.013 
#>  5  59272      22.6        40.6   2.18  3.66 0.595 -1.74 0.041 
#>  6   1628      40.1       114.    1.93  3.21 0.601 -1.7  0.044 
#>  7    659      18.8        39.3   2.26  3.67 0.615 -1.61 0.054 
#>  8    185      41.3       146.    1.92  3.10 0.617 -1.59 0.056 
#>  9  65266       4.77        4.77  2.85  4.59 0.622 -1.56 0.059 
#> 10   1906      55.2       249.    1.79  2.87 0.623 -1.55 0.06  
#> # ... with 685 more rows
```

Distributin of WCD values in test files

    #> Joining, by = "GeneID"

<img src="man/figures/README-pressure-1.png" width="100%" />
