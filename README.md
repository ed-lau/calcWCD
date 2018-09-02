
<!-- README.md is generated from README.Rmd. Please edit that file -->
calcWCD
=======

calcWCD calculates weighted copublication distance between a gene/protein and a disease/term in PubMed literature.

Installation
------------

The released version of calcWCD can be installed from this repo directly:

``` r
devtools::install_github("ed-lau/calcWCD")
```

Example
-------

A basic example using test PMID and Annotation data included in the package:

``` r
require(calcWCD)
w <- wcd(pmids = test_pmids, annot = test_annotation)

w
#> # A tibble: 695 x 8
#>    GeneID Term_count Total_count WCD_n WCD_d   WCD     Z       P
#>     <int>      <dbl>       <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
#>  1   5972     241.        508.    1.15  2.56 0.448 -2.72 0.00323
#>  2   1636     169.        498.    1.30  2.57 0.507 -2.33 0.00984
#>  3    183     154.        445.    1.34  2.62 0.513 -2.29 0.0110 
#>  4   8991      82.3       153.    1.62  3.08 0.524 -2.22 0.0133 
#>  5  59272      22.6        40.6   2.18  3.66 0.595 -1.74 0.0408 
#>  6   1628      40.1       114.    1.93  3.21 0.601 -1.70 0.0443 
#>  7    659      18.8        39.3   2.26  3.67 0.615 -1.61 0.0538 
#>  8    185      41.3       146.    1.92  3.10 0.617 -1.59 0.0555 
#>  9  65266       4.77        4.77  2.85  4.59 0.622 -1.56 0.0591 
#> 10   1906      55.2       249.    1.79  2.87 0.623 -1.55 0.0604 
#> # ... with 685 more rows
```

Distribution of WCD values in the test files:

<img src="man/figures/README-pressure-1.png" width="100%" />

Reference
---------

Edward Lau, Vidya Venkatraman, Cody T Thomas, Jennifer E Van Eyk, Maggie P. Y. Lam. Identifying high-priority proteins across the human diseasome using semantic similarity. bioRxiv 309203; doi: <https://doi.org/10.1101/309203>
