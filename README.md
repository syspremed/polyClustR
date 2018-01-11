polyClustR
================

Introduction
------------

The polyClustR package provides an R implementation of the pipeline described in ***polyCluster: Defining Communities of Reconciled Cancer Subtypes with Biological and Prognostic Significance***; K. Eason, G. Nyamundanda and A. Sadandandam.

Installation & Usage
--------------------

``` r
devtools::install_github('syspremed/polyClustR')

library(polyClustR)

polyCluster('path/to/expression/file', clusterAlg = c('hc', 'km'), ref = 'test_run')
# Output is written to the current working directory.
```

For full details of the arguments required to run `polyCluster`, see `?polyCluster`.
If you have any problems running polyCluster, please let us know using the "Issues" tab above.
