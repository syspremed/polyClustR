polyClustR
================

Introduction
------------

The polyClustR package provides an R implementation of the pipeline described in *polyCluster: Defining Communities of Reconciled Cancer Subtypes with Biological and Prognostic Significance*, K Eason, G Nyamundanda, A Sadandandam.

Installation & Usage
--------------------

``` r
devtools::install_github('syspremed/polyClustR')
library(polyClustR)

polyCluster('path/to/expression/file', clusterAlg = c('hc', 'km'), ref = 'test_run')
```
