polyClustR
================

Introduction
------------

Installation & Usage
--------------------

``` r
devtools::install_github('syspremed/polyClustR')
library(polyClustR)

polyCluster('path/to/expression/file', clusterAlg = c('hc', 'km'), ref = 'test_run')
```
