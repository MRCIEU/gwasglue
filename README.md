# gwasglue

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build status](https://travis-ci.org/MRCIEU/gwasglue.svg?branch=master)](https://travis-ci.org/MRCIEU/gwasglue)
[![Codecov test coverage](https://codecov.io/gh/MRCIEU/gwasglue/branch/master/graph/badge.svg)](https://codecov.io/gh/MRCIEU/gwasglue?branch=master)
<!-- badges: end -->

This R package serves as a conduit between packages that can read or query GWAS summary data, and packages that can analyse GWAS summary data. Here is where it lies in the general ecosystem of GWAS data and analysis:


![schematic](https://drive.google.com/uc?id=1KgjRmqrKA2McoXZrll7GEquenqRMUB5_)


It currently glues data from the following data sources

#### Data sources
- [ieugwasr](https://github.com/mrcieu/ieugwasr)
- [gwasvcf](https://github.com/mrcieu/gwasvcf)

#### Finemapping
- [finemapr](https://github.com/variani/finemapr) 
- [FINEMAP](http://www.christianbenner.com/) - TODO
- [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0) - TODO
- [CAVIAR](https://github.com/fhormoz/caviar) - TODO
- [SuSIE](https://stephenslab.github.io/susie-paper/index.html) - TODO
- [JAM](https://github.com/pjnewcombe/R2BGLiMS) - TODO

#### Colocalisation
- [coloc](https://cloud.r-project.org/web/packages/coloc/index.html) - TODO
- [HEIDI](http://cnsgenomics.com/software/gsmr/) - TODO
- [eCAVIAR](https://github.com/fhormoz/caviar) - TODO
- [S-Predixcan](https://github.com/hakyimlab/MetaXcan) - TODO

#### Mendelian randomization
- [TwoSampleMR](https://github.com/mrcieu/TwoSampleMR)
- [MendelianRandomization](https://cran.r-project.org/web/packages/MendelianRandomization/index.html) - port from TwoSampleMR
- [RadialMR](https://github.com/WSpiller/RadialMR) - port from TwoSampleMR
- [MRPRESSO](https://github.com/rondolab/MR-PRESSO) - port from TwoSampleMR
- [GSMR](http://cnsgenomics.com/software/gsmr/) - TODO
- [MRMix](https://github.com/gqi/MRMix) - TODO

#### Visualisation
- Locus zoom plots e.g. [https://github.com/jrs95/gassocplot] - TODO


## Installation

You can install the development version of gwasglue with:

``` r
devtools::install_github("mrcieu/gwasglue")
```


## Usage

See vignettes etc here: [https://mrcieu.github.io/gwasglue](https://mrcieu.github.io/gwasglue).
