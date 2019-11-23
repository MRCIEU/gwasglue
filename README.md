# gwasglue

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build status](https://travis-ci.org/MRCIEU/gwasglue.svg?branch=master)](https://travis-ci.org/MRCIEU/gwasglue)
[![Codecov test coverage](https://codecov.io/gh/MRCIEU/gwasglue/branch/master/graph/badge.svg)](https://codecov.io/gh/MRCIEU/gwasglue?branch=master)
<!-- badges: end -->

This R package serves as a conduit between packages that can read or query GWAS summary data, and packages that can analyse GWAS summary data. Here is where it lies in the general ecosystem of GWAS data and analysis:


![schematic](https://drive.google.com/uc?id=15CuzYDlNJeRREQ1ZTF7jLsVGc86OnF7Z)


It currently glues data from the following data sources

#### Data sources
- [ieugwasr](https://github.com/mrcieu/ieugwasr)
- [gwasvcf](https://github.com/mrcieu/gwasvcf)

#### Finemapping
- [finemapr](https://github.com/variani/finemapr)
- [FINEMAP](http://www.christianbenner.com/)
- [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0)
- [CAVIAR](https://github.com/fhormoz/caviar)
- [SuSIE](https://stephenslab.github.io/susie-paper/index.html) - TODO
- [JAM](https://github.com/pjnewcombe/R2BGLiMS) - TODO

#### Colocalisation
- [coloc](https://cloud.r-project.org/web/packages/coloc/index.html)
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
- [gassocplot](https://github.com/jrs95/gassocplot)
- Locus zoom plots e.g. [https://github.com/jrs95/gassocplot] - TODO


## Installation

You can install the development version of gwasglue with:

``` r
devtools::install_github("mrcieu/gwasglue")
```


## Usage

See vignettes etc here: [https://mrcieu.github.io/gwasglue](https://mrcieu.github.io/gwasglue).

## Reference datasets

Example GWAS VCF (GIANT 2010 BMI):

- http://fileserve.mrcieu.ac.uk/vcf/IEU-a-2.vcf.gz
- http://fileserve.mrcieu.ac.uk/vcf/IEU-a-2.vcf.gz.tbi

1kg European reference panel for LD:

- http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz

1kg vcf harmonised against human genome reference:

- http://fileserve.mrcieu.ac.uk/vcf/1kg_v3_nomult.vcf.gz
- http://fileserve.mrcieu.ac.uk/vcf/1kg_v3_nomult.vcf.gz.tbi

## Contributing to the resource

For any `<analysis>` package we create a new file called `R/<analysis>.r` which contains two functions:

- `gwasvcf_to_<analysis>`
- `ieugwasr_to_<analysis>`

For an example, see the `R/TwoSampleMR.r` file, which contains the functions `gwasvcf_to_TwoSampleMR` and `ieugwasr_to_TwoSampleMR`.

