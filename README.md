SMARTS (Statistical Methods for Assessement of real-world Treatment
Switching)
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

SMARTS is a tool for assigning pseudo-switching time to patients who
initiate a treatment and stay on the same treatment (continuers), in
comparison to patients who initiate a treatment and switch to another
treatment (switchers). The assigned pseudo-switching time is later used
in analysis to reduce bias caused by time varying confounding.

## Install using devtools

``` r
library(devtools)
install_github("chuangao/SMARTS")
```

## Install from source

``` bash
git clone https://github.com/chuangao/SMARTS
R CMD INSTALL SMARTS
```

## Example usage

``` r
library(SMARTS)
data_tmp <- sim_data(n=500, seed=123,swi_min=1.5,swi_max=4.5, wshape_bf = 1, log_hr_confound = 2)
data <- list(cont = data_tmp %>% filter(swi == 0),
            swi = data_tmp %>% filter(swi == 1))
data_rand <- random_assign(data, nbin=10, seed=123, swi_time='xoyrs', cens_time = 'lot2_eventyrs_long')
```
