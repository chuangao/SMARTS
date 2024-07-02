# SMARTS (Statistical Methods for Assessement of real-world Treatment Switching)

SMARTS is a tool to assign pseudo-switching time to patients who initiate a treatment and stay on the the treatment (continuers), in comparison to patients who initiate a treatment and switch to another treatment (switchers). The assigned pseudo-switching time is used in later analysis for reduce bias caused by time varying confounding.

## Install using devtools

`library(devtools)` <br/>
`install_github("chuangao/SMARTS")` <br/>

## The package can also be rom source (see below) if this is the case.  

## Install from source

`git clone https://github.com/chuangao/SMARTS` <br/>
`R CMD INSTALL SMARTS` <br/>

## Example usage
`library(SMARTS)` <br/>
`data_tmp <- sim_data(n=500, seed=123,swi_min=1.5,swi_max=4.5, wshape_bf = 1, log_hr_confound = 2)` <br/>
`data <- list(cont = data_tmp %>% filter(swi == 0)` <br/>
 `            swi = data_tmp %>% filter(swi == 1))` <br/>
`data_rand <- random_assign(data, nbin=10, seed=123, swi_time='xoyrs', cens_time = 'censyrs')` <br/>