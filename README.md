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
### Mock data
Mock data has been provided for demonstration purposes. The mock data contains the following columns:
id: patient ID
swi: switcher variable indicating whether a patient is a switcher or a continuer
swi_yrs: the years at which the switchers switch; the value is NA for continuers
event1 and event2: two events that occurred sequentially
event1_yrs and event2_yrs: the event years for the two events
fup_yrs: the follow-up end years

``` r
library(SMARTS)
data2 <- list(cont = data %>% filter(swi == 0),
            swi = data %>% filter(swi == 1))
data_rand <- random_assign(data2, nbin=10, seed=123, swi_time='swi_yrs', cens_time = 'fup_yrs')
```
