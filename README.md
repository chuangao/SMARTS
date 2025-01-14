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

## How to use
### Mock data
Mock data has been provided for demonstration purposes. The mock data contains the following columns:  
id: patient ID.   
swi: switcher variable indicating whether a patient is a switcher or a continuer.  
swi_yrs: the years at which the switchers switch; the value is NA for continuers.  
event1 and event2: two events that occurred sequentially.  
event1_yrs and event2_yrs: the event years for the two events.  
fup_yrs: the follow-up end years.  

``` r
library(SMARTS)
library(tidyverse)
library(MatchIt)
library(WeightIt)
head(mock_data)
```

### Fit an ITT model 
To fit the model, the time to event outcome was derived from treatment initiation/switching to follow up end for the continuers/switchers respectively.
``` r
# generate ITT event
data <- mock_data %>% 
  group_by(id) %>% 
  mutate(
    e_count = sum(which(c(event1, event2) ==1))
  ) %>% 
  ungroup()

data <- data %>% mutate(
  itt_event = case_when(
    e_count == 0  ~ 0,
    e_count != 0 ~ 1
  ),
  itt_eventyrs = case_when(
    e_count == 0 ~ fup_yrs,
    e_count == 1 ~ event1_yrs,
    e_count == 2 ~ event2_yrs,
    e_count == 3 ~ event1_yrs
  ))
# fit coxph model
coxph(Surv(time = itt_eventyrs, event = itt_event) ~ as.factor(swi), data = data)
```
### Assign pseudo-switching time
``` r
library(SMARTS)
# convert data to a list of continuers and switchers
data <- list(cont = data %>% filter(swi == 0),
            swi = data %>% filter(swi == 1))
# assign pseudo-switching times
data2 <- random_assign(data, nbin=10, seed=123, swi_time='swi_yrs', cens_time = 'fup_yrs')
```

### check distribution of the switching and the pseudo-switching time
```r
qqplot(data2$assigned$swi$swi_yrs,data2$assigned$cont$swi_yrs)
abline(0,1)
```
### Generate time to event outcome pre-switching and after-switching
``` r
####### regenerate event pre-switching (lot1) and after_switching (lot2) based on the pseudo-switching and switching time
data2 <- rbind(data2$assigned$cont,data2$assigned$swi)
data2 <- data2 %>%
  mutate(
    event1_yrs = ifelse(event1 == 1,event1_yrs,swi_yrs),
    event2_yrs = ifelse(event2 == 1,event2_yrs,fup_yrs),
  )

nm_events <- c("event1","event2")
nm_yrs <- c("event1_yrs","event2_yrs")

data3 <- do.call(rbind,lapply(seq_len(nrow(data2)),function(j){
  
  x <- data2[j,,drop=F]
  
  index_events <- (x[,nm_events] == 1)
  
  index_bf <- as.numeric(x[,nm_yrs]) <= as.numeric(x[,"swi_yrs"])
  index_af <- as.numeric(x[,nm_yrs]) > as.numeric(x[,"swi_yrs"])
  
  if(sum(index_events & index_bf) == 0){
    x[,"lot1_event"] <- 0
    x[,"lot1_eventyrs"] <- x[,"swi_yrs"]
  }
  if(sum(index_events & index_bf) > 0){
    x[,"lot1_event"] <- 1
    x[,"lot1_eventyrs"] <- min(x[,nm_yrs][index_events & index_bf])
  }
  
  if(sum(index_events & index_af) == 0){
    x[,"lot2_event"] <- 0
    x[,"lot2_eventyrs"] <- x[,"event2_yrs"] - x[,"swi_yrs"]
  }
  if(sum(index_events & index_af) > 0){
    x[,"lot2_event"] <- 1
    x[,"lot2_eventyrs"] <- min(x[,nm_yrs][index_events & index_af]) - x[,"swi_yrs"]
  }
  x
}))
```
### Weight the cohorts on pre-switching events, and fit the model on after-switching events
```r
 model <- paste0("swi ~ ",paste(c("swi_yrs","lot1_event"),collapse="+"))

 weight_ate <- WeightIt::weightit(formula(model),
                                   data3,
                                   method = "ps",
                                   estimand = "ATE"
  )
  
  coxph(Surv(time = lot2_eventyrs, event = lot2_event) ~ as.factor(swi), data = data3, weight = weight_ate$weight)
```

