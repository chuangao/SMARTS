#install.packages("CBPS")
#library(CBPS)
# 
# library(simsurv)
# library(MatchIt)
# library(tidyverse)
# library(flexsurv)
# library(survival)
# library(survminer)
# library(WeightIt)
# library(parallel)

#' cut switching time into bins, and record the numbers
#' return a list containing the switchers, continuers, the cut of the bins, number of switchers and continuers in each bin
get_count_per_bin <- function(swi, cont, breaks, swi_time = swi_time){
  #print(swi[1,])
  #print(swi_time)
  #swi_time <- enexpr(!!swi_time)
  
  #print(swi[[swi_time]])
  
  swi$switch_cut <- cut(swi[[swi_time]],breaks, include.lowest = T)
  n_swi <- table(swi$switch_cut)
  n_cont <- n_swi * nrow(cont)/nrow(swi)
  n_cont <- trunc(n_cont)
  n_need <- nrow(cont) - sum(n_cont)
  #set.seed(123)
  index_add <- sample(length(breaks) - 1,n_need)
  n_cont[index_add] <- n_cont[index_add] + 1
  return(list(swi = swi,
              cont=cont, 
              cuts=levels(swi$switch_cut), 
              n_swi = n_swi, 
              n_cont = n_cont,
              #n_assigned_after= 0 * n_cont,
              #n_assigned_pre = 0 * n_cont,
              #n_itr = 0 * n_cont,
              n_needed = n_cont
              #index_swi_keep = vector("list",length(n_cont))
  ))
}

#' for switchers in each bin, assign those switching time to the right number of continuers (among all continuers), then return the assigned continuers, the unassigned continuers for the next round of assign. Early bins is more likely to be assigned than later bins
assign_per_bin <- function(counts, bin, swi_time=swi_time, cens_time = cens_time){
  
  #swi_time <- enexpr(swi_time)
  #cens_time <- enexpr(cens_time)
  
  swi <- counts$swi
  cont <- counts$cont
  
  index_to_sample_from <- with(swi,which(as.numeric(switch_cut)==bin)) ## index of switchers in this bin
  
  n_swi <- counts$n_swi[bin]
  n_cont <- counts$n_cont[bin]
  
  index_used <- which(!is.na(cont[[swi_time]])) ## continuers that are assigned
  index_avail <- which(is.na(cont[[swi_time]])) ## continuers available to assign
  #print(paste0("index_avail ",length(index_avail)))
  
  n_needed <- counts$n_needed
  #print("counts$n_needed in assign per bin")
  #print(counts$n_needed)
  
  index_used_here <- c() ## continuers assigned for this bin
  
  index_swi_useful <- c()  ## index of switchers that can find continuers
  
  #n_itr <- 1
  
  if(is.infinite(-max(cont[[cens_time]][index_avail]))){
    print(length(index_avail))
    stop()
  }
  
  #while(n_needed > 0 &&  !(min(swi[[swi_time]][index_to_sample_from]) > max(cont[[cens_time]][index_avail])) && n_itr <= 10){
  #cat("iterations ", n_itr, "\n")
  #cat("n_needed before ", n_needed, "\n")
  
  #set.seed(123)
  index_sampled <- sample(index_to_sample_from,n_needed[bin],replace=T) ## sampling from switchers
  index_2_assign <- sample(index_avail,n_needed[bin],replace=F) ## switchers to be assigned
  cont[[swi_time]][index_2_assign] <- swi[[swi_time]][index_sampled] ## assign
  
  ## which assignment is valid (index in terms of sampled samples)
  index_met <- which(cont[index_2_assign,][[swi_time]] <= cont[index_2_assign,][[cens_time]]) 
  
  index_swi_useful <- unique(c(index_swi_useful,index_sampled[index_met]))
  
  index_met <- index_2_assign[index_met] ## which assignment is valid (index in terms of all continuers)
  index_unmet <- index_2_assign[!index_2_assign %in% index_met]
  
  cont[[swi_time]][index_unmet] <- NA
  
  index_used_here <- c(index_used_here,index_met) ## continuers used by this bins
  index_used <- c(index_used, index_met)    ## continuers used by all bins
  index_avail <- seq_len(nrow(cont))[-index_used]  
  n_needed <- n_needed[bin] - length(index_used_here)   ## how many continuers can this bin still assign to
  
  #cat("n_needed after ", n_needed, "\n")
  
  #n_itr <- n_itr + 1
  
  #}
  return(list(cont=cont, n_needed = n_needed, index_swi_keep= unique(index_swi_useful)))
}

#' Randomly draw switching time from switchers and assign to continuers
#' @param data a list, data$cont is the continuers, data$swi is the switchers
#' @param nbin number of bins to divide the switchers into
#' @param swi_time the switching time, pre-calculated as the time length from treatment initiation to switching. For the input data, this column only exist in the switchers, a new column with the same name will be created for the continuers after assignment
#' @return a list of assigned and unassigned continuer and switchers
#' \itemize{
#' \item assigned: the continuers that are assigned switching time successfully, and the switchers that are used to assigned the continuers
#' \item unassigned: the continuers that are assigned switching times unsuccessfully, and the switchers that are never used 
#' }
#' @examples 
#' data_tmp <- sim_switch(seed = 123, n = 5000)
#' data <- list(cont = data_tmp %>% filter(swi == 0),
#'              swi = data_tmp %>% filter(swi == 1))
#' data_rand <- random_assign(data, nbin=10, seed=123, swi_time='swi_rs', cens_time = 'fup_yrs')              
#' @export

random_assign <- function(data, nbin=10, seed = 123, swi_time='swi_yrs', cens_time = 'fup_yrs'){
  
  set.seed(seed)
  
  #swi_time <- enexpr(swi_time)
  #cens_time <- enexpr(cens_time)
  
  cont <- data$cont
  cont[[swi_time]] = NA
  
  swi <- data$swi
  
  #swi <- swi[order(swi[[swi_time]]),]
  #swi <- swi[sample(seq_len(nrow(swi)),nrow(swi),replace = FALSE),]
  
  breaks <- quantile(swi[[swi_time]],prob=seq(0, 1, length.out = nbin + 1))
  breaks_bak <- breaks
  
  while(length(unique(breaks)) != length(breaks)){
    breaks <- quantile(swi[[swi_time]],prob=seq(0, 1, length.out = nbin-5 + 1))
    nbin = nbin - 5
  }
  
  if(length(unique(breaks_bak)) != length(breaks)){
    warning(paste0("Too many bins provided, reduced nbin to ",nbin))
  }
  
  nbin <- length(breaks)-1
  
  counts <- get_count_per_bin(swi, cont, breaks, swi_time = swi_time)
  
  #n_fail_cont <- c() ### continuers who failed to assign 
  #swi_used <- c() ## switchers that are used
  
  n_cont_needed_after <- counts$n_cont ## how many quota of continuers are switch bin left to assign to 
  n_cont_needed_pre <- counts$n_cont 
  
  n_itr_total <- 0 * counts$n_cont
  index_swi_keep <- rep(list(numeric(0)),nbin)
  
  while(any(n_itr_total < 20)){
    for(i in rev(seq_along(counts$cuts))){
      #cat("Working on bin ", i, "\n")
      #print(i)
      assign <- assign_per_bin(counts, i, swi_time=swi_time, cens_time = cens_time)
      #print(assign$n_needed)
      counts$cont <- assign$cont
      counts$n_needed[i] <- assign$n_needed 
      
      n_cont_needed_after[i] <- assign$n_needed
      index_swi_keep[[i]] <- c(index_swi_keep[[i]], assign$index_swi_keep)
      if(n_cont_needed_after[i] < n_cont_needed_pre[i]){
        n_itr_total[i] = 0
      }else{
        n_itr_total[i] = n_itr_total[i] + 1
      }
      n_cont_needed_pre[i] <- n_cont_needed_after[i]
      #print("counts$n_needed in random assign")
      #print(counts$n_needed)
    }
    
  }
  
  #index_swi_keep counts$swi$xoyrs
  
  index_swi_keep_all <- unique(unlist(index_swi_keep))
  swi <- counts$swi[index_swi_keep_all,]
  cont <- counts$cont[!is.na(counts$cont[[swi_time]]),]
  
  swi$swi_to_end <- swi[[cens_time]] - swi[[swi_time]]
  cont$swi_to_end <- cont[[cens_time]] - cont[[swi_time]]
  
  return(list(assigned=list(cont=cont, 
                            swi=swi %>% select(-switch_cut)),
              unassigned=list(cont=counts$cont[is.na(counts$cont[[swi_time]]),],
                              swi=counts$swi[-unique(index_swi_keep_all),] %>% select(-switch_cut))))
  
}
