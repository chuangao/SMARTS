
#' cut switching time into bins, and record the numbers
#' @param breaks a vector of quantiles used to break switching times into quantiles
#' @return a list containing the switchers, continuers, the cut of the bins, number of switchers and continuers in each bin
get_count_per_bin <- function(swi, cont, breaks, swi_time = swi_time){
  
  swi$switch_cut <- cut(swi[[swi_time]],breaks, include.lowest = T)
  n_swi <- table(swi$switch_cut)
  n_cont <- n_swi * nrow(cont)/nrow(swi)
  n_cont <- trunc(n_cont)
  n_need <- nrow(cont) - sum(n_cont)
  #set.seed(123)
  index_add <- sample(length(breaks) - 1,n_need)
  n_cont[index_add] <- n_cont[index_add] + 1
  return(list(swi = swi,cont=cont, cuts=levels(swi$switch_cut), n_swi = n_swi, n_cont = n_cont))
}

#' for switchers in each bin, assign those switching time to the right number of continuers (among all continuers), then return the continuers with pseudo-switching time. The unassigned continuers will have switching time of NA. Switchers in the Early bins are more likely to be used to assign than switchers in later bins.
#' @param counts the object returned from get_count_per_bin()
assign_per_bin <- function(counts, bin, swi_time=swi_time, cens_time = cens_time){
  
  swi <- counts$swi
  cont <- counts$cont
  
  index_to_sample_from <- with(swi,which(as.numeric(switch_cut)==bin)) ## index of switchers in this bin
  
  n_swi <- counts$n_swi[bin]
  n_cont <- counts$n_cont[bin]
  
  index_used <- which(!is.na(cont[[swi_time]])) ## index for continuers that are assigned
  index_avail <- which(is.na(cont[[swi_time]])) ## index for continuers available to assign
  print(paste0("index_avail ",length(index_avail)))
  
  n_needed <- n_cont
  
  index_used_here <- c() ## continuers assigned for this bin
  
  index_swi_useful <- c()  ## index of switchers that can find continuers
  
  n_itr <- 1
  
  if(is.infinite(-max(cont[[cens_time]][index_avail]))){
    print(length(index_avail))
    stop()
  }
  
  while(n_needed > 0 &&  !(min(swi[[swi_time]][index_to_sample_from]) > max(cont[[cens_time]][index_avail])) && n_itr <= 10){
  
    #set.seed(123)
    index_sampled <- sample(index_to_sample_from,n_needed,replace=T) ## sampling from switchers
    index_2_assign <- sample(index_avail,n_needed,replace=F) ## switchers to be assigned
    cont[[swi_time]][index_2_assign] <- swi[[swi_time]][index_sampled] ## assign
    
    ## index for continuers that are assigned unsuccessfully
    index_met <- which(cont[index_2_assign,][[swi_time]] <= cont[index_2_assign,][[cens_time]]) 
    
    index_swi_useful <- unique(c(index_swi_useful,index_sampled[index_met]))
    
    index_met <- index_2_assign[index_met] ## which assignment is valid (index in terms of all continuers)
    index_unmet <- index_2_assign[!index_2_assign %in% index_met]
    
    cont[[swi_time]][index_unmet] <- NA
    
    index_used_here <- c(index_used_here,index_met)
    index_used <- c(index_used, index_met)
    index_avail <- seq_len(nrow(cont))[-index_used]
    n_needed <- n_cont - length(index_used_here)
    
    #cat("n_needed after ", n_needed, "\n")
    
    n_itr <- n_itr + 1

  }
  return(list(cont=cont, n_needed = n_needed, index_swi_keep= unique(index_swi_useful)))
}

#' Randomly draw switching time from switchers and assign to continuers
#' @param data a list, data$cont is the continuers, data$swi is the switchers
#' @param nbin number of bins to divide the switchers into
#' @param swi_time the switching time, pre-calculated as the time length from treatment initiation to switching. For the input data, this column only exist in the switchers, a new column with the same name will be created for the continuers after assignment
#' @param seed the random seed 
#' @cens_time the follow up end time. Patients may be continuously followed after the events, this is the time for the patients to stop being followed.

#' @return a list of assigned and unassigned continuer and switchers
#' \itemize{
#' \item assigned: the continuers that are assigned switching time successfully, and the switchers that are used to assigned the continuers
#' \item unassigned: the continuers that are assigned switching times unsuccessfully, and the switchers that are never used 
#' }
#' @examples 
#' data_tmp <- sim_data(n=500, seed=123,swi_min=1.5,swi_max=4.5, wshape_bf = 1, log_hr_confound = 2)
#' data <- list(cont = data_tmp %>% filter(swi == 0),
#'              swi = data_tmp %>% filter(swi == 1))
#' data_rand <- random_assign(data, nbin=10, seed=123, swi_time='xoyrs', cens_time = 'censyrs')              
#' @export

random_assign <- function(data, nbin=10, seed = 123, swi_time='xoyrs', cens_time = 'censyrs'){
 
  set.seed(seed)
  
  cont <- data$cont
  cont[[swi_time]] = NA

  swi <- data$swi
  
  swi <- swi[sample(seq_len(nrow(swi)),nrow(swi),replace = FALSE),]
  
  ## quantiles to break the switching times into bins
  breaks <- quantile(swi[[swi_time]],prob=seq(0, 1, length.out = nbin + 1))
  breaks_bak <- breaks
  
  ## if too many bins, then reduce the number of bins
  while(length(unique(breaks)) != length(breaks)){
    breaks <- quantile(swi[[swi_time]],prob=seq(0, 1, length.out = nbin-5 + 1))
    nbin = nbin - 5
  }
  
  if(length(unique(breaks_bak)) != length(breaks)){
    warning(paste0("Too many bins provided, reduced nbin to ",nbin))
  }
  
  
  nbin <- length(breaks)-1
 
  counts <- get_count_per_bin(swi, cont, breaks, swi_time = swi_time)
  
  n_fail_cont <- c() ### continuers who failed to assign 
  swi_used <- c() ## switchers that are used
  
  ## go through each bin, and assign switching times from each bin to the continuers
  for(i in rev(seq_along(counts$cuts))){
    cat("Working on bin ", i, "\n")
    assign <- assign_per_bin(counts, i, swi_time=swi_time, cens_time = cens_time)
    counts$cont <- assign$cont
    n_fail_cont <- c(n_fail_cont, assign$n_needed)
    swi_used <- c(swi_used,assign$index_swi_keep)
  }
  
  swi <- counts$swi[unique(swi_used),]
  cont <- counts$cont[!is.na(counts$cont[[swi_time]]),]
  
  swi$swi_to_end <- swi[[cens_time]] - swi[[swi_time]]
  cont$swi_to_end <- cont[[cens_time]] - cont[[swi_time]]

  return(list(assigned=list(cont=cont, 
                            swi=swi %>% select(-switch_cut)),
         unassigned=list(cont=counts$cont[is.na(counts$cont[[swi_time]]),],
                         swi=counts$swi[-unique(swi_used),] %>% select(-switch_cut))))
  
}

#data <- random_assign(data, nbin=50, seed=123, swi_time='xoyrs', cens_time = 'censyrs')


#' simulate covariates, including ids and switching variable
#' @param n Total sample size for both continuers and switchers combined
#' @param prob probability of switchers
#' @return a list of data frames, respectively for covariate for continuers, switchers, and continuers and switchers combined
#' @export
sim_cov <- function(n=n,prob = 0.3){
  # Create a data frame with the subject IDs and treatment covariate
  cov <- data.frame(id = 1:n, swi = rbinom(n, 1, prob))
  cov_cont <- cov[cov$swi == 0,]
  n_cont <- sum(cov$swi == 0)
  
  cov_swi <- cov[cov$swi == 1,]
  n_swi <- sum(cov$swi == 1)
  
  cov <- rbind(cov[cov$swi == 0,] %>% mutate(id = 1:n_cont),
               cov[cov$swi == 1,] %>% mutate(id = (n_cont+1):(n_cont+n_swi)))
  list(cov_cont = cov_cont, cov_swi = cov_swi, cov=cov)
}

## now simulate all switchers before switching, by calling sim_switch_one_sample_bf()
sim_switch_n_sample_bf <- function(cov = cov, xoyrs=xoyrs, ...){
  #set.seed(seed)
  res <- do.call(rbind,lapply(seq_len(nrow(cov)),function(i){
    cov_i <- cov[i,,drop=F]
    xoyrs_i <- xoyrs[i]
    obj_i <- sim_switch_one_sample_bf(cov = cov_i, xoyrs = xoyrs_i, ...)
    obj_i
  }))
  res %>% mutate(
    id = cov$id
  ) %>% select(
    -c(eventtime,status)
  )
}

## switchers before switching are censored at switching time, which means each switcher has its own censor time, therefore simulate one by one 
sim_switch_one_sample_bf <- function(cov = cov_i, xoyrs = xoyrs_i, hr_bf= 2, wscale_bf = 0.05, wshape_bf=1.5){
  
  # Simulate switchers
  #set.seed(seed)
  data_bf <- simsurv(lambdas = wscale_bf, 
                     gammas = wshape_bf, 
                     betas = c(swi = log(hr_bf)), 
                     x = cov, 
                     maxt = xoyrs)
  
  data_bf <- data_bf %>% mutate(data_bf,
                                lot1_eventyrs = eventtime,
                                lot1_event = status)
  data_bf
}

## simulate switchers from initiation to end of fup, take the time from switching to end of fup as the LOT2
sim_switch_af <- function(cov = cov, xoyrs = xoyrs, hr_af= 0.7, wscale_bf = 0.05, wshape_bf=1.5, t_fow = 5, log_hr_confound = 1){
  #set.seed(seed)
  data_af <- simsurv(lambdas = wscale_bf, 
                     gammas = wshape_bf, 
                     betas = c(swi = log(hr_af),lot1_event = log_hr_confound), 
                     x = cov, 
                     maxt = t_fow)
  
  # assign a switching time
  data_af$xoyrs <- xoyrs
  
  ## only kept events after switching as lot2 event   
  data_af2 <- do.call(rbind,lapply(seq_len(nrow(data_af)),function(i){
    x = data_af[i,,drop=F]
    if(x[,"status"] == 1){
      if(x[,"eventtime"] < x[,"xoyrs"]){
        x[,"lot2_event"] = 0
        x[,"lot2_eventyrs"] = t_fow
      }else{
        x[,"lot2_event"] = 1
        x[,"lot2_eventyrs"] = x[,"eventtime"]
      }
    }else{
      x[,"lot2_event"] = 0
      x[,"lot2_eventyrs"] = t_fow
    }
    return(x)
  })) %>% as.data.frame()
  
  data_af2 %>% 
    mutate(lot2_eventyrs_long = lot2_eventyrs,
           lot2_eventyrs = lot2_eventyrs - xoyrs
    ) %>%
    select(-c(eventtime, status))
}

#' Simulate data
#' @param n Total sample size for both continuers and switchers
#' @param seed random seed
#' @param swi_min on a scale of 0 which is the treatment initiation time, to follow up end, when does the treatment switching starts to occur
#' @param swi_max on a scale of 0 which is the treatment initiation time, to follow up end, when does the treatment switching stops to occur
#' @param wshape the shape parameter for weibull distribution that is used to simulate the baseline hazard
#' @param wscale the scale parameter for weibull distribution that is used to simulate the baseline hazard
#' @param hr_bf hazard ratio before switching/pseudo-switching
#' @param hr_af hazard ratio after switching/pseudo-switching
#' @param log_hr_confound hazard ratio due to the events before switching/pseudo-switching as confounding factor
#' @param t_fow the maximum follow up end time
#' @param swi_prob probability for the switching variable to take value of 1

#' @return a data frame contining the continuers and switchers, where the variable swi = 0 indicate continuers, swi = 1 indicate switchers 

#' @examples 
#' data <- sim_data(n=500, seed=123,swi_min=1.5,swi_max=4.5, wshape_bf = 1, log_hr_confound = 2)

#' @export
sim_data <- function(n=5000, 
                     seed = 123, 
                     swi_min = 1, 
                     swi_max = 1.5, 
                     wshape = 1.5, 
                     wscale=0.05,
                     hr_bf = 2, 
                     hr_af = 0.5, 
                     log_hr_confound = 1,
                     t_fow = 6,
                     swi_prob = 0.3){
  
  
  set.seed(seed)
  
  # n=n
  # seed = seed
  # swi_min = swi_min
  # swi_max = swi_max
  # 
  # wscale = 0.05  
  # wshape = wshape
  # hr_bf = hr_bf
  # hr_af = hr_af
  # 
  # t_fow = 6
  # 
  # swi_prob = 0.3
  # 
  
  cov_list <- sim_cov(n = n, prob = swi_prob)
  #print(is.data.frame(cov))
  cov_cont <- cov_list$cov_cont
  cov_swi <- cov_list$cov_swi
  
  cov <-cov_list$cov
  n_cont <- sum(cov$swi == 0)
  n_swi <- sum(cov$swi == 1)
  
  
  # Simulate switching time from a uniform distribution
  xoyrs <- runif(n, min=swi_min, max=swi_max)
  xoyrs_cont <- xoyrs[cov$swi ==0]
  xoyrs_swi <- xoyrs[cov$swi ==1]
  
  
  ################################################################################
  ######## simulate data #########################################################
  ################################################################################
  
  set.seed(seed)
  data_cont_bf <- sim_switch_n_sample_bf(cov = cov_cont, xoyrs = xoyrs_cont, hr_bf= 1, wscale = wscale, wshape=wshape)
  
  set.seed(seed)
  data_cont_af <- sim_switch_af(cov = cbind(cov_cont,lot1_event = data_cont_bf$lot1_event), xoyrs = xoyrs_cont, hr_af= 1, wscale = wscale, wshape=wshape, t_fow = t_fow, log_hr_confound = log_hr_confound)
  
  data_cont <- cbind(data_cont_bf,data_cont_af %>% mutate (swi = 0) %>% select(-id))
  
  # data_cont <- sim_cont(cov = cov_cont, xoyrs = xoyrs_cont, wscale = wscale, wshape=wshape, t_fow = t_fow)
  set.seed(seed)
  data_swi_bf <- sim_switch_n_sample_bf(cov = cov_swi, xoyrs = xoyrs_swi, hr_bf= hr_bf, wscale = wscale, wshape=wshape)
  
  set.seed(seed)
  data_swi_af <- sim_switch_af(cov = cbind(cov_swi,lot1_event = data_swi_bf$lot1_event), xoyrs = xoyrs_swi, hr_af= hr_af, wscale = wscale, wshape=wshape, t_fow = t_fow, log_hr_confound = log_hr_confound)
  
  data_swi <- cbind(data_swi_bf,data_swi_af %>% mutate (swi = 1) %>% select(-id))
  
  data <- rbind(
    data.frame(id = 1:n_cont,data_cont %>% select(-id)), 
    data.frame(id = (n_cont+1):n, data_swi %>% select(-id)))
  
  
  data
}

## generate events for use in ITT analysis
gen_itt <- function(data){
  data <- data %>% 
    group_by(id) %>% 
    mutate(
      e_count = sum(which(c(lot1_event, lot2_event) ==1))
    ) %>% 
    ungroup()
  
  data <- data %>% mutate(
    itt1_event = case_when(
      e_count == 0  ~ 0,
      e_count != 0 ~ 1
    ),
    itt1_eventyrs = case_when(
      e_count == 0 ~ lot2_eventyrs_long,
      e_count == 1 ~ lot1_eventyrs,
      e_count == 2 ~ lot2_eventyrs_long,
      e_count == 3 ~ lot1_eventyrs
    ),
    itt2_event = case_when(
      swi == 1 ~ lot2_event,
      swi == 0 ~ itt1_event
    ),
    
    itt2_eventyrs = case_when(
      swi == 1 ~ lot2_eventyrs,
      swi == 0 ~ itt1_eventyrs
    )
  )
  
  data
}

## regenerate events based on the randomly assigned switching time
regen_event <- function(data, xoyrs = "xoyrs"){
  
  nm_events <- c("lot1_event","lot2_event")
  nm_yrs <- c("lot1_eventyrs","lot2_eventyrs_long")
  
  data_out <- lapply(c("cont","swi"),function(i){
    
    datai <- data[[i]]
    datai2 <- do.call(rbind,lapply(seq_len(nrow(datai)),function(j){
      
      x <- datai[j,,drop=F]
      
      index_events <- (x[,nm_events] == 1)
      
      index_bf <- as.numeric(x[,nm_yrs]) <= as.numeric(x[,xoyrs])
      index_af <- as.numeric(x[,nm_yrs]) > as.numeric(x[,xoyrs])
      
      if(sum(index_events & index_bf) == 0){
        x[,"lot1_event_re"] <- 0
        x[,"lot1_eventyrs_re"] <- x[,xoyrs]
      }
      if(sum(index_events & index_bf) > 0){
        x[,"lot1_event_re"] <- 1
        x[,"lot1_eventyrs_re"] <- min(x[,nm_yrs][index_events & index_bf])
      }
      
      if(sum(index_events & index_af) == 0){
        x[,"lot2_event_re"] <- 0
        x[,"lot2_eventyrs_re"] <- x[,"lot2_eventyrs_long"] - x[,xoyrs]
      }
      if(sum(index_events & index_af) > 0){
        x[,"lot2_event_re"] <- 1
        x[,"lot2_eventyrs_re"] <- min(x[,nm_yrs][index_events & index_af]) - x[,xoyrs]
      }
      x
    }))
    datai2
  })
  names(data_out) <- c("cont","swi")
  data_out
}
