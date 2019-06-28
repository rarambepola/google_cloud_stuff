##08/04/2019
##doing causal inference with cases with time lags as variables

#set working directory
setwd(paste0(Sys.getenv("HOME"), "/MDG_v3"))
 
#get z drive path
load("zdrive_path.RData")

# set.seed(27)
library(raster)

changed <- TRUE
run_no <- 3

# load(paste0("Causal_inference/ci_run_", run_no, "_N_1000_w_cases.RData"))
# seed <- 14
# set.seed(seed)
if(!changed) stop("changed seed and run number?")

ptm_start <- proc.time()
###file to implement first step of causal inference
###inferring causal structure between environmental variables
###at various time lags

##open madagascar shapefile
mdg_ecozone <- shapefile(paste0(zdrive_path, "Madagascar_PfPR/Ecozones_shp/Ecozones5.shp"))
mdg_extent <- extent(mdg_ecozone)

####set up covariates/data for step one #####
covariate_names <- c("Rain",
                     "LST_day",
                     "LST_night",
                     "TCB",
                     "EVI",
                     "TSI_pf",
                     "TSI_pv")
N_different_covs <- length(covariate_names)
time_lags <- lapply(1:N_different_covs, function(i) 0:3)

##number of iterations
N_iter <- 10
##number of observations per iteration
N_obs <- 1000
#number to blur
n_blur <- 5

N_sample_points_1 <- (N_iter + 5)*N_obs


for(time_iter in 2:3){
  print(time_iter)
  set.seed(run_no*time_iter)
  
  ###randomly sample locations and times
  ###do all at once to make extracting covariates quicker
  #locations
  locations <- c()
  #accept reject to get points actually in Madagascar
  source("../Madagascar_v2/points_in_country.R")
  #reset working directory
  setwd(paste0(Sys.getenv("HOME"), "/MDG_v3"))


  while(length(locations) < 2*N_sample_points_1){
    locations_all <- cbind(runif(3*N_sample_points_1, mdg_extent[1], mdg_extent[2]),
                           runif(3*N_sample_points_1, mdg_extent[3], mdg_extent[4]))
    locations <- rbind(locations, locations_all[points_in_country(coords=locations_all,
                                                                  country_shapefile=mdg_ecozone), ])
  }

  locations <- locations[1:N_sample_points_1, ]
  years <- sample(2013:2016, N_sample_points_1, replace=TRUE)
  months <- sample(1:12, N_sample_points_1, replace=TRUE)
  times <- months + (years - min(years)) * 12

  environment_df <- data.frame(locs=locations,
                               year=years,
                               month=months,
                               time=times)

  ptm_extract_start <- proc.time()

  # #get covariate values
  source("../Covariate_extraction/monthly.R")
  covariate_values <- monthly_covariates_individual(locations = locations,
                                         years = years,
                                         months = months,
                                         covariates = covariate_names,
                                         resolution = "5km",
                                         timelags = time_lags,
                                         crop_extent = mdg_extent,
                                         n_blur = n_blur,
                                         reverse_order = FALSE,
                                         run_parallel = TRUE,
                                         n_cores = 20)


  #find any NAs
  na_exist_1 <- unique(unlist(lapply(covariate_values, function(obs) which(is.na(obs)))))
  if(length(na_exist_1) > 0){
    covariate_values_clean <- lapply(covariate_values, "[", -na_exist_1)
    environment_df_clean <- environment_df[-na_exist_1, ]
  }else{
    covariate_values_clean <- covariate_values
    environment_df_clean <- environment_df
  }


  
  
  variable_times <- unlist(time_lags)
  N_vars <- length(variable_times)
  
  
  ####set up covariates for step 2 ####
  # load("Causal_inference/hf_obs_clean_iteration_test.RData")
  load(paste0("Causal_inference/hf_obs_clean_iteration_", time_iter, "_test_first_three_months.RData"))
  #make rates
  hf_obs_clean$"rate" <- hf_obs_clean$cases / hf_obs_clean$pop
  # N_timelags <- length(timelags)
  # for(i in timelags){
  #   hf_obs_clean <- cbind(hf_obs_clean[[paste0("cases-", i)]] / hf_obs_clean$pop, hf_obs_clean)
  #   names(hf_obs_clean)[1] <- paste0("rate-", i)
  # }
  N_sample_step_2 <- N_obs
  
  #choose some points to extract covariate values from
  N_sample_points_2 <- N_sample_step_2 * (N_iter + 3)
  sample_index <- sample.int(dim(hf_obs_clean)[1], N_sample_points_2, replace=TRUE)
  hf_sample <- hf_obs_clean[sample_index, ]
  # locations_step_2 <- cbind(hf_obs_clean$x, hf_obs_clean$y)[sample_index, ]
  # years_step_2 <- hf_obs_clean$year[sample_index]
  # months_step_2 <- hf_obs_clean$month[sample_index]
  # times_step_2 <- hf_obs_clean$time[sample_index]
  # rates_step_2 <- hf_obs_clean$rate[sample_index]
  
  #extract covariate values
  covariate_values_step_2 <- monthly_covariates_individual(locations = cbind(hf_sample$x, hf_sample$y),
                                                           years = hf_sample$year,
                                                           months = hf_sample$month,
                                                           covariates = covariate_names,
                                                           resolution = "5km",
                                                           timelags = time_lags,
                                                           crop_extent = mdg_extent,
                                                           n_blur = n_blur,
                                                           reverse_order = FALSE,
                                                           run_parallel = TRUE,
                                                           n_cores = 20)
  obs_list_step_2 <- c(list(hf_sample$rate), covariate_values_step_2)
  names(obs_list_step_2)[1] <- "rate"
  # for(i in timelags){
  #   obs_list_step_2 <- c(list(hf_sample[[paste0("rate-", i)]]), obs_list_step_2)
  #   names(obs_list_step_2)[1] <- paste0("rate-", i)
  # }
  # obs_time_lags <- rev(timelags)
  
  
  
  
  ptm_extract_end <- proc.time()
  
  print(ptm_extract_end - ptm_extract_start)
  
  #find any NAs
  na_exist <- unique(unlist(lapply(obs_list_step_2, function(obs) which(is.na(obs)))))
  if(length(na_exist) > 0){
    obs_list_step_2_clean <- lapply(obs_list_step_2, "[", -na_exist)
    hf_sample_clean <- hf_sample[-na_exist, ]
  }else{
    print("test")
    obs_list_step_2_clean <- obs_list_step_2
    hf_sample_clean <- hf_sample
  }
  
}

for(time_iter in 2:3){
  
  
  ####run pc algorithm ####
  #construct G_0 matrix
  G_0 <- matrix(1, N_vars, N_vars) - diag(N_vars)
  for(i in 1:N_vars){
    G_0[i, which(variable_times > variable_times[i])] <- 0
  }
  
  N_sample <- N_obs
  pc_list <- list()
  pc_list_2 <- list()
  alphas <- c()
  #do iterations
  source("../runpc.R")
  source("../runpc_part2.R")
  for(i in 1:N_iter){
    alpha <- runif(1, 0.9, 0.95)
    alphas[i] <- alpha
    ptm <- proc.time()
    print(i)
    pc_list[[i]] <- runpc(obsDat = covariate_values_clean,
                          locs = cbind(environment_df_clean$locs.1, environment_df_clean$locs.2),
                          times = environment_df_clean$time,
                          alpha = alpha,
                          nSample = N_sample,
                          G_0 = G_0,
                          period=12,
                          plotgraph = FALSE
    )

    print(proc.time() - ptm)
    
    #create new initial adjacency matrix from one inferred in previous step
    #now need to deal with cases with all timelags
    G_0_2 <- matrix(NA, N_vars + 1, N_vars +  1)
    # index_timelags <- 1:N_timelags
    index_rate <- 1
    
    #things can cause rates
    G_0_2[, index_rate] <- 1
    # G_0_2[, c(index_rate, index_timelags)] <- 1
    #rate cannot cause anything
    G_0_2[index_rate, ] <- 0
    
    #add in previously discovered graph between environmental variables
    G_0_2[1 + (1:N_vars), 1 + (1:N_vars)] <- pc_list[[i]][[1]]
    #allow past rates to cause current rate
    # G_0_2[index_timelags, index_rate] <- 1
    #allow rates to cause any future rates
    #but things can't cause rates if they are after in time
    # for(j in 1:N_timelags){
    #   G_0_2[index_timelags[j], which(obs_time_lags < obs_time_lags[j])] <- 1
    #   G_0_2[which(variable_times < obs_time_lags[j]) + 1 + N_timelags, index_timelags[j]] <- 0
    # }
    #nothing causes itself
    diag(G_0_2) <- 0
    
    
    # G_0_2[, 1] <- 1
    # G_0_2[1, ] <- 0
    # G_0_2[1, 1] <- 0
    # G_0_2[1 + (1:N_vars), 1 + (1:N_vars)] <- pc_list[[i]][[1]]
    
    pc_list_2[[i]] <- runpc.second(obsDat = obs_list_step_2_clean,
                                   locs = cbind(hf_sample_clean$x, hf_sample_clean$y),
                                   times = hf_sample_clean$time * 30,
                                   alpha = alpha,
                                   test.index = 1,
                                   G_0 = G_0_2,
                                   nSample=N_sample_step_2,
                                   plotgraph = FALSE
    )
  }
  
  ptm_end <- proc.time()
  print(ptm_end - ptm_start)
  
  # for(k in 1:length(pc_list_2)){
  #   print(all(pc_list[[k]][[1]] == pc_list_2[[k]][[1]][-c(1:(1+N_timelags)), -c(1:(1+N_timelags))]))
  # }
  
  plot.minimal.n(pc_list_2[[1]][[1]], 1, obsNames = names(obs_list_step_2_clean), nDeg = 1)
  
  
  pc_list_out <- lapply(pc_list, "[[", 1)
  pc_list_2_out <- lapply(pc_list_2, "[[", 1)
  save(list=c("pc_list_out", "pc_list_2_out", "alphas"),
       file=paste0("Causal_inference/time_iter_", time_iter, 
                   "_ci_run_", run_no, "_N_", N_obs, "_no_cases_first_three_months.RData")
  )
}

