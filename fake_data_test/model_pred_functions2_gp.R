rmse_vec <- function(x, y, na.rm=T){
  return(sqrt(mean((x-y)^2, na.rm=na.rm)))
}

rowSds <- function(mat, na.rm=T){
  return(apply(mat, 1, sd, na.rm=T))
}

standardise_mat <- function(mat, na.rm=T){
  # print(rate_actual_mat)
  # print(rate_pred_mat)
  # print(dim(mat))
  # print(rowMeans(mat))
  # print(as.matrix(mat))
  mat <- as.matrix(mat)
  return((mat - rowMeans(mat, na.rm=na.rm)) / rowSds(mat, na.rm=na.rm))
}

forward_pred_vec <- function(covs, fit, rep, A, covs_mesh){
  beta0 <- fit$par[names(fit$par) == "beta0"]
  # beta <- fit$par[names(fit$par) == "beta"]
  vals <- fit$par[names(fit$par) == "vals"]
  S <- rep$par.random[names(rep$par.random) == "S"]
  field <- A %*% S
  
  ##for each prediction point, find the nearest mesh point
  pred_mesh_index <- apply(covs, 1, function(pred_row){
    which.min(apply(covs_mesh, 1, function(mesh_row) sum((mesh_row - pred_row)^2)))
  })
  
  pred_gp <- vals[pred_mesh_index]
  
  # log_rate <- beta0 + (covs %*% beta) + field
  log_rate <- beta0 + pred_gp + field
  rate_pred <- exp(log_rate)
  
  return(as.matrix(rate_pred))
}

forward_pred <- function(fit, covs, rep, A, covs_mesh, as_matrix=TRUE){
  if(as_matrix){
    return(do.call("cbind", lapply(covs, forward_pred_vec, fit, rep, A, covs_mesh)))
  }else{
    return(lapply(covs, forward_pred_vec, fit, rep, A, covs_mesh))
  }
}

metrics <- function(rate_pred_mat, rate_actual_mat){
  actual_raw <- as.vector(rate_actual_mat)
  pred_raw <- as.vector(rate_pred_mat)
  

  actual_st <- as.vector(standardise_mat(rate_actual_mat))
  pred_st <- as.vector(standardise_mat(rate_pred_mat))
  
  metric_list <- c()
  metric_list[1] <- cor(actual_st, pred_st, use = "complete")
  metric_list[2] <- cor(actual_raw, pred_raw, use = "complete")
  metric_list[3] <- rmse_vec(actual_st, pred_st, na.rm=T)
  metric_list[4] <- rmse_vec(actual_raw, pred_raw, na.rm=T)
  metric_list[5] <- metric_list[3] / mean(actual_st, na.rm=T)
  metric_list[6] <- metric_list[4] / mean(actual_raw, na.rm=T)
  
  names(metric_list) <- c("cor_st", "cor_raw", "rmse_st", "rmse_raw", "prop_rmse_st", "prop_rmse_raw")
  return(metric_list)
}

metrics_full <- function(rate_pred_mat_list, rate_actual_mat_list,
                         subset_index = NULL){
  if(!is.null(subset_index)){
    rate_pred_mat_list <- lapply(rate_pred_mat_list, function(mat) mat[subset_index, ])
    rate_actual_mat_list <- lapply(rate_actual_mat_list, function(mat) mat[subset_index, ])
  }
  
  metrics_list <- lapply(1:length(rate_pred_mat_list), function(i) metrics(rate_pred_mat_list[[i]],
                                                                           rate_actual_mat_list[[i]]))
  
  metrics_out <- list()
  for(i in 1:length(metrics_list[[1]])){
    metrics_out[[i]] <- sapply(metrics_list, "[", i)
  }
  
  names(metrics_out) <- names(metrics_list[[1]])
  
  return(metrics_out)
}



forward_predict_full <- function(time_iters, 
                                 fit_list, 
                                 cov_mat_list,
                                 cov_mat_fit_list,
                                 rep_list,
                                 mesh_index_list,
                                 n_months_fit, 
                                 n_months_pred,
                                 A,
                                 rate_actual_mat_full,
                                 N_hf_fit,
                                 N_hf_pred,
                                 A_fit,
                                 subset_use=NULL){
  fit_index <- 1:N_hf_fit
  pred_index <- N_hf_fit + (1:N_hf_pred)
  
  # print(dim(cov_mat_list[[1]]))
  if(is.null(subset_use)){
    N_covs_total <- dim(cov_mat_list[[n_months_fit + 1]])[2]
    subset_use <- 1:N_covs_total
  }
  
  if(is.list(subset_use)){
    subset_use_list <- subset_use
  }else{
    subset_use_list <- list()
    for(time_iter in time_iters){
      subset_use_list[[time_iter]] <- subset_use
    }
  }
  
  rate_pred_list <- list()
  rate_actual_list <- list()
  pred_st_list <- list()
  actual_st_list <- list()
  metrics_list <- list()
  
  rate_fit_list <- list()
  

  # for(time_iter in 1:N_time_iters){
  for(time_iter in time_iters){
    print(paste0("time_iter: ", time_iter))
    subset_use <- subset_use_list[[time_iter]]
    N_covs <- length(subset_use)
    pred_covs <- lapply(cov_mat_list, function(mat) mat[, subset_use, drop=FALSE])
    fit_covs <- lapply(cov_mat_fit_list, function(mat) mat[, subset_use, drop=FALSE])
    time_index_fit <- (time_iter - 1) + 1:n_months_fit
    time_index_pred <- (time_iter - 1) + n_months_fit + 1:n_months_pred
    
    cov_fit_all <- do.call("rbind", fit_covs[time_index_fit])
    
    # print(time_index_fit)
    # print(dim(cov_mat_fit_list[[1]]))
    # print(dim(cov_fit_all))
    # print(max(mesh_index_list[[time_iter]]))
    covs_mesh <- cov_fit_all[mesh_index_list[[time_iter]], ]
    

    rate_pred_list[[time_iter]] <- forward_pred(fit_list[[time_iter]], pred_covs[time_index_pred], 
                                                rep_list[[time_iter]], A, covs_mesh)
    
    rate_actual_list[[time_iter]] <- rate_actual_mat_full[, time_index_pred]
    
    
    ##
    rate_fit_list[[time_iter]] <- forward_pred(fit_list[[time_iter]], fit_covs[time_index_fit],
                                               rep_list[[time_iter]], A_fit, covs_mesh)
    
    # metrics_list[[time_iter]] <- metrics(rate_pred_list[[time_iter]],
    #                                      rate_actual_list[[time_iter]])
    
    # actual_st <- as.vector(standardise_mat(rate_actual_list[[time_iter]]))
    # pred_st <- as.vector(standardise_mat(rate_pred_list[[time_iter]]))
    
    actual_st <- standardise_mat(rate_actual_list[[time_iter]])
    pred_st <- standardise_mat(rate_pred_list[[time_iter]])
    
    actual_st_list[[time_iter]] <- actual_st
    pred_st_list[[time_iter]] <- pred_st
  }
  
  metrics_list <- metrics_full(rate_pred_list, rate_actual_list)
  
  metrics_fit <- metrics_full(rate_pred_list, rate_actual_list,
                              fit_index)
  metrics_new <- metrics_full(rate_pred_list, rate_actual_list,
                              pred_index)
  

  ##aggregate to region and district
  load("poly_list.RData")
  source("spatial_aggregation_functions.R")
  district_list <- aggregate_calc(rate_pred_list,
                                  rate_actual_list,
                                  poly_list[[1]],
                                  fit_index,
                                  pred_index)
  region_list <- aggregate_calc(rate_pred_list,
                                  rate_actual_list,
                                  poly_list[[2]],
                                  fit_index,
                                  pred_index)
  output_list <- list(rate_pred_list,
                      rate_actual_list,
                      pred_st_list,
                      actual_st_list,
                      metrics_list,
                      metrics_fit,
                      metrics_new,
                      district_list,
                      region_list,
                      rate_fit_list)
  
  names(output_list) <- c("rate_pred", "rate_actual", "pred_st", 
                          "actual_st", "metrics", "metrics_fit", "metrics_new",
                          "district_list",
                          "region_list",
                          "rate_fit_list")
  return(output_list)
}


aggregate_calc <- function(rate_pred_list, rate_actual_list, poly_list,
                           fit_index,
                           pred_index){
  poly_all <- poly_list[[1]]
  poly_fit <- poly_list[[2]]
  poly_pred <- poly_list[[3]]
  
  rate_pred_agg_all <- list()
  rate_actual_agg_all <- list()
  rate_pred_agg_fit <- list()
  rate_actual_agg_fit <- list()
  rate_pred_agg_pred <- list()
  rate_actual_agg_pred <- list()
  
  pred_std_agg_all  <- list()
  actual_std_agg_all <- list()
  pred_std_agg_fit <- list()
  actual_std_agg_fit <- list()
  pred_std_agg_pred <- list()
  actual_std_agg_pred <- list()
  
  N <- length(rate_pred_list)
  for(i in 1:N){
    rate_pred_agg_all[[i]] <- case_aggregate_matrix(rate_pred_list[[i]], poly_all)
    rate_pred_agg_fit[[i]] <- case_aggregate_matrix(rate_pred_list[[i]][fit_index, ], poly_fit)
    rate_pred_agg_pred[[i]] <- case_aggregate_matrix(rate_pred_list[[i]][pred_index, ], poly_pred)
    
    rate_actual_agg_all[[i]] <- case_aggregate_matrix(rate_actual_list[[i]], poly_all)
    rate_actual_agg_fit[[i]] <- case_aggregate_matrix(rate_actual_list[[i]][fit_index, ], poly_fit)
    rate_actual_agg_pred[[i]] <- case_aggregate_matrix(rate_actual_list[[i]][pred_index, ], poly_pred)
    
    pred_std_agg_all[[i]]  <- standardise_mat(rate_pred_agg_all[[i]])
    actual_std_agg_all[[i]] <- standardise_mat(rate_actual_agg_all[[i]])
    pred_std_agg_fit[[i]] <- standardise_mat(rate_pred_agg_fit[[i]])
    
    actual_std_agg_fit[[i]] <- standardise_mat(rate_actual_agg_fit[[i]])
    pred_std_agg_pred[[i]] <- standardise_mat(rate_pred_agg_pred[[i]])
    actual_std_agg_pred[[i]] <- standardise_mat(rate_actual_agg_pred[[i]])
  }
  
  ##do metrics
  metrics_all <- metrics_full(rate_pred_agg_all, rate_actual_agg_all)
  metrics_fit <- metrics_full(rate_pred_agg_fit, rate_actual_agg_fit)
  metrics_pred <- metrics_full(rate_pred_agg_pred, rate_actual_agg_pred)
  
  rate_pred_list <- list(rate_pred_agg_all, 
                         rate_pred_agg_fit,
                         rate_pred_agg_pred)
  names(rate_pred_list) <- c("all", "fit", "pred")
  
  rate_actual_list <- list(rate_actual_agg_all, 
                         rate_actual_agg_fit,
                         rate_actual_agg_pred)
  names(rate_actual_list) <- c("all", "fit", "pred")
  
  rate_std_actual_list <- list(actual_std_agg_all,
                               actual_std_agg_fit,
                               actual_std_agg_pred)
  names(rate_std_actual_list) <- c("all", "fit", "pred")
  
  rate_std_pred_list <- list(pred_std_agg_all,
                             pred_std_agg_fit,
                             pred_std_agg_pred)
  names(rate_std_pred_list) <- c("all", "fit", "pred")
  
  rate_std_list <- list(rate_std_actual_list, rate_std_pred_list)
  names(rate_std_list) <- c("actual", "pred")

  rate_list <- list(rate_pred_list, rate_actual_list)
  names(rate_list) <- c("rate_pred", "rate_actual")
  
  metrics_list <- list(metrics_all, 
                      metrics_fit,
                      metrics_pred)
  names(metrics_list) <- c("all", "fit", "pred")
  
  
  out_list <- list(rate_list, metrics_list, rate_std_list)
  names(out_list) <- c("rate", "metrics", "rate_std")
  return(out_list)
}
