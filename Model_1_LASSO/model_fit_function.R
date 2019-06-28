
fit_model <- function(time_iters,
                      n_months_fit,
                      response_mat,
                      pops,
                      A,
                      spde,
                      n_s,
                      cov_mat_list,
                      subset_use = NULL,
                      iter.max=300,
                      eval.max=300,
                      getJointPrecision=TRUE,
                      models_silent=TRUE, 
                      get_rep=FALSE){
  fit_list <- list()
  if(rep) rep_list <- list()
  
  #if no subset use all covariates
  if(is.null(subset_use)){
    N_covs_total <- dim(cov_mat_list[[1]])[2]
    subset_use <- 1:N_covs_total
  }
  
  N_covs <- length(subset_use)
  # print(N_covs)
  # 
  # print(subset_use)
  # print(dim(cov_mats_fit[[1]]))
  cov_mat_list <- lapply(cov_mat_list, function(mat) mat[, subset_use, drop=FALSE])
  
  
  
  for(time_iter in time_iters){
    print(paste0("time iter: ", time_iter))
    time_index_fit <- (time_iter - 1) + 1:n_months_fit
    cov_mat_fit <- cov_mat_list[time_index_fit]
    # print(time_index_fit)
    cases_fit <- response_mat[, time_index_fit]
    # plot(cov_mat_fit[[1]][1, ], main=time_iter)

    m <- MakeADFun(
      data = list(X=cov_mat_fit,
                  Y=cases_fit,
                  pops=pops,
                  A=A,
                  spde=spde
      ),
      parameters = list(beta0=runif(1, -1, 1),
                        beta=rep(0, N_covs),
                        S=rep(0, n_s),
                        log_kappa=0.5,
                        log_tau=0.0),
      DLL = "new_model_spatial_field",
      random = c("S"),
      silent=models_silent
    )

    ptm <- proc.time()
    fit <- nlminb(m$par, m$fn, m$gr, control=list(iter.max=iter.max,eval.max=iter.max))
    if(get_rep) rep <- sdreport(m, getJointPrecision = getJointPrecision)
    print(proc.time() - ptm)

    fit_list[[time_iter]] <- fit
    if(get_rep) rep_list[[time_iter]] <- rep
  }
  
  if(!get_rep) return(fit_list)
  return(list(fit_list, rep_list))
}
  
  
  
fit_model_lasso <- function(time_iters,
                      n_months_fit,
                      response_mat,
                      pops,
                      A,
                      spde,
                      n_s,
                      cov_mat_list,
                      beta_scale,
                      subset_use = NULL,
                      iter.max=300,
                      eval.max=300,
                      getJointPrecision=TRUE,
                      models_silent=TRUE){
  fit_list <- list()
  rep_list <- list()
  
  #if no subset use all covariates
  if(is.null(subset_use)){
    N_covs_total <- dim(cov_mat_list[[1]])[2]
    subset_use <- 1:N_covs_total
  }
  
  N_covs <- length(subset_use)
  # print(N_covs)
  # 
  # print(subset_use)
  # print(dim(cov_mats_fit[[1]]))
  cov_mat_list <- lapply(cov_mat_list, function(mat) mat[, subset_use, drop=FALSE])
  
  
  
  for(time_iter in time_iters){
    print(paste0("time iter: ", time_iter))
    time_index_fit <- (time_iter - 1) + 1:n_months_fit
    cov_mat_fit <- cov_mat_list[time_index_fit]
    # print(time_index_fit)
    cases_fit <- response_mat[, time_index_fit]
    # plot(cov_mat_fit[[1]][1, ], main=time_iter)
    
    m <- MakeADFun(
      data = list(X=cov_mat_fit,
                  Y=cases_fit,
                  pops=pops,
                  A=A,
                  spde=spde,
                  beta_scale=beta_scale
      ),
      parameters = list(beta0=runif(1, -1, 1),
                        beta=rep(0, N_covs),
                        S=rep(0, n_s),
                        log_kappa=0.5,
                        log_tau=0.0),
      DLL = "new_model_spatial_field_lasso",
      random = c("S"),
      silent=models_silent
    )
    
    ptm <- proc.time()
    fit <- nlminb(m$par, m$fn, m$gr, control=list(iter.max=iter.max,eval.max=iter.max))
    rep <- sdreport(m, getJointPrecision = getJointPrecision)
    print(proc.time() - ptm)
    
    fit_list[[time_iter]] <- fit
    rep_list[[time_iter]] <- rep
  }
  
  return(list(fit_list, rep_list))
}