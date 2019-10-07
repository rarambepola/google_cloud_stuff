library(doParallel)

# create_K_inv <- function(n_rff, sigma, lambda, covariate_matrix){
#   n_covs <- dim(covariate_matrix)[2]
#   ws <- matrix(rnorm(n_rff*n_covs), n_covs, n_rff)
#   n_points <- dim(covariate_matrix)[1]
#   
#   w <- ws * sigma
#   
#   cov_mat_w <- covariate_matrix %*% w
#   # print(dim(cov_mat_w))
#   z <- matrix(NA, n_points, 2*n_rff)
#   
#   for(i in 1:n_points){
#     for(j in 1:n_rff){
#       z[i, j] <- cos(cov_mat_w[i, j]) / sqrt(n_rff)
#       z[i, j + n_rff] <- sin(cov_mat_w[i, j]) / sqrt(n_rff)
#     }
#   }
#   
#   I_n <- diag(n_points)
#   I_s <- diag(2*n_rff) * lambda
#   
#   ###woodbury matrix identiy
#   K_inv = I_n - (z %*% solve((I_s + t(z) %*% z)) %*% t(z))
# 
#   return(K_inv)
# }

# create_K_inv <- function(sigma, lambda, covariate_matrix){
#   n_covs <- dim(covariate_matrix)[2]
#   n_points <- dim(covariate_matrix)[1]
#   
#   K <- matrix(NA, n_points, n_points)
#   for(i in 1:n_points){
#     for(j in 1:n_points){
#       K[i, j] <- exp(-sum((covariate_matrix[i, ] - covariate_matrix[j, ]))/sigma)
#     }
#   }
#   return(solve(K + lambda*diag(n_points)))
# }

fit_model_mesh <- function(time_iters,
                      n_months_fit,
                      response_mat,
                      pops,
                      A,
                      spde,
                      n_s,
                      cov_mat_list,
                      sigma=0.001,
                      n_rff=10,
                      subset_use = NULL,
                      iter.max=300,
                      eval.max=300,
                      getJointPrecision=TRUE,
                      models_silent=TRUE, 
                      get_rep=FALSE,
                      max_dist=0.1,
                      mesh_n_points_start=100,
                      parallel=FALSE){
  `%myinfix%` <- ifelse(parallel, `%dopar%`, `%do%`)
  
  print("sigma")
  print(sigma)
  # print(log(sqrt(sigma)))
  fit_list <- list()
  if(get_rep) rep_list <- list()
  
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
  
  if(parallel){
    cl <- makeCluster(10)
    registerDoParallel(cl)
  }
  
  index_keep_list <- list()
  
  # time_iter <- 1
  z <- foreach(time_iter = time_iters, .packages = c("TMB", "Rcpp")) %myinfix% {
    # create_K_inv <- function(n_rff, sigma, lambda, covariate_matrix){
    #   n_covs <- dim(covariate_matrix)[2]
    #   ws <- matrix(rnorm(n_rff*n_covs), n_covs, n_rff)
    #   n_points <- dim(covariate_matrix)[1]
    #   
    #   w <- ws * sigma
    #   
    #   cov_mat_w <- covariate_matrix %*% w
    #   # print(dim(cov_mat_w))
    #   z <- matrix(NA, n_points, 2*n_rff)
    #   
    #   for(i in 1:n_points){
    #     for(j in 1:n_rff){
    #       z[i, j] <- cos(cov_mat_w[i, j]) / sqrt(n_rff)
    #       z[i, j + n_rff] <- sin(cov_mat_w[i, j]) / sqrt(n_rff)
    #     }
    #   }
    #   
    #   I_n <- diag(n_points)
    #   I_s <- diag(2*n_rff) * lambda
    #   
    #   K_inv = I_n - (z %*% solve((I_s + t(z) %*% z)) %*% t(z))
    #   
    #   return(K_inv)
    # }
    
    create_K_inv <- function(sigma, lambda, covariate_matrix){
      n_covs <- dim(covariate_matrix)[2]
      n_points <- dim(covariate_matrix)[1]
      
      K <- matrix(NA, n_points, n_points)
      for(i in 1:n_points){
        for(j in 1:n_points){
          # print(sum((covariate_matrix[i, ] - covariate_matrix[j, ])^2))
          # print(exp(-sum((covariate_matrix[i, ] - covariate_matrix[j, ])^2)/sigma))
          K[i, j] <- exp(-sum((covariate_matrix[i, ] - covariate_matrix[j, ])^2)/sigma)
        }
      }
      return(solve(K + lambda*diag(n_points)))
    }
    
    tryCatch(dyn.unload(dynlib("gp_rff_t2_fixed")),
             error = function(e) print(e))
    compile("gp_rff_t2_fixed.cpp")
    dyn.load(dynlib("gp_rff_t2_fixed"))
  # for(time_iter in 1:1){
    # print(paste0("time iter: ", time_iter))
    time_index_fit <- (time_iter - 1) + 1:n_months_fit
    cov_mat_fit <- cov_mat_list[time_index_fit]
    # print(time_index_fit)
    cases_fit <- response_mat[, time_index_fit]
    # plot(cov_mat_fit[[1]][1, ], main=time_iter)


    ##make big covariate matrix of all time points
    cov_mat_all <- do.call("rbind", cov_mat_fit)
    #get unique elmements
    cov_mat <- cov_mat_all
    n_t <- length(cov_mat_fit)
    n_locs <- dim(cov_mat_fit[[1]])[1]
    n_covs <- dim(cov_mat_fit[[1]])[2]

    n_points_u <- dim(cov_mat)[1]


    ##create a mesh
    sourceCpp("mesh_function_fast.cpp")
    index_keep <- sample(1:n_points_u, mesh_n_points_start)
    index_remove <- c()
    
    find_repeat_rows <- function(cov_matrix){
      index_remove <- c()
      n <- dim(cov_matrix)[1]
      for(i in 1:(n-1)){
        row_i <- cov_matrix[i, ]
        index_match <- ((i+1):n)[apply(cov_matrix[((i+1):n), , drop=FALSE], 1, function(row){
          all(row == row_i)
        })]
        index_remove <- c(index_remove, index_match)
      }
      return(index_remove)
    }
    
    
    
    # for(i in 1:(mesh_n_points_start-1)){
    #   index_i <- index_keep[i]
    #   index_match <- (i+1:mesh_n_points_start)
    # }
    
    index_test <- setdiff(1:n_points_u, index_keep)

    index_remove <- find_repeat_rows(cov_mat[index_test, ])
    
    if(length(index_remove) > 0){
      index_test <- index_test[-index_remove]
    }
    index_keep_cpp <- create_mesh_cpp(cov_mat, index_keep-1,
                                      index_test-1, max_dist, FALSE)
    # index_keep_list[[time_iter]] <- index_keep_cpp

    cov_mat_mesh <- cov_mat[index_keep_cpp, ]
    print(paste0("number of mesh points: ", length(index_keep_cpp)))

    #
    # print(n_locs)
    # print(n_t)
    #for each time point, identify each row in the covariate matrix
    #with its row in the big matrix
    locs_use <- matrix(NA, nrow=n_locs, ncol=n_t)
    # print(cov_mat_fit[[1]])
    # print(apply(cov_mat_fit[[1]], 1, function(cov_row){
    #   which(apply(t(t(cov_mat_all) == cov_row), 1, all))}))

    for(i in 1:n_t){
      locs_use[, i] <- apply(cov_mat_fit[[i]], 1, function(cov_row){
        # which(apply(t(t(cov_mat) == cov_row), 1, all)) - 1
        which.min(apply(cov_mat_mesh, 1, function(mesh_row) sum((mesh_row -cov_row)^2))) - 1
      }
      )
    }
    # 
    # ws <- rnorm(n_rff)
    ws <- matrix(rnorm(n_rff*n_covs), n_covs, n_rff)
    # ws <- matrix(rnorm(n_rff*n_covs), n_rff, n_covs)
    #
    # print(cov_mat %*% ws)
    # print(dim())
    # print("ws")
    # print(dim(ws))
    # print("locs_use")
    # print(dim(locs_use))
    # print("cases_fit")
    # print(dim(cases_fit))
    # print("pops")
    # print(length(pops))
    # print("A")
    # print(dim(A))
    # print("n_s")
    # print(n_s)
    # print("cov_mat")
    # print(dim(cov_mat))

    # K_inv <- create_K_inv(n_rff, sigma, lambda=1.0, covariate_matrix = cov_mat_mesh)
    
    K_inv <- create_K_inv(sigma, lambda=0.001, covariate_matrix = cov_mat_mesh)

    m <- MakeADFun(
      data = list(Y=cases_fit,
                  pops=pops,
                  cov_mat=cov_mat_mesh,
                  locs_index=locs_use,
                  K_inv=K_inv,
                  # w=ws,
                  A=A,
                  spde=spde
      ),
      parameters = list(vals=runif(length(index_keep_cpp)),
                        beta0=runif(1),
                        S=runif(n_s),
                        log_kappa=runif(1),
                        log_tau=runif(1)#,
                        # log_sigma=log(sigma),
                        # log_lambda=log(lambda)
      ),
      DLL = "gp_rff_t2_fixed",
      random="S",
      silent=models_silent
    )
  # 
  #   ###old
  #   # m <- MakeADFun(
  #   #   data = list(X=cov_mat_fit,
  #   #               Y=cases_fit,
  #   #               pops=pops,
  #   #               A=A,
  #   #               spde=spde
  #   #   ),
  #   #   parameters = list(beta0=runif(1, -1, 1),
  #   #                     beta=rep(0, N_covs),
  #   #                     S=rep(0, n_s),
  #   #                     log_kappa=0.5,
  #   #                     log_tau=0.0),
  #   #   DLL = "new_model_spatial_field",
  #   #   random = c("S"),
  #   #   silent=models_silent
  #   # )
  #   ###
  # 
    ptm <- proc.time()
    fit <- nlminb(m$par, m$fn, m$gr, control=list(iter.max=iter.max,eval.max=iter.max))
    if(get_rep) rep <- sdreport(m, getJointPrecision = getJointPrecision)
    # print(proc.time() - ptm)

    # fit_list[[time_iter]] <- fit
    # if(get_rep) rep_list[[time_iter]] <- rep

    if(!get_rep){
      list(fit, index_keep_cpp)
    }else{
      list(fit, rep, index_keep_cpp)
    }

  #    
    
    
  }


  # 
  fit_list <- lapply(z, "[[", 1)

  if(get_rep){
    rep_list <- lapply(z, "[[", 2)
    index_keep_list <- lapply(z, "[[", 3)
  }else{
    index_keep_list <- lapply(z, "[[", 2)
  }
  # 
  if(parallel) stopCluster(cl)
  # 
  if(!get_rep) return(list(fit_list, index_keep_list))
  return(list(fit_list, rep_list, index_keep_list))
}
  
  