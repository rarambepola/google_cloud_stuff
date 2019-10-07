path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")
setwd(path_input)
print(getwd())

# 
# #set working directory
# setwd(paste0(Sys.getenv("HOME"), "/v4_MDG/GP"))
# 
# #get z drive path
# load("../zdrive_path.RData")

library(TMB)
library(INLA)

tryCatch(dyn.unload(dynlib("gp_rff_t2")),
         error = function(e) print(e))
compile("gp_rff_t2.cpp")
dyn.load(dynlib("gp_rff_t2"))

# 
# rmse_vec <- function(x, y, na.rm=T){
#   return(sqrt(mean((x-y)^2, na.rm=na.rm)))
# }
# 
# # load("../covariates.RData")
# load("../../v4_MDG_big_files/covariates_all.RData")
# 
# #make mesh and spde object
# mesh <- inla.mesh.2d(loc = coords_fit, cutoff = 0.1, max.edge = 3) 
# plot(mesh)
# points(coords_fit, col="red")
# spde <- (inla.spde2.matern(mesh=mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
# 
# A_pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coords_pred))
# n_s <- nrow(spde$M0)
# mesh_coords <- mesh$loc[, 1:2]
# 
# 
# # cov_subset <- 1:12
# cov_subset <- c(1:4, 20)
# # covs_use <- 1:28
# # covs_use <- c(5:8, 20)
# #add some fake covariates
# 
# n_hf_fit <- dim(all_covs$covmats_normalised[[1]])[1]
# n_hf_pred <- dim(pred_covs$covmats_normalised[[1]])[1]
# 
# models_silent <- FALSE
# n_months_pred <- 24
# 
# # all_cov_cor <- c()
# # causal_cov_cor <- c()
# # lasso_cov_cor <- c()
# # subset_list <- list(1:28,
# #                     c(1:4, 20),
# #                     c(5:8, 20)
# #                     )
# # N_subsets <- length(subset_list)
# subset_names <- c("all", "causal", "lasso")
# 
# 
# # set.seed(10)
# time_points_use <- 1:2
# 
# locs_use <- sample.int(560, 400)
# # locs_use <- 1:50
# n_points <- length(locs_use) * length(time_points_use)
# 
# #get covariate matrix
# cov_mat_all <- c()
# coords_fit_use <- c()
# for(time_point_use in time_points_use){
#   cov_mat_all <- rbind(cov_mat_all, all_covs$covmats_normalised[[time_point_use]][locs_use, cov_subset])
#   coords_fit_use <- rbind(coords_fit_use, coords_fit[locs_use, ])
# }
# A_fit <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coords_fit_use))
# 
# # cov_mat_all <- all_covs$covmats_normalised[[time_point_use]][locs_use, cov_subset]
# cov_mat <- unique(cov_mat_all)
# 
# n_covs <- dim(cov_mat)[2]
# 
# ##find out which locations to associate with row of covariate matrix
# locs_index <- matrix(NA, length(locs_use), length(time_points_use))
# for(j in 1:length(time_points_use)){
#   time_point_use <- time_points_use[j]
#   cov_mat_t <- all_covs$covmats_normalised[[time_point_use]][locs_use, cov_subset]
#   for(i in 1:length(locs_use)){
#     locs_index[i, j] <- which(apply(cov_mat, 1, function(row) all(row == cov_mat_t[i, ]))) - 1
#   }
#   
# }
# 
# 
# 
# # cov_mat <- all_covs$covmats_normalised[[time_point_use]][locs_use, ]
# 
# ##create covariance matrix
# n_points_u <- dim(cov_mat)[1]
# K <- matrix(NA, n_points_u, n_points_u)
# 
# sigma <- 1
# sigma2 <- 1
# for(i in 1:n_points_u){
#   for(j in 1:n_points_u){
#     K[i, j] <- exp(- sum((cov_mat[i, ] - cov_mat[j, ])^2) / sigma)
#   }
# }
# 
# K <- K*sigma2
# 
# ##invert K
# # K_inv <- solve(K)
# 
# # K_inv <- diag(n_points)
# lambda=10
# sigma=10
# D <- 10
# ws <- matrix(rnorm(D, sd=1/sigma), n_covs, D)
# 
# rate_use <- Y_fit[locs_use, time_points_use] / pops_fit[locs_use]
# 
# Y <- as.vector(Y_fit[locs_use, time_points_use])
# Y_t <- Y_fit[locs_use, time_points_use]
# cov_pred <- pred_covs$covmats_normalised[[time_point_use]][, cov_subset]
# 
# 
# # save(list=c("Y_t", "pops_fit", "locs_use", "time_points_use", "cov_mat", "locs_index", "ws",
# #             "A_fit", "spde", "n_points_u", "sigma", "lambda", "models_silent", "A_pred", 
# #             "pops_pred", "Y_pred"),
# #      file="gp_rff_t_pred.RData")

load("gp_rff_t_pred.RData")

sigma <- 0.001
sqrt_sigma <- sqrt(sigma)

m <- MakeADFun(
  data = list(Y=Y_t,
              pops=rep(pops_fit[locs_use], length(time_points_use)),
              cov_mat=cov_mat,
              locs_index=locs_index,
              log_sigma=log(sqrt(sigma)),
              w=ws,
              A=A_fit, 
              spde=spde
  ),
  parameters = list(vals=runif(n_points_u),
                    beta0=runif(1),
                    S=runif(n_s),
                    log_kappa=runif(1),
                    log_tau=runif(1)#,
                    # log_sigma=log(sigma),
                    # log_lambda=log(lambda)
  ),
  DLL = "gp_rff_t2",
  random="S",
  silent=models_silent
)

ptm <- proc.time()
fit <- nlminb(m$par, m$fn, m$gr, control=list(iter.max=300,eval.max=300))
rep <- sdreport(m)
print(proc.time() - ptm)

path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")

save(list=c("fit"),
     file=paste0(path_output, "/gp_rff_t_fit2_sigma_", sigma, ".RData"))

save(list=c("rep"),
     file=paste0(path_output, "/gp_rff_t_rep2_sigma_", sigma, ".RData"))

# 
# plot.new()
# par(mfrow=c(1,1))
# 
# vals <- fit$par[names(fit$par) == "vals"]
# beta0 <- fit$par[names(fit$par) == "beta0"]
# S <- rep$par.random[names(rep$par.random) == "S"]
# 
# # plot(rate_use, exp(beta0 + vals))
# # plot(rate_use, exp(beta0+vals)[as.vector(t(locs_index))])
# plot(Y_t[, 1] / pops_fit[locs_use], exp(beta0+vals[locs_index[, 1] + 1]))
# abline(0, 1)
# plot(Y_t[, 2] / pops_fit[locs_use], exp(beta0+vals[locs_index[, 2] + 1]))
# abline(0, 1)
# # plot(rate_use)
# 
# ##make some predictions
# 
# n_pred <- dim(cov_pred)[1]
# sigma_11 <- matrix(NA, n_pred, n_pred)
# 
# for(i in 1:n_pred){
#   for(j in 1:n_pred){
#     sigma_11[i, j] <- exp(- sum((cov_pred[i, ] - cov_pred[j, ])^2) / sigma)
#   }
# }
# 
# sigma_12 <- matrix(NA, n_pred, n_points_u)
# 
# for(i in 1:n_pred){
#   for(j in 1:n_points_u){
#     sigma_12[i, j] <- exp(- sum((cov_pred[i, ] - cov_mat[j, ])^2) / sigma)
#   }
# }
# 
# sigma_21 <- t(sigma_12)
# 
# field_pred <- A_pred %*% S
# mu_new <- sigma_12 %*% solve((K + lambda*diag(n_points_u)), vals)
# plot(exp(field_pred + beta0 + mu_new), Y_pred[, time_point_use] / pops_pred)
# abline(0,1)
# print(as.vector(cor(exp(as.vector(field_pred) + beta0+mu_new), Y_pred[, time_point_use] / pops_pred, use="complete")))
# 
