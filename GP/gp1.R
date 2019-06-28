path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")
setwd(path_input)
print(getwd())

# #set working directory
# setwd(paste0(Sys.getenv("HOME"), "/v4_MDG/GP"))
# 
# #get z drive path
# load("../zdrive_path.RData")
# 
# library(TMB)
# library(INLA)
# 
# rmse_vec <- function(x, y, na.rm=T){
#   return(sqrt(mean((x-y)^2, na.rm=na.rm)))
# }
# 
# cov_subset <- 1:12
# # covs_use <- c(1:4, 20)
# # covs_use <- 1:28
# # covs_use <- c(5:8, 20)
# #add some fake covariates
# 
# # load("../covariates.RData")
# load("../../v4_MDG_big_files/covariates_all.RData")
# 
# n_hf_fit <- dim(all_covs$covmats_normalised[[1]])[1]
# n_hf_pred <- dim(pred_covs$covmats_normalised[[1]])[1]
# 
# models_silent <- TRUE
# 
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
# time_point_use <- 4
# locs_use <- sample.int(560, 400)
# # locs_use <- 1:50
# n_points <- length(locs_use)
# 
# #get covariate matrix
# cov_mat_all <- all_covs$covmats_normalised[[time_point_use]][locs_use, cov_subset]
# cov_mat <- unique(cov_mat_all)
# 
# ##find out which locations to associate with row of covariate matrix
# locs_index <- c()
# for(i in 1:n_points){
#   locs_index[i] <- which(apply(cov_mat, 1, function(row) all(row == cov_mat_all[i, ])))
# }
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
# K_inv <- solve(K)
# 
# # K_inv <- diag(n_points)
# 
# 
# 
# tryCatch(dyn.unload(dynlib("gp1")),
#          error = function(e) print(e))
# compile("gp1.cpp")
# dyn.load(dynlib("gp1"))
# 
# rate_use <- Y_fit[locs_use, time_point_use] / pops_fit[locs_use]
# 
# cov_pred <- pred_covs$covmats_normalised[[time_point_use]][, cov_subset]
# 
# save(list = c("rate_use", "Y_fit", "time_point_use", "locs_use", "pops_fit",
#               "K_inv", "locs_index", "cov_pred", "n_points_u", "n_points", "cov_mat",
#               "Y_pred", "pops_pred", "K"),
#      file="gp1_prep.RData")

load("gp1_prep.RData")
m <- MakeADFun(
  data = list(Y=Y_fit[locs_use, time_point_use],
              pops=pops_fit[locs_use],
              K_inv=K_inv,
              locs_index=locs_index-1
  ),
  parameters = list(vals=runif(dim(K_inv)[1])
                    # ,beta0=runif(1)
  ),
  DLL = "gp1",
  silent=models_silent
)

ptm <- proc.time()
fit <- nlminb(m$par, m$fn, m$gr, control=list(iter.max=300,eval.max=300))
# rep <- sdreport(m)
print(proc.time() - ptm)

save(list=c("fit"),
     file=paste0(path_output, "/gp1_fit.RData"))

plot.new()
par(mfrow=c(1,1))

vals <- fit$par[names(fit$par) == "vals"]
# beta0 <- fit$par[names(fit$par) == "beta0"]

# plot(rate_use, exp(beta0 + vals))
plot(rate_use, exp(vals)[locs_index])
abline(0, 1)
# plot(rate_use)

##make some predictions

n_pred <- dim(cov_pred)[1]
sigma_11 <- matrix(NA, n_pred, n_pred)

for(i in 1:n_pred){
  for(j in 1:n_pred){
    sigma_11[i, j] <- exp(- sum((cov_pred[i, ] - cov_pred[j, ])^2) / sigma)
  }
}

sigma_12 <- matrix(NA, n_pred, n_points_u)

for(i in 1:n_pred){
  for(j in 1:n_points_u){
    sigma_12[i, j] <- exp(- sum((cov_pred[i, ] - cov_mat[j, ])^2) / sigma)
  }
}

sigma_21 <- t(sigma_12)

mu_new <- sigma_12 %*% solve(K, vals)

plot(exp(mu_new), Y_pred[, time_point_use] / pops_pred)

print(as.vector(cor(exp(mu_new), Y_pred[, time_point_use] / pops_pred, use="complete")))
