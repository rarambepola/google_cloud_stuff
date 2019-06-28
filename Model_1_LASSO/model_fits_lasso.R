# #set working directory
# setwd(paste0(Sys.getenv("HOME"), "/v4_MDG/Model_1"))
# 
# #get z drive path
# load("../zdrive_path.RData")

path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")

library(TMB)
library(INLA)

rmse_vec <- function(x, y, na.rm=T){
  return(sqrt(mean((x-y)^2, na.rm=na.rm)))
}

# covs_use <- c(1:4, 20)
# covs_use <- 1:28
# covs_use <- c(5:8, 20)
#add some fake covariates

load(paste0(path_input, "/covariates.RData"))

n_hf_fit <- dim(all_covs$covmats_normalised[[1]])[1]
n_hf_pred <- dim(pred_covs$covmats_normalised[[1]])[1]

models_silent <- TRUE

n_months_pred <- 24

# all_cov_cor <- c()
# causal_cov_cor <- c()
# lasso_cov_cor <- c()
subset_list <- list(1:28,
                    c(1:4, 20),
                    c(5:8, 20)
                    )
N_subsets <- length(subset_list)
subset_names <- c("all", "causal", "lasso")

residual_cors <- list()
all_cors <- list()


rowSds <- function(mat, na.rm=T){
  return(apply(mat, 1, sd, na.rm=T))
}

standardise_mat <- function(mat, na.rm=T){
  return((mat - rowMeans(mat, na.rm=na.rm)) / rowSds(mat, na.rm=na.rm))
}


#make mesh and spde object
mesh <- inla.mesh.2d(loc = coords_fit, cutoff = 0.1, max.edge = 3) 
plot(mesh)
points(coords_fit, col="red")
spde <- (inla.spde2.matern(mesh=mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
A_fit <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coords_fit))
A_pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coords_pred))
n_s <- nrow(spde$M0)
mesh_coords <- mesh$loc[, 1:2]

tryCatch(dyn.unload(dynlib("new_model_spatial_field_lasso")),
         error = function(e) print(e))
compile(paste0(path_input, "/new_model_spatial_field_lasso.cpp"))
dyn.load(dynlib("new_model_spatial_field_lasso"))

source(paste0(path_input, "/model_fit_function.R"))


##just do all covs for now
time_iters <- 1:10
n_months_fit <- 12
 

fit_rep_list_list <- list()

betas <- c(0.1, 0.075, 0.05, 0.025, 0.01, 0.001, 0.0001)
N_betas <- length(betas)

for(beta_i in 1:N_betas){
  fit_rep_list <- fit_model_lasso(time_iters,
                                   n_months_fit,
                                   Y_fit,
                                   pops_fit,
                                   A_fit,
                                   spde,
                                   n_s,
                                   all_covs$covmats_normalised,
                                   beta_scale = betas[beta_i],
                                   subset_list[[1]],
                                   iter.max=500,
                                   eval.max=500)
  
  save(list = c("fit_rep_list"),
       file = paste0(path_output, "/lasso_step_one_beta_", beta_i, ".RData"))
}


