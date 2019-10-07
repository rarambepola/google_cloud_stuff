path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")
setwd(path_input)
print(getwd())


library(TMB)
library(INLA)
library(geometry)

rmse_vec <- function(x, y, na.rm=T){
  return(sqrt(mean((x-y)^2, na.rm=na.rm)))
}

# covs_use <- c(1:4, 20)
# covs_use <- 1:28
# covs_use <- c(5:8, 20)
#add some fake covariates

# load("../covariates.RData")
load("covariates_all.RData")

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

tryCatch(dyn.unload(dynlib("gp_rff_t2_fixed")),
         error = function(e) print(e))
compile("gp_rff_t2_fixed.cpp")
dyn.load(dynlib("gp_rff_t2_fixed"))


source("model_fit_function_mesh_fixed_hypers_par.R")


# ##just do all covs for now
# time_iters <- 1:13
# n_months_fit <- 12
time_iters <- 1:13
n_months_fit <- 12



library(Rcpp)

ptm <- proc.time()

fit_rep_list <- fit_model_mesh(time_iters,
                          n_months_fit,
                          Y_fit,
                          pops_fit,
                          A_fit,
                          spde,
                          n_s,
                          all_covs$covmats_normalised,
                          subset_use=subset_list[[1]],
                          iter.max=500,
                          eval.max=500,
                          get_rep = TRUE,
                          getJointPrecision = FALSE,
                          models_silent=FALSE,
                          max_dist=2)

print("time to fit: ")
print(proc.time() - ptm)

path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")
save(list = c("fit_rep_list",
              "mesh",
              "spde",
              "A_fit",
              "A_pred"),
     file = paste0(path_output, "/all_covs_mesh_fit_fixed_hypers_par.RData"))

