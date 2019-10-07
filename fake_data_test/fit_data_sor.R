#check whether on google cloud or not
if(Sys.getenv("HOME") == "/root") google_cloud <- TRUE else google_cloud <- FALSE

parallel_on <- TRUE
print(paste0("parallel: ", parallel_on))
#set up directories
if(google_cloud){
  print("fit_data_gp.R")
  path_output <- Sys.getenv("path_output")
  path_input <- Sys.getenv("path_input")
  setwd(path_input)
  print(getwd())
  covariates_folder <- NULL
  output_folder <- paste0(path_output, "/")
}else{
  setwd(paste0(Sys.getenv("HOME"), "/google_cloud_stuff/fake_data_test"))
  covariates_folder <- "../../v4_MDG_big_files/"
  output_folder <- NULL
}


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
load(paste0(covariates_folder, "covariates_all.RData"))

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

tryCatch(dyn.unload(dynlib("sor")),
         error = function(e) print(e))
compile("sor.cpp")
dyn.load(dynlib("sor"))


source("model_fit_function_gp_sor.R")


# ##just do all covs for now
# time_iters <- 1:13
# n_months_fit <- 12
time_iters <- 1:3
n_months_fit <- 12



library(Rcpp)

ptm <- proc.time()

load("fake_data.RData")

fit_rep_list <- fit_model_mesh(time_iters,
                          n_months_fit,
                          fit_data_fake,
                          pops_fit,
                          A_fit,
                          spde,
                          n_s,
                          fit_cov_mats_use,
                          subset_use=1:3,
                          iter.max=500,
                          eval.max=500,
                          get_rep = TRUE,
                          getJointPrecision = FALSE,
                          models_silent=FALSE,
                          max_dist=2.5,
                          parallel = parallel_on
                          sigma=sigma,
                          lambda=lambda)

print("time to fit: ")
print(proc.time() - ptm)

path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")
save(list = c("fit_rep_list",
              "mesh",
              "spde",
              "A_fit",
              "A_pred"),
     file = paste0(output_folder, "fit_gp_sor.RData"))
