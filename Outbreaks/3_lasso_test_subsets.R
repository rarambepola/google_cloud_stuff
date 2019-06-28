# #set working directory
# setwd(paste0(Sys.getenv("HOME"), "/v4_MDG/Model_1_outbreaks"))

#get z drive path
# load("../zdrive_path.RData")
path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")
library(TMB)
library(INLA)

setwd(path_input)
print(getwd())



path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")


# load("../covariates.RData")
load("covariates_all.RData")
load("outbreaks.RData")


n_hf_fit <- dim(all_covs$covmats_normalised[[1]])[1]
n_hf_pred <- dim(pred_covs$covmats_normalised[[1]])[1]

models_silent <- TRUE

n_months_pred <- 24

# all_cov_cor <- c()
# causal_cov_cor <- c()
# lasso_cov_cor <- c()
# subset_list <- list(1:28,
#                     c(1:4, 20),
#                     c(5:8, 20)
#                     )
# N_subsets <- length(subset_list)
# subset_names <- c("all", "causal", "lasso")

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

tryCatch(dyn.unload(dynlib("outbreaks")),
         error = function(e) print(e))
compile("outbreaks.cpp")
dyn.load(dynlib("outbreaks"))

source("model_fit_function.R")
load("lasso_subsets_to_test.RData")
source("model_pred_functions2.R")

time_iters <- 1:13
n_months_fit <- 12

N_iter <- 5

N_pts_total <- n_hf_fit
N_train_pts <- 250

test_train_split <- replicate(N_iter, sample.int(N_pts_total, N_train_pts))

lasso_loss <- list()


for(time_iter in time_iters){
# for(time_iter in 1:2){
  causal_subsets <- lasso_subsets_to_test[[time_iter]]
  N_subsets <- length(causal_subsets)
  iter_loss <- list()

  for(subset_i in 1:N_subsets){
  # for(subset_i in 1:1){
    subset <- causal_subsets[[subset_i]]
    subset_loss <- c()
    for(iter_i in 1:N_iter){
    # for(iter_i in 1:2){
      train_ind <- test_train_split[, iter_i]
      test_ind <- setdiff(1:N_pts_total, train_ind)
      #make mesh and spde object
      mesh <- inla.mesh.2d(loc = coords_fit[train_ind, ], cutoff = 0.1, max.edge = 3)
      # plot(mesh)
      # points(coords_fit, col="red")
      spde <- (inla.spde2.matern(mesh=mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
      A_fit <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coords_fit[train_ind, ]))
      A_pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coords_fit[test_ind, ]))
      n_s <- nrow(spde$M0)
      mesh_coords <- mesh$loc[, 1:2]

      #fix covariate matrices
      covs_fit <- lapply(all_covs$covmats_normalised, function(mat) mat[train_ind, ])
      covs_pred <- lapply(all_covs$covmats_normalised, function(mat) mat[test_ind, ])

      fit_rep <- fit_model_final(time_iter,
                                      6,
                                      outbreaks_fit_full[train_ind, ],
                                      A_fit,
                                      spde,
                                      n_s,
                                      covs_fit,
                                      subset,
                                      iter.max=500,
                                      eval.max=500,
                                      get_rep = TRUE)

      #make predictions
      model_pred <- forward_predict_full(time_iter,
                                         fit_rep[[1]],
                                         covs_pred,
                                         fit_rep[[2]],
                                         6,
                                         6,
                                         A_pred,
                                         outbreaks_fit_full[test_ind, ],
                                         subset_use = subset)
      

      subset_loss[iter_i] <- model_pred$loss[1]

    }
    iter_loss[[subset_i]] <- subset_loss
  }
  # for(subset_i in 1:2){



  lasso_loss[[time_iter]] <- iter_loss

  save(list = c("lasso_loss",
                "mesh",
                "spde",
                "A_fit",
                "A_pred"),
       file = paste0(path_outputs, "/lasso_subset_test_fits_iter_", time_iter, ".RData"))
}

# save(list = c("lasso_loss",
#               "mesh",
#               "spde",
#               "A_fit",
#               "A_pred"),
#      file = "lasso_subset_test_fits.RData")

