# #set working directory
# setwd(paste0(Sys.getenv("HOME"), "/v4_MDG/GP"))
# 
# #get z drive path
# load("../zdrive_path.RData")

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

tryCatch(dyn.unload(dynlib("gp_rff_t2")),
         error = function(e) print(e))
compile("gp_rff_t2.cpp")
dyn.load(dynlib("gp_rff_t2"))


source("model_fit_function.R")


# ##just do all covs for now
# time_iters <- 1:13
# n_months_fit <- 12
time_iters <- 1:2
n_months_fit <- 12


####try and discretise space to have fewer points in GP
##how many points with just two time points? this seems to go okay
##with fixed hyperparameters...
# time_iter <- 1
# n_months_fit <- 2
# cov_mat_list <- all_covs$covmats_normalised
# time_index_fit <- (time_iter - 1) + 1:n_months_fit
# cov_mat_fit <- cov_mat_list[time_index_fit]
# cov_mat_all <- do.call("rbind", cov_mat_fit)
# print("dimension of full covariate matrix (2 months of data): ")
# print(dim(cov_mat_all))
# #get unique elmements
# cov_mat <- unique(cov_mat_all)
# print("dimension of unique cov mat: ")
# print(dim(cov_mat))
## ~ 1000 points

# 
##now with 12 months
n_months_fit <- 12
cov_mat_list <- all_covs$covmats_normalised
time_index_fit <- (time_iter - 1) + 1:n_months_fit
cov_mat_fit <- cov_mat_list[time_index_fit]
cov_mat_all <- do.call("rbind", cov_mat_fit)
print("dimension of full covariate matrix (12 months of data): ")
print(dim(cov_mat_all))
#get unique elmements
cov_mat <- unique(cov_mat_all)
print("dimension of unique cov mat: ")
print(dim(cov_mat))

## ~ 6000 points - makes sense

cov_mat_2d <- cov_mat[, c(1, 5)]
n_dim <- dim(cov_mat_2d)[2]
n_pts <- dim(cov_mat)[1]
plot(cov_mat[, 1], cov_mat[, 2])
# min_vec <- c()
# max_vec <- c()
# n_disc <- 10

min_dist <- 3

cov_mat_out <- cov_mat_2d


##start with convex hull
ptm <- proc.time()
hull_index <- unique(as.vector(convhulln(cov_mat_2d)))
print("finding conxed hull")
print(proc.time() - ptm)
index_keep <- hull_index
index_test <- setdiff(1:n_pts, index_keep)

##randomly add 500 points to set
n_add <- 400
index_add <- sample(index_test, n_add)

index_keep <- c(index_keep, index_add)
index_test <- setdiff(index_test, index_add)

library(Rcpp)
# sourceCpp("mesh_function.cpp")
# ptm <- Sys.time()
# index_keep_cpp <- create_mesh_cpp(cov_mat_2d, index_keep-1,
#                 index_test-1, min_dist)
# print(paste0("cpp function ", Sys.time() - ptm))

sourceCpp("mesh_function_fast.cpp")
covs_use <- 1:28
ptm <- Sys.time()
index_keep_cpp <- create_mesh_cpp(cov_mat[, covs_use], index_keep-1,
                                  index_test-1, min_dist, TRUE)
print(paste0("cpp function ", Sys.time() - ptm))

# 
# d <- delaunayn(cov_mat[index_keep_cpp, 1:2])


# 
# i <- 2
# ptm <- Sys.time()
# while(i < dim(cov_mat_out)[1]){
#   cov_keep <- cov_mat_out[index_keep, , drop=FALSE]
#   cov_rest <- cov_mat_out[index_test, , drop=FALSE]
#   
#   # print(i)
#   # print(length(index_test))
#   # cov_keep <- cov_mat_out[1:i, ]
#   # cov_rest <- cov_mat_out[-(1:i), , drop=FALSE]
#   
#   #check if any of the points are less than min dist away from existing points
#   min_dists <- apply(cov_rest, 1, function(row_test){
#     min(apply(cov_keep, 1, function(row){
#       sqrt(sum((row - row_test)^2))
#       }
#       ), na.rm=T)
#   }
#   )
#   
#   keep_i <- which(min_dists > min_dist)
#   # plot(cov_keep)
#   # plot(cov_mat_2d)
#   # points(cov_rest[keep_i, ], col="green")
#   # points(cov_keep, col="red")
#   if(length(keep_i) == 0){
#     # cov_mat_out <- cov_keep
#     break
#   }else if(length(keep_i) == 1){
#     index_keep <- c(index_keep, index_test[keep_i])
#     break
#   }
#   
#   index_keep <- c(index_keep, index_test[keep_i[1]])
#   index_test <- index_test[keep_i[-1]]
#   # cov_rest <- cov_rest[keep_i, ]
#   # cov_mat_out <- rbind(cov_keep, cov_rest)
# 
#   i <- i + 1
# }
# print(Sys.time() - ptm)
# plot(cov_mat_2d[, 1:2])
# points(cov_mat_out[index_keep, 1:2], col="red", pch=19)

plot(cov_mat[, 1:2])
points(cov_mat[index_keep_cpp, 1:2], col="red", pch=19)

print(length(index_keep_cpp))

ptm <- proc.time()

# d <- delaunayn(cov_mat[index_keep_cpp, covs_use])
# print("time to triangulate")
# print(proc.time() - ptm)




