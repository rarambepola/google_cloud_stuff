#check whether on google cloud or not
if(Sys.getenv("HOME") == "/root") google_cloud <- TRUE else google_cloud <- FALSE

#set up directories
if(google_cloud){
  print("model_predictions_gp.R")
  path_output <- Sys.getenv("path_output")
  path_input <- Sys.getenv("path_input")
  setwd(path_input)
  print(getwd())
  covariates_folder <- NULL
  output_folder <- paste0(path_output, "/")
}else{
  setwd(paste0(Sys.getenv("HOME"), "/v4_MDG/GP"))
  covariates_folder <- "../../v4_MDG_big_files/"
  output_folder <- NULL
}

load(paste0(covariates_folder, "covariates_all.RData"))
# load("../../v4_MDG_big_files/all_covs_fit.RData")
# load("fit_gp_sor_1_0.01.RData")
sigma <- 100
lambda <- 0.1

print("mesh_fixed_par_pred_sor.R")
print(paste0("sigma: ", sigma))
print(paste0("lambda: ", lambda))

load(paste0("fit_gp_sor_", sigma, "_", lambda, ".RData"))
# load("fake_data.RData")

n_hf_fit <- dim(all_covs$covmats_normalised[[1]])[1]
n_hf_pred <- dim(pred_covs$covmats_normalised[[1]])[1]

##combine covs
covs_all_list <- list()
covs_all_list$covmats_normalised <- list()
for(i in 1:length(all_covs$covmats_normalised)){
  covs_all_list$covmats_normalised[[i]] <- rbind(all_covs$covmats_normalised[[i]],
                                                 pred_covs$covmats_normalised[[i]])
}

subset_list <- list(1:28,
                    c(1:4, 20),
                    c(5:8, 20)
)

##make predictions
pred_cov_list <- pred_covs$covmats_normalised

source("model_pred_functions2_gp_sor.R")

n_months_fit <- 12
n_months_pred <- 24
N_iters <- 13
# all_outputs <- forward_predict_full(1:N_iters,
#                                     fit_rep_list[[1]],
#                                     covs_all_list$covmats_normalised,
#                                     fit_rep_list[[2]],
#                                     12,
#                                     n_months_pred,
#                                     rbind(A_fit, A_pred),
#                                     rbind(apply(fit_data_fake, 2, "/", pops_fit),
#                                           apply(pred_data_fake, 2, "/", pops_pred)),
#                                     n_hf_fit,
#                                     n_hf_pred,S
#                                     1:3)

pred_list <- forward_predict_full(1:N_iters, 
                                  fit_rep_list[[1]],
                                  covs_all_list$covmats_normalised,
                                  all_covs$covmats_normalised,
                                  fit_rep_list[[2]],
                                  fit_rep_list[[3]],
                                  n_months_fit,
                                  n_months_pred,
                                  rbind(A_fit, A_pred),
                                  rbind(apply(Y_fit, 2, "/", pops_fit),
                                        apply(Y_pred, 2, "/", pops_pred)),
                                  n_hf_fit,
                                  n_hf_pred,
                                  A_fit,
                                  subset_list[[3]],
                                  sigma=sigma,
                                  lambda=lambda
)

save(list = c("pred_list"),
     file = paste0(output_folder, "pred_gp_sor_", sigma, "_", lambda, ".RData"))

# all_outputs <- forward_predict_full(1:N_iters, 
#                                   fit_rep_list[[1]],
#                                   covs_all_list$covmats_normalised,
#                                   fit_cov_mats_use,
#                                   fit_rep_list[[2]],
#                                   fit_rep_list[[3]],
#                                   n_months_fit,
#                                   12,
#                                   rbind(A_fit, A_pred),
#                                   rbind(apply(fit_data_fake, 2, "/", pops_fit),
#                                         apply(pred_data_fake, 2, "/", pops_pred)),
#                                   n_hf_fit,
#                                   n_hf_pred,
#                                   A_fit,
#                                   subset_use=1:3,
#                                   sigma=sigma,
#                                   lambda=lambda
# )


# save(list = c("all_outputs"),
#      file = paste0(path_output, "final_outputs_gp_sor.RData")
# )
