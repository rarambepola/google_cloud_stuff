path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")
setwd(path_input)
print(getwd())


load("covariates_all.RData")
# load("../../v4_MDG_big_files/all_covs_fit.RData")
load("fit.RData")
load("fake_data.RData")

n_hf_fit <- dim(all_covs$covmats_normalised[[1]])[1]
n_hf_pred <- dim(pred_covs$covmats_normalised[[1]])[1]

##combine covs
covs_all_list <- list()
covs_all_list$covmats_normalised <- list()
for(i in 1:length(all_covs$covmats_normalised)){
  covs_all_list$covmats_normalised[[i]] <- rbind(fit_cov_mats_use[[i]],
                                                 pred_cov_mats_use[[i]])
}

subset_list <- list(1:28,
                    c(1:4, 20),
                    c(5:8, 20)
)

##make predictions
pred_cov_list <- pred_covs$covmats_normalised

source("model_pred_functions2.R")

n_months_fit <- 12
n_months_pred <- 24
N_iters <- 13
all_outputs <- forward_predict_full(1:N_iters,
                                    fit_rep_list[[1]],
                                    covs_all_list$covmats_normalised,
                                    fit_rep_list[[2]],
                                    12,
                                    n_months_pred,
                                    rbind(A_fit, A_pred),
                                    rbind(apply(fit_data_fake, 2, "/", pops_fit),
                                          apply(pred_data_fake, 2, "/", pops_pred)),
                                    n_hf_fit,
                                    n_hf_pred,
                                    1:3)


save(list = c("all_outputs"),
     file = paste0(path_output, "/final_outputs.RData")
)
