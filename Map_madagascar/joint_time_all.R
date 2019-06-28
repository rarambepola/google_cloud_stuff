# #set working directory
# setwd(paste0(Sys.getenv("HOME"), "/Map_madagascar"))
# 
# #get z drive path
# load("zdrive_path.RData")
path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")
setwd(path_input)
print(getwd())

print(paste0(path_output, "/joint_time_all_fits2.RData"))

library(TMB)
library(INLA)
library(raster)
library(ggplot2)
library(stringr)

tryCatch(dyn.unload(dynlib("joint_time_model_all")),
         error = function(e) print(e))
compile("joint_time_model_all.cpp")
dyn.load(dynlib("joint_time_model_all"))

#
# ####pr####
# models_silent <- FALSE
#
# reference_raster <- raster("C:/Users/scro3122/Documents/Map_madagascar/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif")
# mdg_extent <- extent(reference_raster)
#
#
# load("pr_data.RData")
#
# #keep variables we care about
# pr_df <- mdg_pr_new[, c("latitude", "longitude", "pf_pos", "examined", "year_start",
#                         "year_end", "month_start", "month_end")]
# print(sum(is.na(pr_df)))
# coords_pr <- pr_df[, c("longitude", "latitude")]
#
#
# ##for now just keep 2016 vals
# pr_df <- pr_df[pr_df$year_start > 2015, ]
#
# #year of first survey - should be 2011
# # year_start <- min(pr_df$year_start)
# # month_start <- min(pr_df$year_start)
# year_start <- 2015
#
# # month_all <- 1:24
#
# pr_df$"month" <- (pr_df$year_start - year_start) * 12 + pr_df$month_start
# # pr_df$month <- 3 + sample.int(10, length(pr_df$month), replace=TRUE)
#
# # month_clip <- seq(3, 70, by=6)
# # month_clip <- c(10, 40)
#
#
# # month_all <- seq(3, 70, by=6)
#
# # pr_df$"month" <- sapply(pr_df$month, function(m) month_all[which.min(abs(m - month_all))[1]])
#
# ##find months that exist
# month_exist <- unique(pr_df$month)
#
#
# month_all <- c(15, month_exist)
# month_inc_all <- match(month_all + 24, 1:48)
#
# # ##limit times for now
# # n_month_use <- 1
# # month_exist <- month_exist[1:n_month_use]
# year_exist <- year_start + floor(month_exist / 12)
# month_of_year_exist <- ((month_exist - 1) %% 12) + 1
#
#
# #change into half years
# # month_of_year_exist <- sapply(month_of_year_exist, function(month) if(month <6) return(3) else return(9))
# year_month_exist <- paste0(year_exist, ".", str_pad(month_of_year_exist, 2, pad="0"))
# year_month_exist <- unique(year_month_exist)
# # month_exist <- (year_exist - year_start)*12 + month_of_year_exist
# month_exist <- unique(month_exist)
# # n_months <- length(unique(month_exist))
# n_months <- length(month_exist)
# # month_all <- 1:max(month_exist)
# # month_all <- c(1, month_exist)
# n_months_all <- length(month_all)
#
#
# cov_path_static <- c("Z:/mastergrids/Other_Global_Covariates/Accessibility/Weiss/accessibility_to_cities_2015_v1.0.tif",
#                      "Z:/mastergrids/Other_Global_Covariates/Aridity/5km/Aridity_Index.5k.MEAN.tif",
#                      "Z:/mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif"
# )
#
# cov_path_start <- c("Z:/mastergrids/Other_Global_Covariates/Rainfall/CHIRPS/5km/Monthly/chirps-v2-0.",
#                     "Z:/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Monthly/LST_Day_v6.",
#                     "Z:/mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/EVI_v6/5km/Monthly/EVI_v6.")
#
# cov_path_start_all <- c(cov_path_static,
#                         cov_path_start)
#
# cov_path_end <- c(".sum.5km.NN.tif",
#                   ".mean.5km.mean.tif",
#                   ".mean.5km.mean.tif")
#
# cov_path_end_all <- c(cov_path_static,
#                       cov_path_end)
#
# n_covs_static <- length(cov_path_static)
# n_covs_dynamic <- length(cov_path_start)
# n_covs <- n_covs_static + n_covs_dynamic
#
# cov_type <- c(rep("static", n_covs_static),
#               rep("dynamic", n_covs_dynamic))
#
# ##for the months that exists get prevelance data points
# Y_pr <- list()
# N_pr <- list()
# coords_pr_list <- list()
# cov_mat_list <- list()
# for(i in 1:n_months){
#   print(paste0("month ", i, " of ", n_months))
#   print(year_month_exist[i])
#   Y_pr[[i]] <- pr_df$pf_pos[pr_df$month == month_exist[i]]
#   N_pr[[i]] <- pr_df$examined[pr_df$month == month_exist[i]]
#   coords_pr_list[[i]] <- cbind(pr_df$longitude[pr_df$month == month_exist[i]],
#                                pr_df$latitude[pr_df$month == month_exist[i]]
#   )
#   n_locs <- length(Y_pr[[i]])
#
#   cov_mat <- matrix(NA, n_locs, n_covs)
#   for(j in 1:n_covs){
#     if(cov_type[j] == "static"){
#       print(cov_path_start_all[j])
#       r <- raster(cov_path_start_all[j])
#     }else{
#       print(paste0(cov_path_start_all[j],
#                    year_month_exist[i],
#                    cov_path_end_all[j]
#       ))
#       r <- raster(paste0(cov_path_start_all[j],
#                          year_month_exist[i],
#                          cov_path_end_all[j]
#       ))
#     }
#     r <- crop(r, mdg_extent)
#     cov_mat[, j] <- extract(r, coords_pr_list[[i]])
#   }
#   cov_mat_list[[i]] <- cov_mat
# }
#
# ##normalise
# cov_mat_all <- do.call("rbind", cov_mat_list)
# cov_means <- colMeans(cov_mat_all, na.rm=T)
# cov_sds <- apply(cov_mat_all, 2, sd, na.rm=T)
#
# plot(cov_mat_list[[1]][, 1])
# cov_mat_list_n <- lapply(cov_mat_list, function(mat) t((t(mat) - cov_means) / cov_sds))
# plot(cov_mat_list_n[[1]][, 1])
#
#
# ##replace NAs with average
# cov_mat_final <- list()
# for(i in 1:n_months){
#   cov_mat <- cov_mat_list_n[[i]]
#   cov_mat_out <- cov_mat
#   for(j in 1:n_covs){
#     nas <- which(is.na(cov_mat[, j]))
#     if(length(nas) > 0){
#       cov_mat_out[nas, j] <- mean(cov_mat_out[, j], na.rm=T)
#     }
#   }
#   cov_mat_final[[i]] <- cov_mat_out
# }
#
#
# ##get covariate values
# #
# #
# # coords_pr_11 <- coords_pr[pr_df$year_start == 2011, ]
# # coords_pr_16 <- coords_pr[pr_df$year_start == 2016, ]
#
# #create mesh
# mesh <- inla.mesh.2d(loc = coords_pr, cutoff = 0.5, max.edge = 3)
# plot(mesh, asp=1)
# points(coords_pr, col="red")
# spde <- (inla.spde2.matern(mesh=mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
# A <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coords_pr))
# n_s <- nrow(spde$M0)
# mesh_coords <- mesh$loc[, 1:2]
#
# A_pixel <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coordinates(reference_raster)))
#
# A_pr <- list()
# for(i in 1:n_months){
#   A_pr[[i]] <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coords_pr_list[[i]]))
# }
#
# # n_month_use <- 3
# # X_pr <- list(cov_mat_11, cov_mat_16)
# # Y_pr <- list(pr_df$pf_pos[pr_df$year_start == 2011], pr_df$pf_pos[pr_df$year_start == 2016])
# # N_pr <- cbind(pr_df$examined, pr_df$examined)
# # N_pr <- list(pr_df$examined[pr_df$year_start == 2011], pr_df$examined[pr_df$year_start == 2016])
# # A_pr <- list(A_11, A_16)
#
# N_covs <- n_covs
# # month_all <- month_exist
# n_months_all <- length(month_all)
#
#
#
# #####incidence
# ##incidence
# # load("../MDG_v3/Covariates/covs_add_cases.RData")
# # load("../MDG_v3/Covariates/hf_dat_clean.RData")
#
# load("../v4_MDG_big_files/covariates_all.RData")
#
# ##get incidence covariate values
# year_month_inc_all <- paste0(rep(2013:2016, each=12), ".", str_pad(rep(1:12, 4), 2, pad="0"))
# year_month_inc <- year_month_inc_all[month_inc_all]
# n_months_inc <- length(year_month_inc)
# n_hf_fit <- dim(coords_fit)[1]
# cov_mat_list_inc <- list()
# for(i in 1:n_months_inc){
#   print(paste0("month ", i, " of ", n_months_inc))
#   print(year_month_inc[i])
#
#   cov_mat <- matrix(NA, n_hf_fit, n_covs)
#   for(j in 1:n_covs){
#     if(cov_type[j] == "static"){
#       print(cov_path_start_all[j])
#       r <- raster(cov_path_start_all[j])
#     }else{
#       print(paste0(cov_path_start_all[j],
#                    year_month_inc[i],
#                    cov_path_end_all[j]
#       ))
#       r <- raster(paste0(cov_path_start_all[j],
#                          year_month_inc[i],
#                          cov_path_end_all[j]
#       ))
#     }
#     r <- crop(r, mdg_extent)
#     cov_mat[, j] <- extract(r, coords_fit)
#   }
#   cov_mat_list_inc[[i]] <- cov_mat
# }
#
# A_inc <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coordinates(coords_fit)))
#
#
# coords_pixel <- coordinates(reference_raster)
# n_pixel <- dim(coords_pixel)[1]
# ##make all covariates matrix
# cov_mat_list_pixel <- list()
# for(i in 1:n_months_inc){
#   print(paste0("month ", i, " of ", n_months_inc))
#
#   cov_mat <- matrix(NA, n_pixel, n_covs)
#   for(j in 1:n_covs){
#     if(cov_type[j] == "static"){
#       print(cov_path_start_all[j])
#       r <- raster(cov_path_start_all[j])
#     }else{
#       print(paste0(cov_path_start_all[j],
#                    year_month_inc[i],
#                    cov_path_end_all[j]
#       ))
#       r <- raster(paste0(cov_path_start_all[j],
#                          year_month_inc[i],
#                          cov_path_end_all[j]
#       ))
#     }
#     r <- crop(r, mdg_extent)
#     cov_mat[, j] <- extract(r, coords_pixel)
#   }
#   cov_mat_list_pixel[[i]] <- cov_mat
# }
#
# A_inc <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coordinates(coords_fit)))
# A_pixel <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coordinates(coords_pixel)))
#
#
# ##normalise
# # plot(cov_mat_list[[1]][, 1])
# cov_mat_list_inc_n <- lapply(cov_mat_list_inc, function(mat) t((t(mat) - cov_means) / cov_sds))
# cov_mat_list_pixel_n <- lapply(cov_mat_list_pixel, function(mat) t((t(mat) - cov_means) / cov_sds))
# # plot(cov_mat_list_n[[1]][, 1])
#
#
# ##replace NAs with average
# cov_mat_final_inc <- list()
# for(i in 1:n_months_inc){
#   cov_mat <- cov_mat_list_inc_n[[i]]
#   cov_mat_out <- cov_mat
#   for(j in 1:n_covs){
#     nas <- which(is.na(cov_mat[, j]))
#     if(length(nas) > 0){
#       cov_mat_out[nas, j] <- mean(cov_mat_out[, j], na.rm=T)
#     }
#   }
#   cov_mat_final_inc[[i]] <- cov_mat_out
# }
#
#
# ##remove rows with NAs
# cov_mat_final_pixel <- list()
# rows_kept <- list()
# for(i in 1:n_months_all){
#   cov_mat <- cov_mat_list_pixel_n[[i]]
#   keep_row <- which(!is.na(rowSums(cov_mat)))
#   rows_kept[[i]] <- keep_row
#   cov_mat_out <- cov_mat[keep_row, ]
#   # for(j in 1:n_covs){
#   #   nas <- which(is.na(cov_mat[, j]))
#   #   if(length(nas) > 0){
#   #     cov_mat_out[nas, j] <- mean(cov_mat_out[, j], na.rm=T)
#   #   }
#   # }
#   cov_mat_final_pixel[[i]] <- cov_mat_out
# }
#
# ##fit model?

# save.image("joint_time_all.RData")

load("joint_time_all2.RData")

m <- MakeADFun(
  data = list(X_pr=cov_mat_final,
              Y_pr=Y_pr,
              N_pr=N_pr,
              A_pr=A_pr,
              spde=spde,
              times_pr=month_exist,
              times_all=month_all,
              index_times_pr = match(month_exist, month_all) - 1,

              X_inc=cov_mat_final_inc,
              Y_inc=Y_fit[, month_inc_all],
              pops_inc=pops_fit,
              A_inc=A_inc,
              times_inc=month_inc_all,
              X_pixel=cov_mat_final_pixel,
              A_pixel=A_pixel
  ),
  parameters = list(beta0_pr=runif(1, -1, 1),
                    beta_pr=rep(0, n_covs),
                    S_pr=matrix(0, n_months_all, n_s),
                    log_kappa_pr=0.5,
                    log_tau_pr=0.0,
                    log_scale_pr=0.0,

                    beta0_inc=runif(1, -1, 1),
                    beta_inc=rep(0, n_covs),
                    S_inc=matrix(0, n_months_all, n_s),
                    log_kappa_inc=0.5,
                    log_tau_inc=0.0,
                    log_scale_inc=0.0,
                    inc_pr_slope=1.0
                    ),
  random = c("S_pr", "S_inc"),
  DLL = "joint_time_model_all",
  silent=models_silent
)
###
##fit model
ptm <- proc.time()
fit <- nlminb(m$par, m$fn, m$gr, control=list(iter.max=300,eval.max=300))
rep <- sdreport(m)
print(proc.time() - ptm)

# fit <- "test"
# n_s <- 10000
# rep <- matrix(runif(n_s^2), n_s, n_s)

save(list = c("fit", "rep"),
     file=paste0(path_output, "/joint_time_all_fits2.RData"))



# 
# ##predictions
# beta0_pr <- fit$par[names(fit$par) == "beta0_pr"]
# beta0_inc <- fit$par[names(fit$par) == "beta0_inc"]
# beta_pr <- fit$par[names(fit$par) == "beta_pr"]
# beta_inc <- fit$par[names(fit$par) == "beta_inc"]
# 
# S_inc <- matrix(rep$par.random[names(rep$par.random) == "S_inc"], ncol=n_s)
# S_pr <- matrix(rep$par.random[names(rep$par.random) == "S_pr"], ncol=n_s)
# 
# ##make maps
# 
# plot.new()
# par(mfrow=c(2,4))
# for(month_use in 1:4){
#   
#   inc_raster <- reference_raster
#   pr_raster <- reference_raster
#   
#   logit_prob <- as.vector(beta0_pr + (cov_mat_final_pixel[[month_use]] %*% beta_pr) + (A_pixel[rows_kept[[month_use]], ] %*% S_pr[month_use, ]))
#   log_rate <- as.vector(beta0_inc + (cov_mat_final_pixel[[month_use]] %*% beta_inc) + (A_pixel[rows_kept[[month_use]], ] %*% S_inc[month_use, ]))
#   
#   values(inc_raster) <- NA
#   values(inc_raster)[rows_kept[[month_use]]] <- exp(log_rate)
#   values(pr_raster) <- NA
#   values(pr_raster)[rows_kept[[month_use]]] <- exp(logit_prob) / (1 + exp(logit_prob))
#   
#   
#   plot(inc_raster, main=paste0("incidence ", month_use), zlim=c(0, 0.035))
#   plot(pr_raster, main=paste0("prevalence ", month_use), zlim=c(0, 0.3))
# }
# par(mfrow=c(1,1))
# 






