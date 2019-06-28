# #set working directory
# setwd(paste0(Sys.getenv("HOME"), "/v4_MDG/Model_1_outbreaks"))
# 
# #get z drive path
# load("../zdrive_path.RData")


path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")


setwd(path_input)
print(getwd())

beta_scales <- c(0.1, 0.075, 0.05, 0.025, 0.01, 0.001, 0.0001)
N_betas <- length(beta_scales)
# 
# load("outbreaks.RData")

for(beta_i in 5:2){

  print(paste0("beta_test_", beta_i, ".RData"))
  save(list = c("beta"),
       file = paste0("beta_test_", beta_i, ".RData"))

}


for(beta_i in 5:2){
  
  print(paste0(path_output, "/beta_test_", beta_i, ".RData"))
  save(list = c("beta"),
       file = paste0(path_output, "/beta_test_", beta_i, ".RData"))
  
}




