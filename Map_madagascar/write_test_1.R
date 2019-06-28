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

n_s <- 1000000
rep <- matrix(runif(n_s*n_s), n_s, n_s)

save(list = c("rep"),
     file=paste0(path_output, "/write_test1.RData"))


