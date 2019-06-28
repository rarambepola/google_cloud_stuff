# #set working directory
# setwd(paste0(Sys.getenv("HOME"), "/v4_MDG/Model_1"))

#get z drive path
# load("../zdrive_path.RData")

path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")


print(path_output)
print(path_input)
print(paste0(path_output, "/path_output"))