path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")

if(!("FSIC" %in% installed.packages())){
  install.packages(paste0(path_input, "/FSIC_0.1.0.tar.gz"), repos=NULL, type="source")
}

packages_needed <- c("Rcpp", "foreach", "doParallel")
packages_not_installed <- packages_needed[!packages_needed %in% installed.packages()[, "Package"]]

if(length(packages_not_installed) > 0) install.packages(packages_not_installed)


source("runpc.R")
source("runpc_part2.R")

print("done")