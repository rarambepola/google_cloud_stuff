path_output <- Sys.getenv("path_output")

x <- runif(10)

y <- 2*x

save(list = c("x", "y"),
     file = paste0(path_output, "xy.RData"))