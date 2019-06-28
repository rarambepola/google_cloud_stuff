path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")
# setwd(paste0(Sys.getenv("HOME"), "/google_cloud"))

N <- 100
x <- runif(N, max=10)

y <- 2*x + rnorm(N)

if(!("TMB" %in% installed.packages())) install.packages("TMB")

library(TMB)

tryCatch(dyn.unload(dynlib("tmb_test")),
         error = function(e) print(e))
compile(paste0(path_input, "/tmb_test.cpp"))
dyn.load(dynlib("tmb_test"))

m <- MakeADFun(
  data = list(X=x,
              Y=y
  ),
  parameters = list(beta=runif(1)),
  DLL = "tmb_test"
)

ptm <- proc.time()
fit <- nlminb(m$par, m$fn, m$gr, control=list(iter.max=500,eval.max=500))
# rep <- sdreport(m, getJointPrecision = TRUE)
print(proc.time() - ptm)

print(fit)