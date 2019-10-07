setwd(paste0(Sys.getenv("HOME"), "/google_cloud_stuff/fake_data_test"))
library(TMB)

N <- 100
N_test <- 100
x_train <- runif(N, max=5)
x_test <- runif(N_test, max=5)
y_train <- sin(x_train) + rnorm(N, sd=0.3)
y_test <- sin(x_test)

plot(x_train, y_train)


N_subset <- 10

index_subset <- sample.int(N_test, N_subset)

nearest_subset_point <- sapply(x_train, function(x_point) 
  which.min((x_point - x_train[index_subset])^2))

nearest_subset_point_test <- sapply(x_test, function(x_point) 
  which.min((x_point - x_train[index_subset])^2))

#check
i <- 1
points(x_train[i], y_train[i], col="red", pch=19)
points(x_train[index_subset[nearest_subset_point[i]]],
       y_train[index_subset[nearest_subset_point[i]]],
       col="green", pch=19)

points(x_train[index_subset[nearest_subset_point]],
       y_train[index_subset[nearest_subset_point]],
       col="blue")

#make K inv
sigma <- 1
lambda <- 0.01
K <- matrix(NA, N_subset, N_subset)
for(i in 1:N_subset){
  for(j in 1:N_subset){
    K[i, j] <- exp(-(x_train[index_subset[i]]
                     - x_train[index_subset[j]])^2/sigma)
  }
}
K_inv <- solve(K + lambda*diag(N_subset))



model_file_name <- "sod_test"
tryCatch(dyn.unload(dynlib(model_file_name)),
         error = function(e) print(e))
compile(paste0(model_file_name, ".cpp"))
dyn.load(dynlib(model_file_name))


m <- MakeADFun(
  data = list(Y=y_train,
              locs_index=nearest_subset_point - 1,
              K_inv=K_inv
  ),
  parameters = list(vals=runif(N_subset)
  ),
  DLL = "sod_test"
)


fit <- nlminb(m$par, m$fn, m$gr)

fit_vals <- fit$par
fit_vals_all <- fit_vals[nearest_subset_point]
# plot(fit_vals, y_train[index_subset])
# plot(fit_vals_all, y_train)
plot(x_train, fit_vals_all)
# plot(x_train[index_subset], fit_vals_all[index_subset])
points(x_train, y_train, col="red")
# points(x_train[index_subset], y_train[index_subset], col="blue")

plot(fit_vals[nearest_subset_point_test], y_test)
abline(0, 1)

##try subset of regressors
#need to make projection matrix
K_proj_a <- matrix(NA, ncol=N_subset, nrow=N)

for(i in 1:N){
  for(j in 1:N_subset){
    K_proj_a[i, j] <- exp(-(x_train[i]
                            - x_train[index_subset[j]])^2/sigma)
  }
}

K_proj <- K_proj_a %*% K_inv

K_proj_a_test <- matrix(NA, ncol=N_subset, nrow=N_test)

for(i in 1:N){
  for(j in 1:N_subset){
    K_proj_a_test[i, j] <- exp(-(x_test[i]
                            - x_train[index_subset[j]])^2/sigma)
  }
}
K_proj_test <- K_proj_a_test %*% K_inv


model_file_name <- "sor_test"
tryCatch(dyn.unload(dynlib(model_file_name)),
         error = function(e) print(e))
compile(paste0(model_file_name, ".cpp"))
dyn.load(dynlib(model_file_name))


m <- MakeADFun(
  data = list(Y=y_train,
              project_mat=K_proj,
              K_inv=K_inv
  ),
  parameters = list(vals=runif(N_subset)
  ),
  DLL = "sor_test"
)


fit <- nlminb(m$par, m$fn, m$gr)

fit_vals <- fit$par
fit_vals_all <- K_proj %*% fit_vals
plot(fit_vals, y_train[index_subset])
# plot(fit_vals_all, y_train)
plot(x_train, fit_vals_all)
# plot(x_train[index_subset], fit_vals_all[index_subset])
points(x_train, y_train, col="red")
# points(x_train[index_subset], y_train[index_subset], col="blue")

plot(K_proj_test %*% fit_vals, y_test)
abline(0, 1)
