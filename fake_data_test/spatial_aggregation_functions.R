

#function to aggregate cases
case_aggregate_matrix <- function(case_matrix, district_index){
  # print(dim(case_matrix))
  # print(length(district_index))
  n_t <- dim(case_matrix)[2]
  unq_dists <- na.omit(unique(district_index))
  N_dists <- length(unq_dists)
  aggregate_cases <- matrix(NA, ncol=n_t, nrow=N_dists)

  for(i in 1:N_dists){
    unq_dist <- unq_dists[i]
    # print(length(colMeans(case_matrix[which(district_index == unq_dist), ,
    #                                   drop=FALSE],
    #                       na.rm=T)))
    # print(length(aggregate_cases[i, ]))
    # print(n_t)
    aggregate_cases[i, ] <- colMeans(case_matrix[which(district_index == unq_dist), ,
                                                 drop=FALSE], 
                                     na.rm=T)
    
  }
  return(aggregate_cases)
}