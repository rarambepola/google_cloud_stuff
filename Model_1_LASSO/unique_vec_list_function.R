vec_list_sort <- function(vec_list){
  vec_lengths <- sapply(vec_list, length)
  vec_order <- order(vec_lengths)
  return(vec_list[vec_order])
}

unique_vec_list <- function(vec_list, sort=TRUE){
  out_list <- list()
  out_list[[1]] <- vec_list[[1]]
  for(i in 2:length(vec_list)){
    test_vec <- vec_list[[i]]
    already_in_list <- FALSE
    
    for(out_vec in out_list){
      if(length(out_vec) == length(test_vec)){
        if(all(sort(test_vec) == sort(out_vec))) already_in_list <- TRUE
      }
    }
    if(!already_in_list) out_list[[length(out_list) + 1]] <- vec_list[[i]]
  }
  
  if(sort){
    out_list <- out_list[order(sapply(out_list, length))]
  }
  return(out_list)
}