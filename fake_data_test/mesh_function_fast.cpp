#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <math.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;

double norm(NumericVector v){
  int n=v.size();
  double v_norm = 0;
  
  for(int i=0; i<n; i++){
    v_norm += v(i) * v(i);
  }
  
  return(sqrt(v_norm));
}


// [[Rcpp::export]]
IntegerVector create_mesh_cpp(NumericMatrix cov_mat, IntegerVector index_keep,
                             IntegerVector index_test, double dist_cutoff, 
                             bool print_output=false){
  int n_keep = index_keep.size();
  int n_test = index_test.size();
  
  if(print_output){
    std::cout << "printing output" << "\n";
  }
  
  // (cov_mat.row(0) - cov_mat.row(1));
  

 
  // int stop_i = 10;
  // int count_i = 0;
  // while((n_test > 1) & (count_i < stop_i)){
  while(n_test > 1){
    double dist=0;
    IntegerVector index_keep_new = index_keep;
    
    //don't use pushback for index test - it's slow?
    IntegerVector index_test_new(n_test);
    int index_test_size = 0;
    n_keep = index_keep.size();
    bool found_one = FALSE;
    for(int i=0; i<n_test; i++){
      for(int j=0; j<n_keep; j++){
        dist = norm(cov_mat.row(index_test[i]) - cov_mat.row(index_keep_new[j]));
        if(dist < dist_cutoff){
          dist = -1;
          break;
        }
      }
      if(dist != -1){
        if(found_one){
          index_test_new[index_test_size] = index_test[i];
          index_test_size += 1;
        }else{
          index_keep_new.push_back(index_test[i]);
          found_one = TRUE;
        }
      }
    }
    
    IntegerVector index_test_new2(index_test_size);
    
    for(int i=0; i<index_test_size; i++){
      index_test_new2[i] = index_test_new[i];
    }
    index_test_new = index_test_new2;
    
    index_test = index_test_new;
    index_keep = index_keep_new;
    n_test = index_test.size();
    if(print_output) std::cout << n_test << "\n";
    
    if(n_test == 1){
      index_keep.push_back(index_test[0]);
    }
  }
  
  return(index_keep + 1);
}
