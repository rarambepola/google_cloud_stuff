#include <TMB.hpp>
#include <math.h>
#include <stdio.h>

template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace Eigen;
  
  DATA_VECTOR(Y);
  DATA_MATRIX(project_mat);
  DATA_MATRIX(K_inv);
  
  PARAMETER_VECTOR(vals);
  
  Type f=0;
  int n = Y.size();
  Type sd = 0.1;
  
  vector<Type> vals_all = project_mat * vals;
  for(int j=0; j<n; j++){
    f -= dnorm(Y(j), vals_all(j), sd, true);
  }
  f += (vals * (K_inv * vals)).sum();

  return(f);
}
