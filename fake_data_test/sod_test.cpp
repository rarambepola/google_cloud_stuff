#include <TMB.hpp>
#include <math.h>
#include <stdio.h>

template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace Eigen;
  
  DATA_VECTOR(Y);
  DATA_IVECTOR(locs_index);
  DATA_MATRIX(K_inv);
  
  PARAMETER_VECTOR(vals);
  
  Type f=0;
  int n = Y.size();
  Type sd = 0.1;
  for(int j=0; j<n; j++){
    f -= dnorm(Y(j), vals(locs_index(j)), sd, true);
  }
  f += (vals * (K_inv * vals)).sum();

  return(f);
}
