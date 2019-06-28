#include <TMB.hpp>
#include <math.h>
#include <stdio.h>



template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(X);
  DATA_VECTOR(Y);
  PARAMETER(beta);

  int n = X.size();
  
  Type f=0;
  
  for(int i=0; i<n; i++){
    f -= dnorm(Y(i), X(i) * beta, (Type) 1.0, true);
  }
  return(f);
  
}
