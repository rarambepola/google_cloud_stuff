#include <TMB.hpp>
#include <math.h>
#include <stdio.h>

const double pi = 3.141592653589793238462643383279502884;

//define container for list of matrices - LOM type
template <class Type>
  struct LOM_t : vector<matrix<Type> > {
    //constructor
    LOM_t(SEXP x){ /*x is list from R*/
        (*this).resize(LENGTH(x)); /*resize pointer to be the right size*/
        for(int i=0; i<LENGTH(x); i++){ 
          SEXP m = VECTOR_ELT(x, i); /*get ith element of list */
            (*this)(i) = asMatrix<Type>(m); /*put it as ith element*/
        }
    }
  };

// template <class Type>
  // Type f_sc(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template <class Type>
  Type dist1d(Type t1, Type t2){
    return((t1 - t2) * (t1 - t2));
  }

//great circle distance
template <class Type>
  Type distgc(Type t1, Type t2, Type period){
    //rescale to make period = 2PI
    Type dist = abs(t1 - t2) * 2 * pi / period;
    //get great circle distance
    while(dist > 2*pi){
      dist = dist - 2*pi;
    }
    if(dist > pi){
      dist = 2*pi - dist;
    }
    dist = dist * period / (2 * pi);
    return(dist);
  }

template <class Type>
  Type cov(Type t1, Type t2, Type scale2=1, Type period=12){
    // return(exp(-(dist1d(t1, t2) / scale1) - (distgc(t1, t2, period) / scale2)));
    // return(exp(-distgc(t1, t2, period) / scale2));
    return(exp(-dist1d(t1, t2) / scale2));
  }


template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  DATA_STRUCT(X, LOM_t);
  DATA_MATRIX(Y);
  DATA_SPARSE_MATRIX(A);
  DATA_STRUCT(spde, spde_t);
  DATA_SCALAR(beta_scale);
  
  PARAMETER(beta0);
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(S);
  PARAMETER(log_kappa);
  PARAMETER(log_tau);

  
  Type f=0;
  int n_s = Y.rows();
  int n_t = Y.cols();
  // int n_field = S.size();
  
  vector<Type> field = A * S;
  
  Type log_kappa_mean = 0.0;
  Type log_kappa_sd = 0.01;
  Type log_tau_mean = 1.0;
  Type log_tau_sd = 0.01;
  
  
  f -= dnorm(log_kappa, log_kappa_mean, log_kappa_sd, true);
  f -= dnorm(log_tau, log_tau_mean, log_tau_sd, true);
  
  // int n_hf = hf_effect.size();
  
  

  //priors
  Type beta_mean = 0.0;
  Type beta_sd = 1.0;
  int n_beta = beta.size();
  
  f -= dnorm(beta0, beta_mean, beta_sd, true);
  for(int i=0; i<n_beta; i++){
    // f -= dnorm(beta[i], beta_mean, beta_sd, true);
    // f += abs(beta[i] - beta_mean) / beta_scale;
    if(beta[i] - beta_mean > 0){
      f += (beta[i] - beta_mean) / beta_scale;
    }else{
      f-= (beta[i] - beta_mean) / beta_scale;
    }
  }
  
  
  
  matrix<Type> Xb;
  Type prob_pred;
  
  
  for(int i=0; i<n_t; i++){
    Xb = X(i) * beta + field;
    for(int j=0; j<n_s; j++){
      prob_pred = Xb(j) + beta0;
      f -= dbinom_robust(Y(j, i), (Type) 1, prob_pred, true);
    }
  }
  
  SparseMatrix<Type> Q = Q_spde(spde, exp(log_kappa));
  f += SCALE(GMRF(Q), 1 / exp(log_tau))(S);
    
  return(f);
  
  }
