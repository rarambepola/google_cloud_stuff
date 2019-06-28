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
  
  DATA_VECTOR(Y);
  DATA_VECTOR(pops);
  DATA_MATRIX(K_inv);
  DATA_IVECTOR(locs_index)
  
  PARAMETER_VECTOR(vals);
  // PARAMETER(beta0);
  
  
  
  int n = Y.size();
  
  Type f=0;
  for(int i=0; i<n; i++){
    // f -= dpois(Y(i), exp(beta0+vals(i)), true);
    f -= dpois(Y(i), exp(vals(locs_index(i))) * pops(i), true);
  }
  
  // f -= dnorm(beta0, (Type) 0, (Type) 0.001, true);
  
  // f += (vals * vals).sum();
  f += (vals * (K_inv * vals)).sum();
  
  // DATA_STRUCT(X, LOM_t);
  // DATA_MATRIX(Y);
  // DATA_VECTOR(pops);
  // DATA_SPARSE_MATRIX(A);
  // DATA_STRUCT(spde,spde_t);
  // 
  // PARAMETER(beta0);
  // PARAMETER_VECTOR(beta);
  // PARAMETER_MATRIX(S);
  // PARAMETER(log_kappa);
  // PARAMETER(log_tau);
  // PARAMETER(log_scale);
  // 
  // 
  // Type scale = exp(log_scale);
  // Type kappa = exp(log_kappa);
  // SparseMatrix<Type> Q = Q_spde(spde, kappa);
  // 
  // Type f=0;
  // int n_s = Y.rows();
  // int n_t = Y.cols();
  // int n_field = S.cols();
  // 
  // 
  // Type log_kappa_mean = 0.0;
  // Type log_kappa_sd = 0.01;
  // Type log_tau_mean = 1.0;
  // Type log_tau_sd = 0.01;
  // 
  // 
  // f -= dnorm(log_kappa, log_kappa_mean, log_kappa_sd, true);
  // f -= dnorm(log_tau, log_tau_mean, log_tau_sd, true);
  // 
  // 
  // 
  // //priors
  // Type beta_mean = 0.0;
  // Type beta_sd = 1.0;
  // int n_beta = beta.size();
  // 
  // f -= dnorm(beta0, beta_mean, beta_sd, true);
  // for(int i=0; i<n_beta; i++){
  //   f -= dnorm(beta[i], beta_mean, beta_sd, true);
  // }
  // 
  // // int n_months = 12;
  // matrix<Type> field = (A * S.transpose());
  // 
  // matrix<Type> Xb;
  // Type case_pred;
  // 
  // 
  // for(int i=0; i<n_t; i++){
  //   Xb = X(i) * beta;
  //   for(int j=0; j<n_s; j++){
  //     case_pred = exp(Xb(j) + field(j, i) + beta0) * pops(j);
  //     f -= dpois(Y(j, i), case_pred, true);
  //   }
  // }
  // 
  // //need array to pass to function
  // array<Type> S_2(n_t, n_field);
  // 
  // for(int i=0; i<n_t; i++){
  //   for(int j=0; j<n_field; j++){
  //     S_2(i, j) = S(i, j);
  //   }
  // }
  // matrix<Type> Sigma_t(n_t, n_t);
  // 
  // 
  // Type scale_mean = - 1.0;
  // Type scale_sd = 0.1;
  // 
  // 
  // 
  // f -= dnorm(log_scale, scale_mean, scale_sd, true);
  // 
  // 
  // for(int i=0; i<n_t; i++){
  //   for(int j=0; j<n_t; j++){
  //     Sigma_t(i, j) = cov((Type) i, (Type) j, scale);
  //   }
  // }
  // 
  // 
  // f += SEPARABLE(SCALE(GMRF(Q), 1 / exp(log_tau)), MVNORM(Sigma_t))(S_2);
  // 
  // // REPORT(Q);
  // // REPORT(log_tau);
  // // REPORT(Sigma_beta);
  return(f);
  
  }
