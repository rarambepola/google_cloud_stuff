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


template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  DATA_MATRIX(Y);
  DATA_VECTOR(pops);
  DATA_MATRIX(cov_mat);
  DATA_IMATRIX(locs_index);
  DATA_SCALAR(log_sigma);
  // DATA_SCALAR(lambda);
  DATA_MATRIX(w);
  DATA_SPARSE_MATRIX(A);
  DATA_STRUCT(spde, spde_t);
  
  PARAMETER_VECTOR(vals);
  PARAMETER(beta0);
  PARAMETER_VECTOR(S);
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  // PARAMETER(log_sigma);
  // PARAMETER(log_lambda);
  
  
  //spatial field
  vector<Type> field = A * S;
  
  Type log_kappa_mean = 0.0;
  Type log_kappa_sd = 0.01;
  Type log_tau_mean = 1.0;
  Type log_tau_sd = 0.01;
  
  // Type log_sigma=0.0;
  Type log_lambda=0.0;
  
  Type f=0;
  f -= dnorm(log_kappa, log_kappa_mean, log_kappa_sd, true);
  f -= dnorm(log_tau, log_tau_mean, log_tau_sd, true);
  
  // f -= dnorm(log_lambda, (Type) 0.0, (Type) 1.0);
  // f -= dnorm(log_sigma, (Type) 5.0, (Type) 0.1);
  
  Type sigma=exp(log_sigma);
  Type lambda=exp(log_lambda);
  // int n_covs = cov_mat.cols();
  
  int n = Y.rows();
  int n_t = Y.cols();
  int n_points = cov_mat.rows();
  int D = w.cols();
  
  w = w*sigma;
  
  //construct Z
  matrix<Type> Z(n_points, 2*D);
  matrix<Type> cov_mat_w = cov_mat * w;

  for(int i=0; i<n_points; i++){
    for(int j=0; j<D; j++){
      Z(i, j) = cos(cov_mat_w(i, j)) / sqrt(D);
      Z(i, j + D) = sin(cov_mat_w(i, j)) / sqrt(D);
    }
  }
  
  //construct K_inv
  matrix<Type> I_n(n_points, n_points);
  matrix<Type> I_s(2*D, 2*D);

  // I_n.Identity();

  for(int i=0; i<n_points; i++){
    for(int j=0; j<n_points; j++){
      I_n(i, j) = 0;
    }
    I_n(i, i) = 1;
  }

  for(int i=0; i<2*D; i++){
    for(int j=0; j<2*D; j++){
      I_s(i, j) = 0;
    }
    I_s(i, i) = lambda;
  }

  matrix<Type> K_inv = I_n - (Z * (I_s + Z.transpose() * Z).householderQr().solve(Z.transpose()));
  K_inv = K_inv / lambda;
  

  
  
  f -= dnorm(beta0, (Type) 0.0, (Type) 1.0, true);
  
  for(int i=0; i<n_t; i++){
    for(int j=0; j<n; j++){
      f -= dpois(Y(j, i), exp(field(j) + beta0 + vals(locs_index(j, i))) * pops(i), true);
    }
  }
  // for(int i=0; i<n; i++){
  //   // f -= dpois(Y(i), exp(beta0+vals(i)), true);
  //   f -= dpois(Y(i), exp(field(i) + beta0 + vals(locs_index(i))) * pops(i), true);
  // }
  
  // f -= dnorm(beta0, (Type) 0, (Type) 0.001, true);
  
  // f += (vals * vals).sum();
  
  
  f += (vals * (K_inv * vals)).sum();
  
  SparseMatrix<Type> Q = Q_spde(spde, exp(log_kappa));
  f += SCALE(GMRF(Q), 1 / exp(log_tau))(S);
  
  
  return(f);
  
  }
