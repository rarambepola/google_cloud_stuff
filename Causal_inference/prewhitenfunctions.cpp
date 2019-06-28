//file with all the functions needed to prewhiten time-space data
//these are NOT the regression functions used for HSIC/FSIC 
//only calculating residuals, no prediction to keep it simple

#include <RcppEigen.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
using namespace Rcpp;
using namespace Eigen;

const double pi = 3.141592653589793238462643383279502884;

// [[Rcpp::depends(RcppEigen)]]

// //2D radial basis kernel (for distances in 2d space)
// //input vectors (x1, y1) and (x2, y2)
// double rbf2d(double x1, double y1, double x2, double y2, double sigma=1){
//   return(exp(-(pow(x1 - x2, 2) + pow(y1 - y2, 2)) / (2*sigma)));
// }
// 
// //1D radial basis kernel (for non-periodic time distance)
// double rbf1d(double t1, double t2, double sigma=1){
//   return(exp(-abs(t1-t2) / (2*sigma)));
// }
// 
// //powered exponential kernel on S1 (circle)
// //see https://arxiv.org/pdf/1111.7077.pdf
// //alpha should be in (0, 1]
// double powexp(double t1, double t2, double period, double c=1, double alpha=1){
//   //rescale to make period = 2PI
//   double dist = abs(t1 - t2) * 2 * pi / period;
//   //get great circle distance
//   if(dist > pi){
//     dist = 2*pi - dist;
//   }
//   return(exp(-pow(dist/c, alpha)));
// }
// 
// 
// //"spiral" kernel - combination of 1d rbd and circle powexp
// double spiral(double t1, double t2, double period, double c=1, double alpha=1, double sigma=1){
//   return(powexp(t1,t2,period,c,alpha) * rbf1d(t1,t2,sigma));
// }

//2d eucl distance sqaured
double dist2dsq(double x1, double y1, double x2, double y2){
  return(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

//1d eucl distance
double dist1d(double t1, double t2){
  return(abs(t1 - t2));
}

//great circle distance
double distgc(double t1, double t2, double period){
  //rescale to make period = 2PI
  double dist = abs(t1 - t2) * 2 * pi / period;
  //get great circle distance
  while(dist > 2*pi){
    dist = dist - 2*pi;
  }
  if(dist > pi){
    dist = 2*pi - dist;
  }
  return(dist);
}


//2D radial basis kernel (for distances in 2d space)
//input vectors (x1, y1) and (x2, y2)
double rbf2d(double dist, double sigma=1){
  return(exp(- dist / (2*sigma)));
}

//1D radial basis kernel (for non-periodic time distance)
double rbf1d(double dist, double sigma=1){
  return(exp(- dist / (2*sigma)));
}

//powered exponential kernel on S1 (circle)
//see https://arxiv.org/pdf/1111.7077.pdf
//alpha should be in (0, 1]
double powexp(double dist, double c=1, double alpha=0.5){
  return(exp(-pow(dist/c, alpha)));
}


//function to perform regression
// [[Rcpp::export]]
NumericVector prewhiten(NumericVector vals, NumericVector x, NumericVector y, NumericVector t,
                        double period=365, double alpha=0.5){
  int n = vals.size();
  int n_dist = n * (n - 1) / 2;
  
  
  //create distance matrix and distance vector
  MatrixXd dist_mat_eucl2(n, n);
  MatrixXd dist_mat_eucl1(n, n);
  MatrixXd dist_mat_gc(n, n);
  std::vector<double> dist_vec_eucl2(n_dist);
  std::vector<double> dist_vec_eucl1(n_dist);
  std::vector<double> dist_vec_gc(n_dist);
  double dist_eucl2, dist_eucl1, dist_gc;
  
  int count = 0;
  for(int i=0; i<n; i++){
    for(int j=0; j<i; j++){
      dist_eucl2 = dist2dsq(x[i], y[i], x[j], y[j]);
      dist_eucl1 = dist1d(t[i], t[j]);
      dist_gc = distgc(t[i], t[j], period=period);
      
      dist_mat_eucl2(i, j) = dist_eucl2;
      dist_mat_eucl1(i, j) = dist_eucl1;
      dist_mat_gc(i, j) = dist_gc;
      
      //std::cout<<distance<<"\n";
      dist_vec_eucl2[count] = dist_eucl2;
      dist_vec_eucl1[count] = dist_eucl1;
      dist_vec_gc[count] = dist_gc;
      count += 1;
    }
  }
  
  std::sort(dist_vec_eucl2.begin(), dist_vec_eucl2.end());
  std::sort(dist_vec_eucl1.begin(), dist_vec_eucl1.end());
  std::sort(dist_vec_gc.begin(), dist_vec_gc.end());
  
  //get median
  double med_2, med_1, med_gc;
  if(n_dist % 2 == 0){
    med_2 = (dist_vec_eucl2[n_dist / 2] + dist_vec_eucl2[(n_dist/2) - 1]) / 2;
    med_1 = (dist_vec_eucl1[n_dist / 2] + dist_vec_eucl1[(n_dist/2) - 1]) / 2;
    med_gc = (dist_vec_gc[n_dist / 2] + dist_vec_gc[(n_dist/2) - 1]) / 2;
  }else{
    med_2 = dist_vec_eucl2[(n_dist - 1) / 2];
    med_1 = dist_vec_eucl1[(n_dist - 1) / 2];
    med_gc = dist_vec_gc[(n_dist - 1) / 2];
  }
  double scale_2 = med_2;
  double scale_1 = med_1;
  double scale_gc = med_gc;
  
  if(scale_2==0){
    scale_2 = 1;
  }
  if(scale_1==0){
    scale_1 = 1;
  }
  if(scale_gc==0){
    scale_gc = 1;
  }
  
  //create matrix K
  MatrixXd K = MatrixXd(n, n);
  
  for(int i=0; i<n; i++){
    for(int j=0; j<i; j++){
      K(i, j) = rbf2d(dist_mat_eucl2(i, j), scale_2) *
      //  rbf1d(dist_mat_eucl1(i, j), scale_1) * powexp(dist_mat_gc(i, j), scale_gc, alpha);
          rbf1d(dist_mat_eucl1(i, j), scale_1) * powexp(dist_mat_gc(i, j), 100.0, alpha);
      K(j, i) = K(i, j);
    }
    K(i, i) = 1;
  }
  
  VectorXd vals_vec(n);
  for(int i=0; i<n; i++){
    vals_vec[i] = vals[i];
  }
  
  double lambda = 0.7;
  MatrixXd final_K = K + lambda * MatrixXd::Identity(n, n);
  
  VectorXd vals_pred = K * final_K.llt().solve(vals_vec); 
  
  
  NumericVector res(n);
  for(int i=0; i<n; i++){
    res[i] = vals_pred[i] - vals[i];
  }

  return(res);
}


