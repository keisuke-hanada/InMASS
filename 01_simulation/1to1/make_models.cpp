#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


arma::rowvec colVars(const arma::mat& X) {
  arma::rowvec m = arma::mean(X,0);
  arma::rowvec v = arma::sum(arma::square(X.each_row() - m),0)/(X.n_rows-1);
  return v;
}


arma::vec generate_covariate(int n_k, std::string dist_type, double treatment_effect) {
  arma::vec out(n_k, arma::fill::none);
  
  if (dist_type == "normal") {

    for (int i=0; i<n_k; i++) {
      double z = R::rnorm(0.0,1.0);
      out[i] = z + treatment_effect;
    }
    
  } else if (dist_type == "chi2") {

    for (int i=0; i<n_k; i++) {
      double z1 = R::rnorm(0.0,1.0);
      double z2 = R::rnorm(0.0,1.0);
      double chi3 = (z1*z1 + z2*z2)/2.0;
      out[i] = chi3 + treatment_effect - 1.0;
    }
    
  } else {

    for (int i=0; i<n_k; i++) {
      double z = R::rnorm(0.0,1.0);
      out[i] = z + treatment_effect;
    }
  }
  
  return out;
}

// [[Rcpp::export]]
List make_models(int n, double sigma, int strata, int nsim, 
                 std::string dist_x2 = "normal", std::string dist_x3 = "normal") {
  // 固定パラメータ
  arma::vec beta = {1.0, 2.0, -1.0, 0.0, 0.5}; // (Intercept, x1k, x2k, x3k, x1k*x2k)
  

  List params = List::create(
    Named("n")=n,
    Named("formula")="yik ~ 1 + x1k + x2k + x1k:x2k",
    Named("beta")=beta,
    Named("sigma")=sigma,
    Named("strata")=strata,
    Named("nsim")=nsim,
    Named("dist_x2")=dist_x2,
    Named("dist_x3")=dist_x3
  );
  

  arma::ivec n_strata_vec(strata);
  {
    NumericVector rnd = runif(strata, (double)n/2.0, 2.0*(double)n);
    for (int i=0; i<strata; i++) {
      n_strata_vec[i] = 2 * (int)std::round(rnd[i]);
    }
  }
  

  arma::vec treatment_strata_x2(strata, arma::fill::none);
  arma::vec treatment_strata_x3(strata, arma::fill::none);
  for (int i=0; i<strata; i++) {

    treatment_strata_x2[i] = 4.0*((double)i)/((double)strata-1.0)-1.0;
    treatment_strata_x3[i] = 0.0;
  }
  

  int total_n_ipd = n * nsim;
  NumericVector IPD_x1k(total_n_ipd), IPD_x2k(total_n_ipd), IPD_x3k(total_n_ipd), IPD_yik(total_n_ipd), IPD_nsim(total_n_ipd);
  
  int sum_n_strata = 0;
  for (int i=0; i<strata; i++) sum_n_strata += n_strata_vec[i];
  int total_n_strata_ipd = sum_n_strata * nsim;
  int total_n_strata_ad = strata * 4 * nsim;
  
  NumericVector SIPD_x1k(total_n_strata_ipd), SIPD_x2k(total_n_strata_ipd), SIPD_x3k(total_n_strata_ipd),
  SIPD_yik(total_n_strata_ipd), SIPD_strata(total_n_strata_ipd), SIPD_nsim(total_n_strata_ipd);
  
  NumericVector SAD_x1k(total_n_strata_ad), SAD_x2k(total_n_strata_ad), SAD_x3k(total_n_strata_ad),
  SAD_yik(total_n_strata_ad);
  IntegerVector SAD_strata(total_n_strata_ad), SAD_n(total_n_strata_ad);
  CharacterVector SAD_var(total_n_strata_ad);
  NumericVector SAD_nsim(total_n_strata_ad);
  

  IntegerVector arm_allocation(n);
  for (int i=0; i<n/2; i++) {
    arm_allocation[i*2] = 0;
    arm_allocation[i*2+1] = 1;
    
  }
  
  int ipd_pos = 0;
  int sipd_pos = 0;
  int sad_pos = 0;
  
  for (int sim_i=1; sim_i<=nsim; sim_i++) {
    arma::vec global_x2 = generate_covariate(n, dist_x2, 0.0);
    arma::vec global_x3 = generate_covariate(n, dist_x3, 0.0);
    
    arma::mat dat01_mat(n,3,arma::fill::none);
    // x1k: arm
    for (int i=0; i<n; i++) dat01_mat(i,0) = (double)arm_allocation[i];
    dat01_mat.col(1) = global_x2;
    dat01_mat.col(2) = global_x3;
    

    arma::vec linpred = beta[0] + dat01_mat * beta.subvec(1,3);

    for (int i=0; i<n; i++) linpred[i] += arm_allocation[i] * global_x2[i] * beta[4] + R::rnorm(0.0, sigma);
    
    for (int i=0; i<n; i++) {
      IPD_x1k[ipd_pos] = dat01_mat(i,0);
      IPD_x2k[ipd_pos] = dat01_mat(i,1);
      IPD_x3k[ipd_pos] = dat01_mat(i,2);
      IPD_yik[ipd_pos] = linpred[i];
      IPD_nsim[ipd_pos] = sim_i;
      ipd_pos++;
    }
    

    for (int k=0; k<strata; k++) {
      int n_k = n_strata_vec[k];
      int half_n_k = n_k/2;
      

      arma::vec arm_k(n_k, arma::fill::none);
      for (int i=0; i<half_n_k; i++) {
        arm_k[i] = 0.0;
        arm_k[half_n_k+i] = 1.0;
      }
      

      arma::vec x2_k = generate_covariate(n_k, dist_x2, treatment_strata_x2[k]);
      arma::vec x3_k = generate_covariate(n_k, dist_x3, treatment_strata_x3[k]);
      
      arma::mat dat11_mat(n_k,3,arma::fill::none);
      dat11_mat.col(0) = arm_k;
      dat11_mat.col(1) = x2_k;
      dat11_mat.col(2) = x3_k;
      

      arma::vec linpred_str = beta[0] + dat11_mat * beta.subvec(1,3);
      for (int i=0; i<n_k; i++) linpred_str[i] += arm_k[i] * x2_k[i] * beta[4] + R::rnorm(0.0, sigma);
      
      for (int i=0; i<n_k; i++) {
        SIPD_x1k[sipd_pos] = dat11_mat(i,0);
        SIPD_x2k[sipd_pos] = dat11_mat(i,1);
        SIPD_x3k[sipd_pos] = dat11_mat(i,2);
        SIPD_yik[sipd_pos] = linpred_str[i];
        SIPD_strata[sipd_pos] = k+1;
        SIPD_nsim[sipd_pos] = sim_i;
        sipd_pos++;
      }
      
      arma::uvec idx0 = arma::find(dat11_mat.col(0)==0.0);
      arma::uvec idx1 = arma::find(dat11_mat.col(0)==1.0);
      
      arma::mat arm0_X = dat11_mat.rows(idx0);
      arma::mat arm1_X = dat11_mat.rows(idx1);
      
      arma::vec arm0_y = linpred_str.elem(idx0);
      arma::vec arm1_y = linpred_str.elem(idx1);
      
      arma::rowvec m0_x = arma::mean(arm0_X,0);
      arma::rowvec m1_x = arma::mean(arm1_X,0);
      double m0_y = arma::mean(arm0_y);
      double m1_y = arma::mean(arm1_y);
      
      arma::rowvec v0_x = colVars(arm0_X);
      arma::rowvec v1_x = colVars(arm1_X);
      double v0_y = arma::as_scalar(arma::sum(arma::square(arm0_y - m0_y)))/(half_n_k-1);
      double v1_y = arma::as_scalar(arma::sum(arma::square(arm1_y - m1_y)))/(half_n_k-1);
      

      SAD_x1k[sad_pos] = 0.0; SAD_x2k[sad_pos] = m0_x[1]; SAD_x3k[sad_pos] = m0_x[2]; SAD_yik[sad_pos] = m0_y;
      SAD_strata[sad_pos] = k+1; SAD_nsim[sad_pos] = sim_i; SAD_var[sad_pos] = "mean"; SAD_n[sad_pos] = half_n_k;
      sad_pos++;
      

      SAD_x1k[sad_pos] = 1.0; SAD_x2k[sad_pos] = m1_x[1]; SAD_x3k[sad_pos] = m1_x[2]; SAD_yik[sad_pos] = m1_y;
      SAD_strata[sad_pos] = k+1; SAD_nsim[sad_pos] = sim_i; SAD_var[sad_pos] = "mean"; SAD_n[sad_pos] = half_n_k;
      sad_pos++;
      
      // var arm0
      SAD_x1k[sad_pos] = 0.0; SAD_x2k[sad_pos] = v0_x[1]; SAD_x3k[sad_pos] = v0_x[2]; SAD_yik[sad_pos] = v0_y;
      SAD_strata[sad_pos] = k+1; SAD_nsim[sad_pos] = sim_i; SAD_var[sad_pos] = "var"; SAD_n[sad_pos] = half_n_k;
      sad_pos++;
      
      // var arm1
      SAD_x1k[sad_pos] = 1.0; SAD_x2k[sad_pos] = v1_x[1]; SAD_x3k[sad_pos] = v1_x[2]; SAD_yik[sad_pos] = v1_y;
      SAD_strata[sad_pos] = k+1; SAD_nsim[sad_pos] = sim_i; SAD_var[sad_pos] = "var"; SAD_n[sad_pos] = half_n_k;
      sad_pos++;
    }
  }
  
  DataFrame target_ipd = DataFrame::create(
    Named("x1k")=IPD_x1k,
    Named("x2k")=IPD_x2k,
    Named("x3k")=IPD_x3k,
    Named("yik")=IPD_yik,
    Named("nsim")=IPD_nsim
  );
  
  DataFrame strata_ipd = DataFrame::create(
    Named("x1k")=SIPD_x1k,
    Named("x2k")=SIPD_x2k,
    Named("x3k")=SIPD_x3k,
    Named("yik")=SIPD_yik,
    Named("strata")=SIPD_strata,
    Named("nsim")=SIPD_nsim
  );
  
  DataFrame strata_ad = DataFrame::create(
    Named("x1k")=SAD_x1k,
    Named("x2k")=SAD_x2k,
    Named("x3k")=SAD_x3k,
    Named("yik")=SAD_yik,
    Named("strata")=SAD_strata,
    Named("nsim")=SAD_nsim,
    Named("var")=SAD_var,
    Named("n")=SAD_n
  );
  
  return List::create(
    Named("target_ipd")=target_ipd,
    Named("strata_ipd")=strata_ipd,
    Named("strata_ad")=strata_ad,
    Named("params")=params
  );
}
