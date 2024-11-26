#include <RcppArmadillo.h>
using namespace Rcpp;

//' calculates euclidean distance
//'
//' @param x a vector of x values.
//' @param y a vector of y values.
//'
//' @returns the distance between the two vectors.
//' @keywords internal
// [[Rcpp::export]]
double vec_dist(const NumericVector& x, const NumericVector& y)
{
  double out = sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2));
  return(out);
}


//' calculates full product for one grid point
//'
//' @param xgrid a vector of grid values for x.
//' @param ygrid a vector of grid values for y.
//' @param tgrid a t value.
//' @param data a matrix of data.
//' @param params a vector of parameters.
//'
//' @returns returns the product.
//' @keywords internal
// [[Rcpp::export]]
double full_product(const double xgrid, const double ygrid, const double tgrid,
                    const NumericMatrix& data, const NumericVector& params)
{
  double alpha2 = params[0];
  double beta2 = params[1];
  int n = data.nrow();
  double product_result = 1;
  double r;
  double phi_interaction;
  NumericVector xygrid = NumericVector::create(xgrid, ygrid);
  for(int j = 0; j < n; ++j)
  {
    if (data(j, 0) >= tgrid) break;
    NumericVector xydata = NumericVector::create(data(j, 1), data(j, 2));
    r = vec_dist(xygrid, xydata);
    phi_interaction = (r <= alpha2) * pow((r / alpha2), beta2) + (r > alpha2);
    product_result = product_result * phi_interaction;
  }
  return(product_result);
}


//' calculates c_theta
//'
//' @param xgrid a vector of grid values for x.
//' @param ygrid a vector of grid values for y.
//' @param tgrid a t value.
//' @param data a matrix of data.
//' @param params a vector of parameters.
//' @param bounds a vector of bounds for time, x, and y.
//'
//' @returns returns the product.
//' @keywords internal
// [[Rcpp::export]]
double C_theta2_i(const NumericVector& xgrid,
                  const NumericVector& ygrid,
                  const double tgrid,
                  const NumericMatrix& data,
                  const NumericVector& params,
                  const NumericVector& bounds)
{
  int K = xgrid.size();
  int L = ygrid.size();
  double i_result = 0;
  for(int k = 0; k < K; ++k)
  {
    for(int l = 0; l < L; ++l)
    {
      i_result +=  full_product(xgrid[k], ygrid[l], tgrid, data, params);
    }
  }
  return((i_result * bounds[1] * bounds[2])/(double(K) * double(L)));
}

//' calculates sum of values < t
//'
//' @param obs_t a vector of observed t values.
//' @param eval_t a t value.
//' @param y a vector of values.
//'
//' @returns the conditional sum.
//' @keywords internal
// [[Rcpp::export]]
double conditional_sum(const NumericVector& obs_t,
                       const double eval_t,
                       const NumericVector& y)
{
  double cond_sum = 0;
  int n = obs_t.size();
  for(int i = 0; i < n; i++)
  {
    if (obs_t[i] >= eval_t) break;
    cond_sum = cond_sum + y[i];
  }
  return(cond_sum);
}


//' calculates sum of values < t
//'
//' @param obs_t a vector of observed t values.
//' @param eval_t a t value.
//' @param y a vector of values.
//'
//' @returns the conditional sum.
//' @keywords internal
// [[Rcpp::export]]
double conditional_sum_logical(const NumericVector& obs_t,
                               const double eval_t,
                               const LogicalVector& y)
{
  double cond_sum = 0;
  int n = obs_t.size();
  for(int i = 0; i < n; i++)
  {
    if (obs_t[i] >= eval_t) break;
    cond_sum = cond_sum + y[i];
  }
  return(cond_sum);
}



//' calculates euclidean distance between a vector and a matrix
//'
//' @param eval_u a vector of x and y coordinates.
//' @param x_col a vector of x coordinates.
//' @param y_col a vector of y coordinates.
//'
//' @returns a vector of distances between a vector and each row of a matrix.
//' @keywords internal
// [[Rcpp::export]]
NumericVector vec_to_mat_dist(const NumericVector& eval_u,
                              const NumericVector& x_col,
                              const NumericVector& y_col)
{
  int n = x_col.size();
  NumericVector out(n);
  for(int i = 0; i < n; i++)
  {
    out[i] = sqrt(pow(eval_u[0] - x_col[i] , 2) + pow(eval_u[1] - y_col[i], 2));
  }
  return(out);
}

//' calculates distance in one dim
//'
//' @param eval_t a t value.
//' @param obs_t a vector of t values.
//'
//' @returns distance between a single t and the vector of all t values.
//' @keywords internal
// [[Rcpp::export]]
NumericVector dist_one_dim(const double eval_t, const NumericVector& obs_t)
{
  int n = obs_t.size();
  NumericVector out(n);
  for(int i = 0; i < n; i++)
  {
    out[i] = eval_t - obs_t[i];
  }
  return(out);
}

//' calculates part 1-1 full
//'
//' @param data a matrix of locations and times.
//' @param params a vector of parameters.
//'
//' @returns full likelihood evaluation for part 1.
//' @keywords internal
// [[Rcpp::export]]
double part_1_1_full(const NumericMatrix& data,
                     const NumericVector& params)
{
  double alpha1 = params[0];
  double beta1  = params[1];
  double gamma1 = params[2];
  int n = data.nrow();
  double temporal_result = 0;
  for (int i = 1; i < n; ++i)
  {
    temporal_result += alpha1 + beta1 * data(i, 0) - gamma1 * conditional_sum(data(_, 0), data(i, 0), rep(1.0, n)) ;
  }
  return(temporal_result);
}


//' calculates part 1-2 full
//'
//' @param data a matrix of locations and times.
//' @param params a vector of parameters.
//'
//' @returns full likelihood evaluation for part 2.
//' @keywords internal
// [[Rcpp::export]]
double part_1_2_full(const NumericMatrix& data,
                     const NumericVector& params)
{
  double alpha2 = params[0];
  double beta2 = params[1];
  int n = data.nrow();
  double part_1_2_result = 0;
  for(int i = 1; i < n; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      double r = sqrt(pow( (data(i, 1) - data(j, 1)), 2)
                    + pow( (data(i, 2) - data(j, 2)), 2));
      part_1_2_result += log(( r <= alpha2) * pow((r / alpha2), beta2) + (r > alpha2));
    }
  }
  return(part_1_2_result);
}


//' calculates part 1-3
//'
//' @param xgrid a vector of grid values for x.
//' @param ygrid a vector of grid values for y.
//' @param tgrid a t value.
//' @param data a matrix of times and locations.
//' @param params a vector of parameters.
//' @param bounds a vector of time, x, and y bounds.
//'
//' @returns full likelihood evaluation for part 3.
//' @keywords internal
// [[Rcpp::export]]
double part_1_3_full(const NumericVector& xgrid,
                     const NumericVector& ygrid,
                     const NumericVector& tgrid,
                     const NumericMatrix& data,
                     const NumericVector& params,
                     const NumericVector& bounds)
{
  double part_1_3_result = 0;
  int n = tgrid.size();
  for(int i = 1; i < n; ++i)
  {
    part_1_3_result += log(C_theta2_i(xgrid, ygrid, tgrid[i], data, params, bounds));
  }
  return(part_1_3_result);
}


//' calculates part 1-4
//'
//' @param data a matrix of times and locations.
//' @param params a vector of parameters.
//'
//' @returns full likelihood evaluation for part 4.
//' @keywords internal
// [[Rcpp::export]]
double part_1_4_full(const NumericMatrix& data,
                     const NumericVector& params)
{
  double alpha3 = params[0];
  double beta3  = params[1];
  double gamma3 = params[2];
  int n = data.nrow();
  double g_result = 0;
  double r, lag;
  for (int i = 1; i < n; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      r = sqrt(pow(data(i, 1) - data(j, 1) , 2) + pow(data(i, 2) - data(j, 2), 2));
      lag = data(i, 0) - data(j, 0);
      g_result += (r <= beta3) * (lag >= gamma3);
    }
  }
  return(-alpha3 * g_result);
}


//' calculates part 1 of the likelihood
//'
//' @param xgrid a vector of grid values for x.
//' @param ygrid a vector of grid values for y.
//' @param tgrid a t value.
//' @param data a matrix of times and locations.
//' @param params a vector of parameters.
//' @param bounds a vector of bounds for time, x, and y.
//'
//' @returns full evaluation of first part of likelihood.
//' @keywords internal
// [[Rcpp::export]]
double part_1_full(const NumericVector& xgrid,
                   const NumericVector& ygrid,
                   const NumericVector& tgrid,
                   const NumericMatrix& data,
                   const NumericVector& params,
                   const NumericVector& bounds)
{
  double alpha1 = params[0];
  double beta1  = params[1];
  double gamma1 = params[2];
  double alpha2 = params[3];
  double beta2  = params[4];
  double alpha3 = params[5];
  double beta3  = params[6];
  double gamma3 = params[7];
  double part_1_result;
  part_1_result = part_1_1_full(data, NumericVector::create(alpha1, beta1, gamma1))
                + part_1_2_full(data, NumericVector::create(alpha2, beta2))
                - part_1_3_full(xgrid, ygrid, tgrid, data,
                                NumericVector::create(alpha2, beta2),
                                bounds)
                + part_1_4_full(data, NumericVector::create(alpha3, beta3, gamma3));
  return(part_1_result);
}


//' calculates part 2 of the likelihood
//'
//' @param xgrid a vector of grid values for x.
//' @param ygrid a vector of grid values for y.
//' @param tgrid a vector of grid values for t.
//' @param data a matrix of times and locations.
//' @param params a vector of parameters.
//' @param bounds a vector of bounds for time, x, and y.
//'
//' @returns full evaluation of second part of likelihood.
//' @keywords internal
// [[Rcpp::export]]
double part_2_full(const NumericVector& xgrid,
                   const NumericVector& ygrid,
                   const NumericVector& tgrid,
                   const NumericMatrix& data,
                   const NumericVector& params,
                   const NumericVector& bounds)
{
  double alpha1 = params[0];
  double beta1  = params[1];
  double gamma1 = params[2];
  double alpha2 = params[3];
  double beta2  = params[4];
  double alpha3 = params[5];
  double beta3  = params[6];
  double gamma3 = params[7];

  int n = data.nrow();
  int K = xgrid.size();
  int L = ygrid.size();
  int M = tgrid.size();
  double i_result = 0;
  for(int m = 0; m < M; ++m)
  {
    for(int k = 0; k < K; ++k)
    {
      for(int l = 0; l < L; ++l)
      {
        // NumericVector twoDist = vec_to_mat_dist(NumericVector::create(xgrid[k], ygrid[l]), data(_, Range(1, 2)));
        NumericVector twoDist = vec_to_mat_dist(NumericVector::create(xgrid[k], ygrid[l]), data(_, 1), data(_, 2));
        NumericVector oneDist = dist_one_dim(tgrid[m], data(_, 0));
        i_result += (exp(alpha1 + beta1 * tgrid[m] - gamma1 * conditional_sum(data(_, 0), tgrid[m], rep(1.0, n)))
                      * ( full_product(xgrid[k], ygrid[l], tgrid[m], data, NumericVector::create(alpha2, beta2))
                        /C_theta2_i(xgrid,  ygrid,  tgrid[m], data, NumericVector::create(alpha2, beta2), bounds)
                        )
                      * exp(-alpha3 * conditional_sum_logical(data(_, 0), tgrid[m],( (twoDist <= beta3) * (oneDist >= gamma3) ))));
      }
     }
   }
  return(i_result * (bounds[0] * bounds[1] * bounds[2])/(double(M) * double(K) * double(L)));
}

//' calculates full self-correcting log-likelihood
//'
//' @param xgrid a vector of grid values for x.
//' @param ygrid a vector of grid values for y.
//' @param tgrid a vector of grid values for t.
//' @param tobs a vector of observed values for t.
//' @param data a matrix of times and locations.
//' @param params a vector of parameters.
//' @param bounds a vector of bounds for time, x, and y.
//'
//' @returns evaluation of full log-likelihood.
//' @keywords internal
// [[Rcpp::export]]
double full_sc_lhood(const NumericVector& xgrid,
                     const NumericVector& ygrid,
                     const NumericVector& tgrid,
                     const NumericVector& tobs,
                     const NumericMatrix& data,
                     const NumericVector& params,
                     const NumericVector& bounds)
{
  double full_likeli;

  full_likeli = part_1_full(xgrid, ygrid, tobs, data, params, bounds)
              - part_2_full(xgrid, ygrid, tgrid, data, params, bounds);

  return(full_likeli);
}



//' calculates spatial interaction
//'
//' @param Hist a matrix of points.
//' @param newp a new point vector.
//' @param params a vector of parameters.
//'
//' @returns calculated probability of new point.
//' @keywords internal
// [[Rcpp::export]]
double spat_interaction(const NumericMatrix& Hist,
                        const NumericVector& newp,
                        const NumericVector& params) {
  double alpha2 = params[0];
  double beta2  = params[1];
  double result = 1;
  int I = Hist.nrow();
  for(int j = 0; j < I; ++j){
    double r = sqrt(pow(newp[0] - Hist(j, 0), 2)
                  + pow(newp[1] - Hist(j, 1), 2));
    result *= ((r <= alpha2) * std::pow( (r / alpha2), beta2) + (r > alpha2));
  }
  return(result);
}


//' calculates spatio-temporal interaction
//'
//' @param data a matrix of times and locations.
//' @param params a vector of parameters.
//'
//' @returns a vector of interaction probabilities for every point.
//' @keywords internal
// [[Rcpp::export]]
NumericVector interaction_st(const NumericMatrix& data,
                             const NumericVector& params) {
  double alpha4 = params[0];
  double beta4  = params[1];
  double gamma4 = params[2];
  int n = data.nrow();
  NumericVector g_result(n);
  g_result[0] = 0;
  double sumval;
  double r;
  double lag;
  double indic;
  for (int i = 1; i < n; ++i)
  {
    sumval = 0;
    for(int j = 0; j < i; ++j)
    {
      r = sqrt(pow(data(i, 1) - data(j, 1) , 2) + pow(data(i, 2) - data(j, 2) , 2));
      lag = data(i, 0) - data(j, 0);
      indic = (r <= beta4) * (lag >= gamma4);
      sumval += indic;
    }
    g_result[i] = sumval;
  }
  return exp(-alpha4 * g_result);

}


//' calculates temporal likelihood
//'
//' @param params a vector of parameters (alpha_1, beta_1, gamma_1).
//' @param eval_t a t value.
//' @param obs_t a vector of t values.
//'
//' @returns evaluation of full temporal likelihood.
//' @keywords internal
// [[Rcpp::export]]
double temporal_sc(const NumericVector& params,
                   const double eval_t,
                   const NumericVector& obs_t) {
  double alpha1 = params[0];
  double beta1 =  params[1];
  double gamma1 = params[2];
  int N_row = obs_t.size();
  NumericVector oneDist = dist_one_dim(eval_t, obs_t);
  NumericVector ones = rep(1.0, N_row);
  double N_t = conditional_sum(obs_t, eval_t, ones);
  double result = exp(alpha1 + beta1 * eval_t - gamma1 * N_t);
  return(result);
}



//' Simulate the temporal component of the self-correcting model
//'
//' @param Tmin minimum time value.
//' @param Tmax maximum time value.
//' @param params a vector of parameters (alpha_1, beta_1, gamma_1).
//'
//' @return a vector of thinned and unthinned temporal samples.
//' @keywords internal
// [[Rcpp::export]]
NumericVector sim_temporal_sc(double Tmin = 0,
                              double Tmax = 1,
                              const NumericVector& params = NumericVector::create(0, 0, 0)) {
  double alpha1 = params[0];
  double beta1  = params[1];

  double Lmax = exp(alpha1 + beta1 * Tmax);
  double N_max = Lmax * Tmax;

  // Generate a sample size from a Poisson distribution
  int sample_size = R::rpois(N_max);

  // Generate uniformly distributed samples
  NumericVector sample = Rcpp::runif(sample_size, Tmin, Tmax);

  // Sort the sample in place
  std::sort(sample.begin(), sample.end());

  // Initialize history (hist) as an empty std::vector
  std::vector<double> hist;
  hist.push_back(0.0);  // Add an initial value

  // Generate a uniformly distributed vector U for acceptance-rejection
  NumericVector U = Rcpp::runif(sample_size, 0, 1);

  // Create output vector t_sim and lambda_star
  NumericVector t_sim(sample_size);
  NumericVector lambda_star(sample_size);

  // Main loop
  for (int i = 0; i < sample_size; ++i) {
    lambda_star[i] = temporal_sc(params, sample[i], NumericVector(wrap(hist)));

    // Acceptance-Rejection condition
    if (U[i] < (lambda_star[i] / Lmax)) {
      t_sim[i] = sample[i];
      hist.push_back(sample[i]);  // Update history
    } else {
      t_sim[i] = NA_REAL;  // Mark rejected samples as NA
    }
  }

  return t_sim;  // Return the final simulated times
}


//' Simulate the spatial component of the self-correcting model
//'
//' @param M_n a vector of (x,y)-coordinates for largest point.
//' @param params a vector of parameters (alpha_2, beta_2).
//' @param nsim_t number of points to simulate.
//' @param xy_bounds vector of lower and upper bounds for the domain (2 for x, 2 for y).
//'
//' @return a matrix of point locations in the (x,y)-plane.
//' @keywords internal
//[[Rcpp::export]]
NumericMatrix sim_spatial_sc(const NumericVector& M_n,
                             const NumericVector& params,
                             int nsim_t,
                             const NumericVector& xy_bounds){
  arma::mat Loc(1,2);
  Loc(0,0) = M_n(0);
  Loc(0,1) = M_n(1);
  int lengLoc = 1;
  while(lengLoc < nsim_t){
    arma::rowvec new_point(2);
    new_point(0) = R::runif(xy_bounds[0], xy_bounds[1]);
    new_point(1) = R::runif(xy_bounds[2], xy_bounds[3]);

    if( R::runif(0, 1) <= spat_interaction(Rcpp::NumericMatrix(Rcpp::wrap(Loc)), Rcpp::NumericVector(Rcpp::wrap(new_point)), params) ){
      Loc = arma::join_cols(Loc, new_point);
      lengLoc += 1;
    }
  }
  return(Rcpp::NumericMatrix(Rcpp::wrap(Loc)));
}
