#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <cmath>
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
  const double alpha1 = params[0];
  const double beta1  = params[1];
  const double gamma1 = params[2];

  const int n = data.nrow();
  for (int i = 1; i < n; ++i) {
    if (data(i, 0) < data(i - 1, 0)) stop("data must be sorted by time (nondecreasing).");
  }

  double temporal_result = 0.0;

  // first_idx tracks the first index of the current time value
  int first_idx = 0;
  for (int i = 1; i < n; ++i) {
    if (data(i, 0) > data(i - 1, 0)) first_idx = i; // new (strictly larger) time value
    temporal_result += alpha1 + beta1 * data(i, 0) - gamma1 * static_cast<double>(first_idx);
  }

  return temporal_result;
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
  // unpack params
  const double alpha1 = params[0];
  const double beta1  = params[1];
  const double gamma1 = params[2];
  const double alpha2 = params[3];
  const double beta2  = params[4];
  const double alpha3 = params[5];
  const double beta3  = params[6];
  const double gamma3 = params[7];

  // basic guards (avoid NaNs from log(alpha2) etc.)
  if (alpha2 <= 0.0) return std::numeric_limits<double>::infinity();
  if (beta3  <  0.0) return std::numeric_limits<double>::infinity();
  if (gamma3 <  0.0) return std::numeric_limits<double>::infinity();

  const int n = data.nrow();
  const int K = xgrid.size();
  const int L = ygrid.size();
  const int M = tgrid.size();
  const int G = K * L;

  // assume data and tgrid are sorted (your current code implicitly assumes this via "break" logic)
  for (int i = 1; i < n; ++i) {
    if (data(i, 0) < data(i - 1, 0)) stop("data must be sorted by time (nondecreasing).");
  }
  for (int m = 1; m < M; ++m) {
    if (tgrid[m] < tgrid[m - 1]) stop("tgrid must be sorted (nondecreasing).");
  }

  // flatten (x,y) grid to g = 0..G-1
  std::vector<double> gx(G), gy(G);
  {
    int g = 0;
    for (int k = 0; k < K; ++k) {
      for (int l = 0; l < L; ++l) {
        gx[g] = xgrid[k];
        gy[g] = ygrid[l];
        ++g;
      }
    }
  }

  // pull data columns into std::vectors (faster than repeated data(i, j) inside tight loops)
  std::vector<double> tobs(n), xobs(n), yobs(n);
  for (int i = 0; i < n; ++i) {
    tobs[i] = data(i, 0);
    xobs[i] = data(i, 1);
    yobs[i] = data(i, 2);
  }

  // precompute distances from each event j to each grid point g:
  // store by event (contiguous in memory for the incremental event updates)
  const size_t GN = static_cast<size_t>(G) * static_cast<size_t>(n);
  std::vector<double> d2(GN);
  std::vector<double> log_r(GN); // log(sqrt(d2)) = 0.5*log(d2); -Inf if d2==0

  for (int j = 0; j < n; ++j) {
    const double xj = xobs[j];
    const double yj = yobs[j];
    const size_t base = static_cast<size_t>(j) * static_cast<size_t>(G);
    for (int g = 0; g < G; ++g) {
      const double dx = gx[g] - xj;
      const double dy = gy[g] - yj;
      const double r2 = dx*dx + dy*dy;
      d2[base + g] = r2;
      if (r2 > 0.0) {
        log_r[base + g] = 0.5 * std::log(r2);
      } else {
        log_r[base + g] = -std::numeric_limits<double>::infinity();
      }
    }
  }

  // constants for comparisons / log-phi updates
  const double a2_2   = alpha2 * alpha2;
  const double log_a2 = std::log(alpha2);
  const double b3_2   = beta3 * beta3;

  // state over grid points, updated incrementally as t increases
  std::vector<double> log_prod(G, 0.0); // log(prod phi)
  std::vector<int>    count(G, 0);      // theta3 neighborhood counts

  int hist_idx = 0; // number of events with t_j < t
  int lag_idx  = 0; // number of events with t_j <= t - gamma3

  double sum_over_t = 0.0;

  for (int m = 0; m < M; ++m) {
    const double t = tgrid[m];

    // bring history up to t: add events with tobs < t into log_prod
    while (hist_idx < n && tobs[hist_idx] < t) {
      const size_t base = static_cast<size_t>(hist_idx) * static_cast<size_t>(G);
      for (int g = 0; g < G; ++g) {
        const double r2 = d2[base + g];
        if (r2 <= a2_2) {
          // log(phi) = beta2 * (log(r) - log(alpha2)), with log(r)=0.5*log(r2)
          log_prod[g] += beta2 * (log_r[base + g] - log_a2);
        }
        // else phi=1 => log(phi)=0
      }
      ++hist_idx;
    }

    // update lagged-event set for theta3 counts: add events with tobs <= t - gamma3
    const double cutoff = t - gamma3;
    while (lag_idx < n && tobs[lag_idx] <= cutoff) {
      const size_t base = static_cast<size_t>(lag_idx) * static_cast<size_t>(G);
      for (int g = 0; g < G; ++g) {
        if (d2[base + g] <= b3_2) count[g] += 1;
      }
      ++lag_idx;
    }

    // ratio = (Σ_g exp(log_prod[g]) * exp(-alpha3*count[g])) / (Σ_g exp(log_prod[g]))
    double max_den = -std::numeric_limits<double>::infinity();
    double max_num = -std::numeric_limits<double>::infinity();

    for (int g = 0; g < G; ++g) {
      const double lp = log_prod[g];
      if (lp > max_den) max_den = lp;
      const double ln = lp - alpha3 * static_cast<double>(count[g]);
      if (ln > max_num) max_num = ln;
    }

    double den = 0.0;
    double num = 0.0;
    for (int g = 0; g < G; ++g) {
      den += std::exp(log_prod[g] - max_den);
      num += std::exp((log_prod[g] - alpha3 * static_cast<double>(count[g])) - max_num);
    }

    double ratio = 0.0;
    if (den > 0.0 && num > 0.0) {
      ratio = std::exp(max_num - max_den) * (num / den);
    }

    // temporal part uses N(t) = #{t_i < t} which is hist_idx by construction
    const double temporal = std::exp(alpha1 + beta1 * t - gamma1 * static_cast<double>(hist_idx));

    sum_over_t += temporal * ratio;
  }

  // IMPORTANT: this matches your *existing* estimator under the assumption of a uniform grid.
  // (bounds[1]*bounds[2]) cancels analytically, leaving only bounds[0]/M scaling.
  return (bounds[0] / static_cast<double>(M)) * sum_over_t;
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


//' calculates fast full self-correcting log-likelihood
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
double full_sc_lhood_fast(const NumericVector& xgrid,
                          const NumericVector& ygrid,
                          const NumericVector& tgrid,   // integration grid (M)
                          const NumericVector& tobs,    // observed times (n)
                          const NumericMatrix& data,    // (t,x,y), sorted by t
                          const NumericVector& params,  // (a1,b1,g1,a2,b2,a3,b3,g3)
                          const NumericVector& bounds)  // (bt,bx,by)
{
  // unpack params
  const double alpha1 = params[0];
  const double beta1  = params[1];
  const double gamma1 = params[2];
  const double alpha2 = params[3];
  const double beta2  = params[4];
  const double alpha3 = params[5];
  const double beta3  = params[6];
  const double gamma3 = params[7];

  const int n = data.nrow();
  const int K = xgrid.size();
  const int L = ygrid.size();
  const int M = tgrid.size();
  const int nObs = tobs.size();
  const int G = K * L;

  if (n < 2) return 0.0;
  if (G <= 0 || M <= 0 || nObs <= 0) stop("Empty grid or tobs/tgrid.");

  // require sorted by time (your current code assumes this via break-on-time logic)
  for (int i = 1; i < n; ++i) {
    if (data(i, 0) < data(i - 1, 0)) stop("data must be sorted by time (nondecreasing).");
  }
  for (int i = 1; i < M; ++i) {
    if (tgrid[i] < tgrid[i - 1]) stop("tgrid must be sorted (nondecreasing).");
  }
  for (int i = 1; i < nObs; ++i) {
    if (tobs[i] < tobs[i - 1]) stop("tobs must be sorted (nondecreasing).");
  }

  const double bt = bounds[0];
  const double bx = bounds[1];
  const double by = bounds[2];
  const double area = bx * by;

  // ------------------------------------------------------------
  // Fast part 1-1 (temporal log term): handle ties like strict "< t"
  // ------------------------------------------------------------
  double part_1_1 = 0.0;
  int first_idx = 0;
  for (int i = 1; i < n; ++i) {
    if (data(i, 0) > data(i - 1, 0)) first_idx = i; // new time block
    part_1_1 += alpha1 + beta1 * data(i, 0) - gamma1 * static_cast<double>(first_idx);
  }

  // ------------------------------------------------------------
  // Fast part 1-2: sum log(phi(||x_i-x_j||)) for j<i
  // Uses log form to avoid sqrt/pow/log in most cases.
  // ------------------------------------------------------------
  double part_1_2 = 0.0;
  if (alpha2 > 0.0 && beta2 != 0.0) {
    const double a2_2 = alpha2 * alpha2;
    const double log_a2 = std::log(alpha2);

    for (int i = 1; i < n; ++i) {
      const double xi = data(i, 1), yi = data(i, 2);
      for (int j = 0; j < i; ++j) {
        const double dx = xi - data(j, 1);
        const double dy = yi - data(j, 2);
        const double r2 = dx*dx + dy*dy;
        if (r2 <= a2_2) {
          if (r2 == 0.0) {
            // log(0^beta2) => -Inf for beta2>0 (rare unless duplicate points)
            return -std::numeric_limits<double>::infinity();
          }
          part_1_2 += beta2 * (0.5 * std::log(r2) - log_a2);
        }
      }
    }
  }
  // else alpha2<=0 or beta2==0 => phi ≡ 1 so part_1_2 = 0

  // ------------------------------------------------------------
  // Fast part 1-4: -alpha3 * sum_{i>j} 1[dist<=beta3 & lag>=gamma3]
  // ------------------------------------------------------------
  double part_1_4 = 0.0;
  if (alpha3 != 0.0 && beta3 >= 0.0 && gamma3 >= 0.0) {
    const double b3_2 = beta3 * beta3;
    double g_sum = 0.0;

    for (int i = 1; i < n; ++i) {
      const double ti = data(i, 0), xi = data(i, 1), yi = data(i, 2);
      for (int j = 0; j < i; ++j) {
        const double lag = ti - data(j, 0);
        if (lag >= gamma3) {
          const double dx = xi - data(j, 1);
          const double dy = yi - data(j, 2);
          const double r2 = dx*dx + dy*dy;
          if (r2 <= b3_2) g_sum += 1.0;
        }
      }
    }
    part_1_4 = -alpha3 * g_sum;
  }

  // ------------------------------------------------------------
  // Precompute grid flattening and event arrays
  // ------------------------------------------------------------
  std::vector<double> gx(G), gy(G);
  {
    int g = 0;
    for (int k = 0; k < K; ++k) {
      for (int l = 0; l < L; ++l) {
        gx[g] = xgrid[k];
        gy[g] = ygrid[l];
        ++g;
      }
    }
  }

  std::vector<double> te(n), xe(n), ye(n);
  for (int i = 0; i < n; ++i) {
    te[i] = data(i, 0);
    xe[i] = data(i, 1);
    ye[i] = data(i, 2);
  }

  // event-major storage for incremental updates
  const size_t GN = static_cast<size_t>(G) * static_cast<size_t>(n);
  std::vector<double> d2(GN);
  std::vector<double> log_r(GN);

  for (int j = 0; j < n; ++j) {
    const double xj = xe[j], yj = ye[j];
    const size_t base = static_cast<size_t>(j) * static_cast<size_t>(G);
    for (int g = 0; g < G; ++g) {
      const double dx = gx[g] - xj;
      const double dy = gy[g] - yj;
      const double r2 = dx*dx + dy*dy;
      d2[base + g] = r2;
      log_r[base + g] = (r2 > 0.0) ? (0.5 * std::log(r2))
        : (-std::numeric_limits<double>::infinity());
    }
  }

  // ------------------------------------------------------------
  // Unified sweep to compute:
  //   part_1_3 = sum_{i=2..n} log C_theta2(t_i)
  //   part_2   = (bt/M) * sum_{m=1..M} lambda(t_m) * E_h[ g(t_m, X) ]
  //
  // Original part_2 integrates:
  //   sum_{k,l} (prod/C) * g * (area/G) = (sum prod*g) / (sum prod)
  // so we never need C inside part_2. (This is the ratio-of-sums.)
  // Your current slow form does (prod/C) explicitly :contentReference[oaicite:4]{index=4}.
  // ------------------------------------------------------------

  struct Eval {
    double t;
    int type;   // 0 = tobs logC (part_1_3), 1 = tgrid integral (part_2)
    int idx;
  };

  std::vector<Eval> evals;
  evals.reserve((nObs > 1 ? (nObs - 1) : 0) + M);

  // tobs: skip the first like your existing part_1_3_full :contentReference[oaicite:5]{index=5}
  for (int i = 1; i < nObs; ++i) evals.push_back({tobs[i], 0, i});

  // tgrid: include all integration points
  for (int m = 0; m < M; ++m) evals.push_back({tgrid[m], 1, m});

  std::sort(evals.begin(), evals.end(),
            [](const Eval& a, const Eval& b) {
              if (a.t < b.t) return true;
              if (a.t > b.t) return false;
              return a.type < b.type;
            });

  std::vector<double> log_prod(G, 0.0);
  std::vector<int>    count(G, 0);

  int hist_idx = 0; // events with te < t
  int lag_idx  = 0; // events with te <= t - gamma3

  const double b3_2 = beta3 * beta3;

  // phi constants
  const bool use_phi = (alpha2 > 0.0 && beta2 != 0.0);
  const double a2_2 = use_phi ? (alpha2 * alpha2) : 0.0;
  const double log_a2 = use_phi ? std::log(alpha2) : 0.0;

  double part_1_3 = 0.0;
  double integral_sum = 0.0;

  size_t p = 0;
  while (p < evals.size()) {
    const double t = evals[p].t;

    // update history for prod: include events with te < t
    while (hist_idx < n && te[hist_idx] < t) {
      if (use_phi) {
        const size_t base = static_cast<size_t>(hist_idx) * static_cast<size_t>(G);
        for (int g = 0; g < G; ++g) {
          const double r2 = d2[base + g];
          if (r2 <= a2_2) {
            log_prod[g] += beta2 * (log_r[base + g] - log_a2);
          }
        }
      }
      ++hist_idx;
    }

    // update lagged set for theta3 counts: include events with te <= t - gamma3
    const double cutoff = t - gamma3;
    while (lag_idx < n && te[lag_idx] <= cutoff) {
      const size_t base = static_cast<size_t>(lag_idx) * static_cast<size_t>(G);
      for (int g = 0; g < G; ++g) {
        if (d2[base + g] <= b3_2) count[g] += 1;
      }
      ++lag_idx;
    }

    // Compute logC(t) and ratio(t) once for this time t
    double max_den = -std::numeric_limits<double>::infinity();
    double max_num = -std::numeric_limits<double>::infinity();

    for (int g = 0; g < G; ++g) {
      const double lp = log_prod[g];
      if (lp > max_den) max_den = lp;
      const double ln = lp - alpha3 * static_cast<double>(count[g]);
      if (ln > max_num) max_num = ln;
    }

    double den = 0.0;
    double num = 0.0;
    for (int g = 0; g < G; ++g) {
      den += std::exp(log_prod[g] - max_den);
      num += std::exp((log_prod[g] - alpha3 * static_cast<double>(count[g])) - max_num);
    }

    const double log_sum_prod = max_den + std::log(den);
    const double logC = log_sum_prod + std::log(area) - std::log(static_cast<double>(G));

    double ratio = 0.0;
    if (den > 0.0 && num > 0.0) {
      ratio = std::exp(max_num - max_den) * (num / den);
    }

    // apply all evals at this same t
    while (p < evals.size() && evals[p].t == t) {
      if (evals[p].type == 0) {
        part_1_3 += logC;
      } else {
        // N(t) = #{te < t} = hist_idx
        const double temporal = std::exp(alpha1 + beta1 * t - gamma1 * static_cast<double>(hist_idx));
        integral_sum += temporal * ratio;
      }
      ++p;
    }
  }

  const double part_1 = part_1_1 + part_1_2 - part_1_3 + part_1_4;
  const double part_2 = (bt / static_cast<double>(M)) * integral_sum;

  return part_1 - part_2;
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
