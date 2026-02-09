#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <cmath>
#include <unordered_map>
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
                          const NumericVector& tgrid,   // integration grid (M), sorted
                          const NumericVector& tobs,    // observed times (nObs), sorted
                          const NumericMatrix& data,    // (t,x,y), sorted by t and aligned with tobs
                          const NumericVector& params,  // (a1,b1,g1,a2,b2,a3,b3,g3)
                          const NumericVector& bounds)  // (bt,bx,by)
{
  const double alpha1 = params[0];
  const double beta1  = params[1];
  const double gamma1 = params[2];
  const double alpha2 = params[3];
  const double beta2  = params[4];
  const double alpha3 = params[5];
  const double beta3  = params[6];
  const double gamma3 = params[7];

  const int nx = xgrid.size();
  const int ny = ygrid.size();
  const int M  = tgrid.size();
  const int nObs = tobs.size();
  const int n = data.nrow();

  if (n < 2 || nObs < 2 || M < 1 || nx < 1 || ny < 1) {
    return -std::numeric_limits<double>::infinity();
  }
  if (n != nObs) stop("full_sc_lhood_fast: data.nrow() must equal length(tobs).");

  const double bt = bounds[0];
  const double bx = bounds[1];
  const double by = bounds[2];
  const double area = bx * by;
  const int G = nx * ny;

  // ---- Part 1_1 (tie-handling matches old) ----
  double part_1_1 = 0.0;
  int first_idx = 0;
  for (int i = 1; i < n; ++i) {
    if (data(i, 0) > data(i - 1, 0)) first_idx = i;
    part_1_1 += alpha1 + beta1 * data(i, 0) - gamma1 * static_cast<double>(first_idx);
  }

  // ---- Part 1_2 (matches old) ----
  double part_1_2 = 0.0;
  const bool use_phi = (alpha2 > 0.0 && beta2 != 0.0);
  const double a2_2 = use_phi ? (alpha2 * alpha2) : 0.0;
  const double log_a2 = use_phi ? std::log(alpha2) : 0.0;

  if (use_phi) {
    for (int i = 1; i < n; ++i) {
      const double xi = data(i, 1), yi = data(i, 2);
      for (int j = 0; j < i; ++j) {
        const double dx = xi - data(j, 1);
        const double dy = yi - data(j, 2);
        const double r2 = dx*dx + dy*dy;
        if (r2 <= a2_2) {
          if (r2 == 0.0) return -std::numeric_limits<double>::infinity();
          part_1_2 += beta2 * (0.5 * std::log(r2) - log_a2);
        }
      }
    }
  }

  // ---- Part 1_4 (matches old) ----
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

  // ---- Sweep for part_1_3 and part_2 (ratio-of-sums) ----
  std::vector<double> log_prod(static_cast<size_t>(G), 0.0);
  std::vector<int>    count(static_cast<size_t>(G), 0);

  int hist_idx = 0;
  int lag_idx  = 0;

  const double b3_2 = beta3 * beta3;

  double part_1_3 = 0.0;
  double integral_sum = 0.0;

  int ig = 0;
  int io = 1;

  // cache current (logC, ratio) and only recompute if state changed
  bool state_changed = true;
  double cached_logC = 0.0;
  double cached_ratio = 0.0;

  auto include_event_into_log_prod = [&](int e) {
    const double xe = data(e, 1);
    const double ye = data(e, 2);

    int g = 0;
    for (int ix = 0; ix < nx; ++ix) {
      const double dx = xgrid[ix] - xe;
      const double dx2 = dx * dx;
      for (int iy = 0; iy < ny; ++iy, ++g) {
        const double dy = ygrid[iy] - ye;
        const double r2 = dx2 + dy * dy;

        if (use_phi && r2 <= a2_2) {
          // IMPORTANT: match old behavior — hard fail if r2 == 0 inside cutoff
          if (r2 == 0.0) {
            return false; // signal -Inf likelihood
          }
          log_prod[static_cast<size_t>(g)] += beta2 * (0.5 * std::log(r2) - log_a2);
        }
      }
    }
    return true;
  };

  auto lag_event_into_count = [&](int e) {
    const double xe = data(e, 1);
    const double ye = data(e, 2);

    int g = 0;
    for (int ix = 0; ix < nx; ++ix) {
      const double dx = xgrid[ix] - xe;
      const double dx2 = dx * dx;
      for (int iy = 0; iy < ny; ++iy, ++g) {
        const double dy = ygrid[iy] - ye;
        const double r2 = dx2 + dy * dy;
        if (r2 <= b3_2) count[static_cast<size_t>(g)] += 1;
      }
    }
  };

  auto recompute_logC_and_ratio = [&]() {
    double max_den = -std::numeric_limits<double>::infinity();
    double max_num = -std::numeric_limits<double>::infinity();

    for (int g = 0; g < G; ++g) {
      const double lp = log_prod[static_cast<size_t>(g)];
      if (lp > max_den) max_den = lp;
      const double ln = lp - alpha3 * static_cast<double>(count[static_cast<size_t>(g)]);
      if (ln > max_num) max_num = ln;
    }

    double den = 0.0;
    double num = 0.0;
    for (int g = 0; g < G; ++g) {
      den += std::exp(log_prod[static_cast<size_t>(g)] - max_den);
      num += std::exp((log_prod[static_cast<size_t>(g)] - alpha3 * static_cast<double>(count[static_cast<size_t>(g)])) - max_num);
    }

    if (den <= 0.0) {
      // sum prod is 0 => logC=-Inf, ratio=0
      cached_logC = -std::numeric_limits<double>::infinity();
      cached_ratio = 0.0;
      return;
    }

    const double log_sum_prod = max_den + std::log(den);
    cached_logC = log_sum_prod + std::log(area) - std::log(static_cast<double>(G));

    if (num > 0.0) {
      cached_ratio = std::exp(max_num - max_den) * (num / den);
    } else {
      cached_ratio = 0.0;
    }
  };

  while (ig < M || io < nObs) {
    const double next_g = (ig < M)    ? tgrid[ig] : R_PosInf;
    const double next_o = (io < nObs) ? tobs[io]  : R_PosInf;
    const double t = (next_g < next_o) ? next_g : next_o;

    // include events with te < t
    while (hist_idx < n && data(hist_idx, 0) < t) {
      if (!include_event_into_log_prod(hist_idx)) {
        return -std::numeric_limits<double>::infinity();
      }
      ++hist_idx;
      state_changed = true;
    }

    // lagged events with te <= t - gamma3
    const double cutoff = t - gamma3;
    while (lag_idx < n && data(lag_idx, 0) <= cutoff) {
      lag_event_into_count(lag_idx);
      ++lag_idx;
      state_changed = true;
    }

    if (state_changed) {
      recompute_logC_and_ratio();
      state_changed = false;
    }

    // obs evals at time t
    while (io < nObs && tobs[io] == t) {
      part_1_3 += cached_logC;
      ++io;
    }

    // grid evals at time t
    while (ig < M && tgrid[ig] == t) {
      const double temporal = std::exp(alpha1 + beta1 * t - gamma1 * static_cast<double>(hist_idx));
      integral_sum += temporal * cached_ratio;
      ++ig;
    }
  }

  const double part_1 = part_1_1 + part_1_2 - part_1_3 + part_1_4;
  const double part_2 = (bt / static_cast<double>(M)) * integral_sum;

  return part_1 - part_2;
}


//' calculates spatial interaction
//'
//' @param hist a matrix of points (x,y), n x 2.
//' @param newp a new point vector (x,y), length 2.
//' @param params a vector of parameters (alpha2, beta2).
//'
//' @returns calculated probability of new point.
//' @keywords internal
// [[Rcpp::export]]
double spat_interaction(const Rcpp::NumericMatrix& hist,
                        const Rcpp::NumericVector& newp,
                        const Rcpp::NumericVector& params) {

  const double alpha2 = params[0];
  const double beta2  = params[1];

  if (alpha2 <= 0.0 || beta2 == 0.0) return 1.0;

  const int I = hist.nrow();
  const double a2_2 = alpha2 * alpha2;
  const double expo = 0.5 * beta2;

  const double x = newp[0];
  const double y = newp[1];

  double result = 1.0;

  for (int j = 0; j < I; ++j) {
    const double dx = x - hist(j, 0);
    const double dy = y - hist(j, 1);
    const double r2 = dx*dx + dy*dy;

    if (r2 <= a2_2) {
      result *= std::pow(r2 / a2_2, expo);
      if (result <= 0.0) break;
    }
  }

  return result;
}


//' Fast spatio-temporal interaction for the self-correcting model
//'
//' Computes $g_i = exp(-alpha3 * sum_\{j<i\} 1[ ||x_i-x_j|| <= beta3 AND (t_i - t_j) >= gamma3 ])$
//' for i = 1..n, with g_0 = exp(0) = 1.
//'
//'
//' @param data NumericMatrix with columns (time, x, y). Assumed sorted by time ascending.
//' @param params NumericVector length 3: (alpha3, beta3, gamma3)
//' @return NumericVector length n of exp(-alpha3 * counts)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector interaction_st_fast(const Rcpp::NumericMatrix& data,
                                        const Rcpp::NumericVector& params) {
  const int n = data.nrow();
  Rcpp::NumericVector out(n);
  if (n == 0) return out;

  const double alpha3 = params[0];
  const double beta3  = params[1];
  const double gamma3 = params[2];

  out[0] = 1.0;
  if (n == 1) return out;

  if (beta3 <= 0.0) {
    for (int i = 1; i < n; ++i) out[i] = 1.0;
    return out;
  }

  const double b2 = beta3 * beta3;

  // ---- small-n path: lag pointer + direct scan (usually faster up to a few hundred) ----
  if (n <= 300) {
    int lag_end = 0;
    for (int i = 1; i < n; ++i) {
      const double ti = data(i, 0);
      const double xi = data(i, 1);
      const double yi = data(i, 2);

      const double cutoff = ti - gamma3;
      while (lag_end < i && data(lag_end, 0) <= cutoff) ++lag_end;

      double count = 0.0;
      for (int j = 0; j < lag_end; ++j) {
        const double dx = xi - data(j, 1);
        const double dy = yi - data(j, 2);
        const double r2 = dx*dx + dy*dy;
        if (r2 <= b2) count += 1.0;
      }
      out[i] = std::exp(-alpha3 * count);
    }
    return out;
  }

  // ---- large-n path: hashed grid ----
  const double cell = beta3;

  auto key_of = [](int cx, int cy) -> long long {
    return (static_cast<long long>(cx) << 32) ^ static_cast<unsigned long long>(cy);
  };

  std::unordered_map<long long, std::vector<int>> cells;
  cells.reserve(static_cast<size_t>(n) * 2);

  int lag_ptr = 0;

  for (int i = 1; i < n; ++i) {
    const double ti = data(i, 0);
    const double xi = data(i, 1);
    const double yi = data(i, 2);

    const double cutoff = ti - gamma3;
    while (lag_ptr < i && data(lag_ptr, 0) <= cutoff) {
      const double xj = data(lag_ptr, 1);
      const double yj = data(lag_ptr, 2);
      const int cx = static_cast<int>(std::floor(xj / cell));
      const int cy = static_cast<int>(std::floor(yj / cell));
      cells[key_of(cx, cy)].push_back(lag_ptr);
      ++lag_ptr;
    }

    if (lag_ptr == 0) {
      out[i] = 1.0;
      continue;
    }

    const int cxi = static_cast<int>(std::floor(xi / cell));
    const int cyi = static_cast<int>(std::floor(yi / cell));

    double count = 0.0;
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        const auto it = cells.find(key_of(cxi + dx, cyi + dy));
        if (it == cells.end()) continue;
        const std::vector<int>& idxs = it->second;
        for (int idx : idxs) {
          const double dxp = xi - data(idx, 1);
          const double dyp = yi - data(idx, 2);
          const double r2  = dxp*dxp + dyp*dyp;
          if (r2 <= b2) count += 1.0;
        }
      }
    }

    out[i] = std::exp(-alpha3 * count);
  }

  return out;
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
Rcpp::NumericVector sim_temporal_sc(double Tmin = 0,
                                    double Tmax = 1,
                                    const Rcpp::NumericVector& params = Rcpp::NumericVector::create(0, 0, 0)) {
  const double alpha1 = params[0];
  const double beta1  = params[1];
  const double gamma1 = params[2];

  // Upper bound for intensity (simple, safe bound)
  const double Lmax = std::exp(alpha1 + beta1 * Tmax);
  const double N_max = Lmax * (Tmax - Tmin);

  const int sample_size = R::rpois(N_max);
  if (sample_size <= 0) return Rcpp::NumericVector(0);

  Rcpp::NumericVector sample = Rcpp::runif(sample_size, Tmin, Tmax);
  std::sort(sample.begin(), sample.end());

  Rcpp::NumericVector U = Rcpp::runif(sample_size, 0, 1);
  Rcpp::NumericVector t_sim(sample_size);

  // store accepted times (monotone increasing)
  std::vector<double> hist;
  hist.reserve(static_cast<size_t>(sample_size) + 1);
  hist.push_back(Tmin); // sentinel so N(t) = hist.size()-1

  for (int i = 0; i < sample_size; ++i) {
    const double t = sample[i];

    // N(t) = number of accepted times strictly less than t
    const double N_t = static_cast<double>(hist.size() - 1);

    const double lambda_t = std::exp(alpha1 + beta1 * t - gamma1 * N_t);

    if (U[i] < (lambda_t / Lmax)) {
      t_sim[i] = t;
      hist.push_back(t);
    } else {
      t_sim[i] = NA_REAL;
    }
  }

  return t_sim;
}


//' Simulate the spatial component of the self-correcting model (faster)
//'
//' @param M_n a vector of (x,y)-coordinates for anchor/first point.
//' @param params a vector of parameters (alpha_2, beta_2).
//' @param nsim_t number of points to simulate.
//' @param xy_bounds vector: (ax, bx, ay, by).
//'
//' @return a matrix nsim_t x 2 of point locations (x,y).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix sim_spatial_sc(const Rcpp::NumericVector& M_n,
                                   const Rcpp::NumericVector& params,
                                   int nsim_t,
                                   const Rcpp::NumericVector& xy_bounds) {

  // Always return a matrix with 2 cols (even if 0 rows) to avoid downstream subscript issues.
  Rcpp::NumericMatrix Loc(std::max(nsim_t, 0), 2);
  if (nsim_t <= 0) return Loc;

  const double ax = xy_bounds[0], bx = xy_bounds[1];
  const double ay = xy_bounds[2], by = xy_bounds[3];

  const double alpha2 = params[0];
  const double beta2  = params[1];

  // First point fixed at M_n
  Loc(0, 0) = M_n[0];
  Loc(0, 1) = M_n[1];

  int filled = 1;

  // If alpha2 <= 0 or beta2 == 0, interaction is identically 1 under your formula,
  // so everything is accepted: just draw iid uniform points.
  if (alpha2 <= 0.0 || beta2 == 0.0) {
    while (filled < nsim_t) {
      Loc(filled, 0) = R::runif(ax, bx);
      Loc(filled, 1) = R::runif(ay, by);
      ++filled;
    }
    return Loc;
  }

  const double a2_2 = alpha2 * alpha2;
  const double expo = 0.5 * beta2;

  while (filled < nsim_t) {
    const double x = R::runif(ax, bx);
    const double y = R::runif(ay, by);

    double prob = 1.0;

    // product over existing points
    for (int j = 0; j < filled; ++j) {
      const double dx = x - Loc(j, 0);
      const double dy = y - Loc(j, 1);
      const double r2 = dx*dx + dy*dy;

      if (r2 <= a2_2) {
        prob *= std::pow(r2 / a2_2, expo);
        if (prob == 0.0) break;
      }
    }

    if (R::runif(0.0, 1.0) <= prob) {
      Loc(filled, 0) = x;
      Loc(filled, 1) = y;
      ++filled;
    }
  }

  return Loc;
}
