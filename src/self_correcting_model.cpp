#include <Rcpp.h>
using namespace Rcpp;

//' calculates euclidean distance
//'
//' @param x a vector x
//' @param y a vector y
//'
//' @returns distance between the two vectors
//' @export
// [[Rcpp::export]]
double pdistCVec(NumericVector x, NumericVector y)
{
  double  out = sqrt(std::pow(x[0] - y[0] , 2) + std::pow(x[1] - y[1], 2));
  return out;
}


//' calculates full product
//'
//' @param xgrid a vector of grid values for x
//' @param ygrid a vector of grid values for y
//' @param tgrid a t value
//' @param data a matrix of data
//' @param params a vector of parameters
//'
//' @returns returns the product
//' @export
// [[Rcpp::export]]
double  prodFullCpp(double xgrid, double ygrid, double tgrid,
                    NumericMatrix data,  NumericVector params)
{
  double alpha2 = params[0];
  double beta2 = params[1];
  int n = data.nrow();
  double prodrslt=1;
  double r,interactPhi;
  NumericVector xygrid = NumericVector::create(xgrid, ygrid);
  for(int j = 0; j < n; ++j)
  {
    if (data(j, 0) >= tgrid) break;
    NumericVector xydata = NumericVector::create(data(j, 1), data(j, 2));
    r = pdistCVec(xygrid, xydata);
    interactPhi = (r <= alpha2) * pow((r / alpha2), beta2) + (r > alpha2);
    prodrslt = prodrslt * interactPhi;
  }
  return prodrslt;
}


//' calculates c_theta
//'
//' @param xgrid a vector of grid values for x
//' @param ygrid a vector of grid values for y
//' @param tgrid a t value
//' @param data a matrix of data
//' @param params a vector of parameters
//' @param bounds a vector of bounds for time, x, and y
//'
//' @returns returns the product
//' @export
// [[Rcpp::export]]
double Ctheta2i(NumericVector xgrid, NumericVector ygrid, double tgrid,
                NumericMatrix data, NumericVector params,
                NumericVector bounds)
{
  int K = xgrid.size();
  int L = ygrid.size();
  double Iresult = 0;
  for(int k = 0; k < K; ++k)
  {
    for(int l = 0; l < L; ++l)
    {
      Iresult +=  prodFullCpp(xgrid[k], ygrid[l], tgrid, data, params);
    }
  }
  return ((Iresult * bounds[1] * bounds[2])/(double(K) * double(L)));
}

//' calculates sum of values < t
//'
//' @param obst a vector of observed t
//' @param evalt a t value
//' @param y a vector of values
//'
//' @returns the sum
//' @export
// [[Rcpp::export]]
double CondSumCpp(NumericVector obst, double evalt, NumericVector y)
{
  double CondSum = 0;
  int n = obst.size();
  for(int i = 0; i < n; i++)
  {
    if (obst[i] >= evalt) break;
    CondSum = CondSum + y[i];
  }
  return CondSum;
}


//' calculates sum of values < t
//'
//' @param obst a vector of observed t
//' @param evalt a t value
//' @param y a vector of values
//'
//' @returns the sum
//' @export
// [[Rcpp::export]]
double CondSumCppR(NumericVector obst, double evalt, LogicalVector y)
{
  double CondSum=0;
  int n = obst.size();
  for(int i = 0; i < n; i++)
  {
    if (obst[i] >= evalt) break;
    CondSum = CondSum + y[i];
  }
  return CondSum;
}



//' calculates euclidean distance
//'
//' @param evalu a vector
//' @param obsu a matrix
//'
//' @returns distance between a vector and each row of a matrix
//' @export
// [[Rcpp::export]]
NumericVector pdistC(NumericVector evalu, NumericMatrix obsu)
{
  int n = obsu.nrow();
  NumericVector out(n);
  for(int i = 0; i < n; i++)
  {
    out[i] = sqrt(std::pow(evalu[0] - obsu(i, 0) , 2) + std::pow(evalu[1] - obsu(i, 1), 2));
  }
  return out;
}

//' calculates distance in one dim
//'
//' @param evalt a t value
//' @param obst a vector of t
//'
//' @returns distance between a t and all t
//' @export
// [[Rcpp::export]]
NumericVector rdistC(double evalt, NumericVector obst)
{
  int n = obst.size();
  NumericVector out(n);
  for(int i = 0; i < n; i++)
  {
    out[i] = evalt - obst[i];
  }
  return out;
}


//' calculates part 2
//'
//' @param xgrid a vector of grid values for x
//' @param ygrid a vector of grid values for y
//' @param tgrid a vector of grid values for t
//' @param data a matrix of times and locations
//' @param params a vector of parameters
//' @param bounds a vector of bounds for time, x, and y
//'
//' @returns second part of likelihood
//' @export
// [[Rcpp::export]]
double Part2FullCpp(NumericVector xgrid, NumericVector ygrid,
                    NumericVector tgrid, NumericMatrix data, NumericVector params,
                    NumericVector bounds)
{
  double alpha1 = params[0];
  double beta1  = params[1];
  double gamma1 = params[2];
  double alpha2 = params[3];
  double beta2  = params[4];
  double alpha3 = params[5];
  double beta3  = params[6];
  double gamma3 = params[7];
  Environment base( "package:base" ) ;
  Function repcpp = base["rep.int"];
  int n = data.nrow();
  int K = xgrid.size();
  int L = ygrid.size();
  int M = tgrid.size();
  double Iresult = 0;
  for(int m = 0; m < M; ++m)
  {
    for(int k = 0; k < K; ++k)
    {
      for(int l = 0; l < L; ++l)
      {
        NumericVector twoDist = pdistC(NumericVector::create(xgrid[k], ygrid[l]), data(_, Range(1, 2)));
        NumericVector oneDist = rdistC(tgrid[m], data(_, 0));
        Iresult +=( exp(alpha1 + beta1 * tgrid[m] - gamma1 * CondSumCpp(data(_, 0), tgrid[m], repcpp(1, n)))
                      * (
                          prodFullCpp(xgrid[k], ygrid[l], tgrid[m], data, NumericVector::create(alpha2, beta2))
                      /Ctheta2i(xgrid,  ygrid,  tgrid[m], data, NumericVector::create(alpha2, beta2), bounds)
                      )
                      * exp(-alpha3 * CondSumCppR(data(_, 0), tgrid[m],( (twoDist <= beta3) * (oneDist >= gamma3) ))));
      }
    }
  }
  return (Iresult * (bounds[0] * bounds[1] * bounds[2])/(double(M) * double(K) * double(L)));
}


//' calculates part 1-1 full
//'
//' @param data a matrix of locations and times
//' @param paramt a vector of parameters
//'
//' @returns full likelihood for part 1
//' @export
// [[Rcpp::export]]
double Part1_1FullCpp(NumericMatrix data, NumericVector paramt)
{
  double alpha1 = paramt[0];
  double beta1  = paramt[1];
  double gamma1 = paramt[2];
  Environment base( "package:base" );
  Function repcpp = base["rep.int"];
  int n=data.nrow();
  double temrslt = 0;
  for (int i = 1; i < n; ++i)
  {
    temrslt += alpha1 + beta1 * data(i, 0) - gamma1 * CondSumCpp(data(_, 0), data(i, 0), repcpp(1, n)) ;
  }
  return temrslt;
}


//' calculates part 1-2 full
//'
//' @param data a matrix of locations and times
//' @param params a vector of parameters
//'
//' @returns full likelihood for part 2
//' @export
// [[Rcpp::export]]
double Part1_2FullCpp(NumericMatrix data, NumericVector params)
{
  double alpha2 = params[0];
  double beta2 = params[1];
  int n = data.nrow();
  double p12result = 0;
  for(int i = 1; i < n; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      double r = sqrt(std::pow( (data(i, 1) - data(j, 1)), 2)
                        + std::pow( (data(i, 2) - data(j, 2)), 2));
      p12result += log(( r <= alpha2) * std::pow((r / alpha2), beta2) + (r > alpha2));
    }
  }
  return p12result;
}


//' calculates part 1-3
//'
//' @param xgrid a vector of grid values for x
//' @param ygrid a vector of grid values for y
//' @param tgrid a t value
//' @param data a matrix of times and locations
//' @param params a vector of parameters
//' @param bounds a vector of time, x, and y bounds
//'
//' @returns full likelihood for part 3
//' @export
// [[Rcpp::export]]
double Part1_3FullCpp(NumericVector xgrid, NumericVector ygrid,
                      NumericVector tgrid, NumericMatrix data, NumericVector params,
                      NumericVector bounds)
{
  double p13result = 0;
  int n = tgrid.size();
  for(int i = 1; i < n; ++i)
  {
    p13result += log(Ctheta2i(xgrid, ygrid, tgrid[i], data, params, bounds));
  }
  return p13result;
}


//' calculates part 1-4
//'
//' @param data a matrix of times and locations
//' @param paramg a vector of parameters
//'
//' @returns full likelihood for part 4
//' @export
// [[Rcpp::export]]
double Part1_4FullCpp(NumericMatrix data, NumericVector paramg)
{
  double alpha4 = paramg[0];
  double beta4  = paramg[1];
  double gamma4 = paramg[2];
  int n = data.nrow();
  double grslt = 0;
  double r, lag;
  for (int i = 1; i < n; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      r = sqrt(std::pow(data(i, 1) - data(j, 1) , 2) + std::pow(data(i, 2) - data(j, 2), 2));
      lag = data(i, 0) - data(j, 0);
      grslt += (r <= beta4) * (lag >= gamma4);
    }
  }
  return -alpha4 * grslt;
}


//' calculates part 1
//'
//' @param xgrid a vector of grid values for x
//' @param ygrid a vector of grid values for y
//' @param tgrid a t value
//' @param data a matrix of times and locations
//' @param params a vector of parameters
//' @param bounds a vector of bounds for time, x, and y
//'
//' @returns first part of likelihood
//' @export
// [[Rcpp::export]]
double Part1FullCpp(NumericVector xgrid, NumericVector ygrid, NumericVector tgrid,
                    NumericMatrix data, NumericVector params, NumericVector bounds)
{
  double alpha1 = params[0];
  double beta1  = params[1];
  double gamma1 = params[2];
  double alpha2 = params[3];
  double beta2  = params[4];
  double alpha3 = params[5];
  double beta3  = params[6];
  double gamma3 = params[7];
  double p1rslt;
  p1rslt = Part1_1FullCpp(data, NumericVector::create(alpha1, beta1, gamma1))
           + Part1_2FullCpp(data, NumericVector::create(alpha2, beta2))
           - Part1_3FullCpp(xgrid, ygrid, tgrid, data,
                            NumericVector::create(alpha2, beta2), bounds)
           + Part1_4FullCpp(data, NumericVector::create(alpha3, beta3, gamma3));
  return (p1rslt);
}


//' calculates full self-correcting log-likelihood
//'
//' @param xgrid a vector of grid values for x
//' @param ygrid a vector of grid values for y
//' @param tgrid a vector of grid values for t
//' @param tobs a vector of observed values for t
//' @param data a matrix of times and locations
//' @param params a vector of parameters
//' @param bounds a vector of bounds for time, x, and y
//'
//' @returns full log-likelihood
//' @export
// [[Rcpp::export]]
double full_sc_lhood(NumericVector xgrid, NumericVector ygrid, NumericVector tgrid,
                     NumericVector tobs, NumericMatrix data, NumericVector params,
                     NumericVector bounds)
{
  double full_likeli;

  full_likeli = Part1FullCpp(xgrid, ygrid, tobs, data, params, bounds)
                - Part2FullCpp(xgrid, ygrid, tgrid, data, params, bounds);

  return(full_likeli);
}



//' calculates spatial interaction
//'
//' @param Hist a matrix of points
//' @param newp a new point vector
//' @param par a vector of parameters
//'
//' @returns probability of new point
//' @export
// [[Rcpp::export]]
double interactionCpp(NumericMatrix Hist, NumericVector newp,
                      NumericVector par) {
  double alpha2 = par[0];
  double beta2  = par[1];
  double rslt = 1;
  int I = Hist.nrow();
  for(int j = 0; j < I; ++j){
    double r = sqrt(std::pow(newp[0] - Hist(j, 0), 2)
                      + std::pow(newp[1] - Hist(j, 1), 2));
    rslt *= ((r <= alpha2) * std::pow( (r / alpha2), beta2) + (r > alpha2));
  }
  return rslt;
}


//' calculates spatio-temporal interaction
//'
//' @param data a matrix of times and locations
//' @param paramg a vector of parameters
//'
//' @returns interaction probabilities for every point
//' @export
// [[Rcpp::export]]
NumericVector interactionCpp_st(NumericMatrix data, NumericVector paramg) {
  double alpha4 = paramg[0];
  double beta4  = paramg[1];
  double gamma4 = paramg[2];
  int n = data.nrow();
  NumericVector grslt(n);
  grslt[0] = 0;
  double sumval;
  double r, lag;
  double indic;
  for (int i = 1; i < n; ++i)
  {
    sumval = 0;
    for(int j = 0; j < i; ++j)
    {
      r = sqrt(std::pow(data(i, 1) - data(j, 1) , 2) + std::pow(data(i, 2) - data(j, 2) , 2));
      lag = data(i, 0) - data(j, 0);
      indic = (r <= beta4) * (lag >= gamma4);
      sumval += indic;
    }
    grslt[i] = sumval;
  }
  return exp(-alpha4 * grslt);

}
