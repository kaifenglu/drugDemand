#include <Rcpp.h>
using namespace Rcpp;

#ifndef __UTILITIES__
#define __UTILITIES__

double norm_rej(const double a, const double b);
double unif_rej(const double a, const double b);
double halfnorm_rej(const double a, const double b);
double exp_rej(const double a, const double b);
double rtnormcpp(const double mean,
                 const double sd,
                 const double lower,
                 const double upper);

#endif // __UTILITIES__
