#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;

//' @title Dosing Date Imputation for Ongoing Patients
//' @description Imputes the dosing dates after cutoff for ongoing
//' patients with dosing records.
//'
//' @param usubjid The unique subject ID.
//' @param V The last dosing visit date relative to randomization.
//' @param C The cutoff date relative to randomization.
//' @param D The discontinuation date relative to randomization.
//' @param model_ki The model for the number of skipped
//'   visits between two consecutive drug dispensing visits.
//'   Options include "constant", "poisson", "zero-inflated poisson",
//'   and "negative binomial".
//' @param theta_ki The model parameters for the number of skipped visits
//'   between two consecutive drug dispensing visits.
//' @param model_ti The model for the gap time between two consecutive
//'   drug dispensing visits. Options include "least squares"
//'   and "least absolute deviations".
//' @param theta_ti The model parameters for the gap time between
//'   two consecutive drug dispensing visits.
//'
//' @return A data frame with two variables:
//'
//' * \code{usubjid}: The unique subject ID.
//'
//' * \code{day}: The dosing visit date relative to randomization.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' set.seed(314)
//'
//' f_dose_ongoing_cpp(
//'   usubjid = "A001", V = 297, C = 329, D = 569,
//'   model_ki = "zero-inflated poisson", theta_ki = c(0.4, 2.5),
//'   model_ti = "least squares", theta_ti = c(21, 2.3))
//'
//' @export
// [[Rcpp::export]]
DataFrame f_dose_ongoing_cpp(const StringVector usubjid,
                             const NumericVector V,
                             const NumericVector C,
                             const NumericVector D,
                             const std::string model_ki,
                             const NumericVector theta_ki,
                             const std::string model_ti,
                             const NumericVector theta_ti) {

  std::string mki = model_ki;
  std::for_each(mki.begin(), mki.end(), [](char & c) {
    c = std::tolower(c);
  });

  std::string mti = model_ti;
  std::for_each(mti.begin(), mti.end(), [](char & c) {
    c = std::tolower(c);
  });

  StringVector subject = Rcpp::StringVector(0);
  NumericVector day = Rcpp::NumericVector(0);

  int n = usubjid.size();
  double z, ki, Ti, Vi, u;

  for (int i=0; i<n; i++) { // loop by subject
    // initialize Vi
    Vi = V(i);

    // generate the number of skipped visits
    if (mki == "constant") {
      ki = theta_ki(0);
    } else if (mki == "poisson") {
      ki = R::rpois(theta_ki(0));
    } else if (mki == "zero-inflated poisson") {
      z = R::rbinom(1, theta_ki(0));
      if (z == 1) {
        ki = 0; // extra zeros
      } else {
        ki = R::rpois(theta_ki(1));
      }
    } else if (mki == "negative binomial") {
      ki = R::rnbinom(theta_ki(0), theta_ki(1));
    } else {
      ki = 0;
      stop("incorrect model for ki");
    }

    // gap time to the first drug dispensing visit after data cut
    if (mti == "least squares") { // draw from truncated normal distribution
      Ti = rtnormcpp((ki+1)*theta_ti(0), theta_ti(1),
                     C(i) - Vi, R_PosInf);
    } else { // draw from truncated Laplace distribution
      u = (C(i) - Vi - (ki+1)*theta_ti(0))/theta_ti(1);
      u = (u>0 ? 0.5*exp(-u) : 1-0.5*exp(u)) * R::runif(0,1);
      Ti = u<0.5 ? -log(2*u) : log(2*(1-u));
      Ti = Ti*theta_ti(1) + (ki+1)*theta_ti(0);
    }
    Ti = std::max(std::round(Ti), C(i) - Vi + 1.0);

    while (Vi + Ti < D(i)) {
      // next dispensing visit
      Vi = Vi + Ti;

      // add the new dispensing visit information
      subject.push_back(usubjid(i));
      day.push_back(Vi + 1);

      // repeat for the next dispensing visit
      if (mki == "constant") {
        ki = theta_ki(0);
      } else if (mki == "poisson") {
        ki = R::rpois(theta_ki(0));
      } else if (mki == "zero-inflated poisson") {
        z = R::rbinom(1, theta_ki(0));
        if (z == 1) {
          ki = 0; // extra zeros
        } else {
          ki = R::rpois(theta_ki(1));
        }
      } else if (mki == "negative binomial") {
        ki = R::rnbinom(theta_ki(0), theta_ki(1));
      } else {
        ki = 0;
        stop("incorrect model for ki");
      }

      // gap time to the next drug dispensing visit
      if (mti == "least squares") {
        // draw from truncated normal distribution
        Ti = rtnormcpp((ki+1)*theta_ti(0), theta_ti(1), 0, R_PosInf);
      } else { // draw from truncated Laplace distribution
        u = -(ki+1)*theta_ti(0)/theta_ti(1);
        u = (1-0.5*exp(u)) * R::runif(0,1);
        Ti = u<0.5 ? -log(2*u) : log(2*(1-u));
        Ti = Ti*theta_ti(1) + (ki+1)*theta_ti(0);
      }
      Ti = std::max(std::round(Ti), 1.0);
    }
  }

  DataFrame df = DataFrame::create(
    _["usubjid"] = subject,
    _["day"] = day);

  return df;
}



//' @title Dosing Date Imputation for New Patients
//' @description Imputes the dosing dates for new patients and ongoing
//' patients with no dosing records.
//'
//' @param usubjid The unique subject ID.
//' @param V Initialized to 0 and corresponds to the randomization visit.
//' @param C The cutoff date relative to randomization.
//' @param D The discontinuation date relative to randomization.
//' @param model_k0 The model for the number of skipped
//'   visits between randomization and the first drug dispensing visit.
//'   Options include "constant", "poisson", "zero-inflated poisson",
//'   and "negative binomial".
//' @param theta_k0 The model parameters for the number of skipped
//'   visits between randomization and the first drug dispensing visit.
//' @param model_t0 The model for the gap time between randomization
//'   and the first drug dispensing visit when there is no visit skipping.
//'   Options include "constant", "exponential", "weibull",
//'   "log-logistic", and "log-normal".
//' @param theta_t0 The model parameters for the gap time between
//'   randomization and the first drug dispensing visit when there is
//'   no visit skipping.
//' @param model_t1 The model for the gap time between randomization
//'   and the first drug dispensing visit when there is visit skipping.
//'   Options include "least squares", and "least absolute deviations".
//' @param theta_t1 The model parameters for the gap time between
//'   randomization and the first drug dispensing visit when there is
//'   visit skipping.
//' @param model_ki The model for the number of skipped
//'   visits between two consecutive drug dispensing visits.
//'   Options include "constant", "poisson", "zero-inflated poisson",
//'   and "negative binomial".
//' @param theta_ki The model parameters for the number of skipped
//'   visits between two consecutive drug dispensing visits.
//' @param model_ti The model for the gap time between two consecutive
//'   drug dispensing visits. Options include "least squares"
//'   and "least absolute deviations".
//' @param theta_ti The model parameters for the gap time between
//'   two consecutive drug dispensing visits.
//'
//' @return A data frame with two variables:
//'
//' * \code{usubjid}: The unique subject ID.
//'
//' * \code{day}: The dosing visit date relative to randomization.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' set.seed(529)
//'
//' f_dose_new_cpp(
//'   usubjid = "Z001", V = 0, C = 87, D = 985,
//'   model_k0 = "zero-inflated poisson", theta_k0 = c(0.6, 1.1),
//'   model_t0 = "log-logistic", theta_t0 = c(-1.0, 0.7),
//'   model_t1 = "least squares", theta_t1 = c(21.5, 1.9),
//'   model_ki = "zero-inflated poisson", theta_ki = c(0.1, 0.4),
//'   model_ti = "least squares", theta_ti = c(21, 2.3))
//'
//' @export
// [[Rcpp::export]]
DataFrame f_dose_new_cpp(const StringVector usubjid,
                         const NumericVector V,
                         const NumericVector C,
                         const NumericVector D,
                         const std::string model_k0,
                         const NumericVector theta_k0,
                         const std::string model_t0,
                         const NumericVector theta_t0,
                         const std::string model_t1,
                         const NumericVector theta_t1,
                         const std::string model_ki,
                         const NumericVector theta_ki,
                         const std::string model_ti,
                         const NumericVector theta_ti) {

  std::string mk0 = model_k0;
  std::for_each(mk0.begin(), mk0.end(), [](char & c) {
    c = std::tolower(c);
  });

  std::string mt0 = model_t0;
  std::for_each(mt0.begin(), mt0.end(), [](char & c) {
    c = std::tolower(c);
  });

  std::string mt1 = model_t1;
  std::for_each(mt1.begin(), mt1.end(), [](char & c) {
    c = std::tolower(c);
  });

  std::string mki = model_ki;
  std::for_each(mki.begin(), mki.end(), [](char & c) {
    c = std::tolower(c);
  });

  std::string mti = model_ti;
  std::for_each(mti.begin(), mti.end(), [](char & c) {
    c = std::tolower(c);
  });

  StringVector subject = Rcpp::StringVector(0);
  NumericVector day = Rcpp::NumericVector(0);

  int n = usubjid.size();
  double z, ki, Ti, Vi, p, u;

  for (int i=0; i<n; i++) { // loop by subject
    // initialize Vi
    Vi = V(i);

    // generate the number of skipped visits between randomization and
    // the first drug dispensing visit
    if (mk0 == "constant") {
      ki = theta_k0(0);
    } else if (mk0 == "poisson") {
      ki = R::rpois(theta_k0(0));
    } else if (mk0 == "zero-inflated poisson") {
      z = R::rbinom(1, theta_k0(0));
      if (z == 1) {
        ki = 0; // extra zeros
      } else {
        ki = R::rpois(theta_k0(1));
      }
    } else if (mk0 == "negative binomial") {
      ki = R::rnbinom(theta_k0(0), theta_k0(1));
    } else {
      ki = 0;
      stop("incorrect model for k0");
    }

    // generate the gap time between randomization and the first
    // drug dispensing visit
    if (ki == 0) {
      // simulate the gap time between randomization and the first
      // drug dispensing visit when there is no visit skipping
      if (C(i) >= 0) {
        // draw from truncated distributions for ongoing patients
        if (mt0 == "constant") {
          if (theta_t0(0) > C(i) - Vi) {
            Ti = theta_t0(0);
          } else {
            Ti = C(i) - Vi;
            stop("error in drawing T0");
          }
        } else if (mt0 == "exponential") {
          Ti = R::rexp(theta_t0(0)) + C(i) - Vi;
        } else if (mt0 == "weibull") {
          p = R::pweibull(C(i) - Vi, theta_t0(0), theta_t0(1), 0, 0);
          Ti = R::qweibull(R::runif(0,1)*p, theta_t0(0), theta_t0(1),
                           0, 0);
        } else if (mt0 == "log-logistic") {
          p = R::plogis(log(C(i) - Vi), theta_t0(0), theta_t0(1), 0, 0);
          Ti = exp(R::qlogis(R::runif(0,1)*p, theta_t0(0), theta_t0(1),
                             0, 0));
        } else if (mt0 == "log-normal") {
          Ti = exp(rtnormcpp(theta_t0(0), theta_t0(1), log(C(i) - Vi),
                             R_PosInf));
        } else {
          Ti = C(i) - Vi;
          stop("incorrect model for T0");
        }
        Ti = std::max(std::floor(Ti), C(i) - Vi + 1.0);
      } else {
        // draw from regular distributions for new patients
        if (mt0 == "constant") {
          Ti = theta_t0(0);
        } else if (mt0 == "exponential") {
          Ti = R::rexp(theta_t0(0));
        } else if (mt0 == "weibull") {
          Ti = R::rweibull(theta_t0(0), theta_t0(1));
        } else if (mt0 == "log-logistic") {
          Ti = exp(R::rlogis(theta_t0(0), theta_t0(1)));
        } else if (mt0 == "log-normal") {
          Ti = R::rlnorm(theta_t0(0), theta_t0(1));
        } else {
          Ti = 0;
          stop("incorrect model for T0");
        }
        Ti = std::max(std::floor(Ti), 0.0);
      }
    } else {
      // simulate the gap time between randomization and the first
      // drug dispensing visit when there is visit skipping
      if (C(i) >= 0) { // draw from truncated normal for ongoing patients
        if (mt1 == "least squares") {
          // draw from truncated normal distribution
          Ti = rtnormcpp(ki*theta_t1(0), theta_t1(1),
                         C(i) - Vi, R_PosInf);
        } else { // draw from truncated Laplace distribution
          u = (C(i) - Vi - ki*theta_t1(0))/theta_t1(1);
          u = (u>0 ? 0.5*exp(-u) : 1-0.5*exp(u)) * R::runif(0,1);
          Ti = u<0.5 ? -log(2*u) : log(2*(1-u));
          Ti = Ti*theta_t1(1) + ki*theta_t1(0);
        }
        Ti = std::max(std::round(Ti), C(i) - Vi + 1.0);
      } else { // draw from regular normal or Laplace for new patients
        if (mt1 == "least squares") {
          Ti = R::rnorm(ki*theta_t1(0), theta_t1(1));
        } else {
          u = R::runif(0,1);
          Ti = u<0.5 ? -log(2*u) : log(2*(1-u));
          Ti = Ti*theta_t1(1) + ki*theta_t1(0);
        }
        Ti = std::max(std::round(Ti), 1.0);
      }
    }

    // generate drug dispensing visits
    while (Vi + Ti < D(i)) {
      // next dispensing visit
      Vi = Vi + Ti;

      // add the new dispensing visit information
      subject.push_back(usubjid(i));
      day.push_back(Vi + 1);

      // generate the number of skipped visits to the next
      // drug dispensing visit
      if (mki == "constant") {
        ki = theta_ki(0);
      } else if (mki == "poisson") {
        ki = R::rpois(theta_ki(0));
      } else if (mki == "zero-inflated poisson") {
        z = R::rbinom(1, theta_ki(0));
        if (z == 1) {
          ki = 0; // extra zeros
        } else {
          ki = R::rpois(theta_ki(1));
        }
      } else if (mki == "negative binomial") {
        ki = R::rnbinom(theta_ki(0), theta_ki(1));
      } else {
        ki = 0;
        stop("incorrect model for ki");
      }

      // gap time to the next drug dispensing visit
      if (mti == "least squares") { // draw from truncated normal distribution
        Ti = rtnormcpp((ki+1)*theta_ti(0), theta_ti(1), 0, R_PosInf);
      } else { // draw from truncated Laplace distribution
        u = -(ki+1)*theta_ti(0)/theta_ti(1);
        u = (1-0.5*exp(u)) * R::runif(0,1);
        Ti = u<0.5 ? -log(2*u) : log(2*(1-u));
        Ti = Ti*theta_ti(1) + (ki+1)*theta_ti(0);
      }
      Ti = std::max(std::round(Ti), 1.0);
    }
  }

  DataFrame df = DataFrame::create(
    _["usubjid"] = subject,
    _["day"] = day);

  return df;
}


