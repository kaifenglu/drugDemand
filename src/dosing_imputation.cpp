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
//' @param model_ki The model for the number of skipped visits between
//'   two consecutive drug dispensing visits.
//' @param theta_ki The model parameters for the number of skipped visits
//'   between two consecutive drug dispensing visits.
//' @param muT The regression coefficient for the linear model for the
//'   gap time between two consecutive drug dispensing visits.
//' @param sigmaT The residual standard deviation for the linear model
//'   for the gap time between two consecutive drug dispensing visits.
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
//'   model_ki = "zip", theta_ki = c(0.4, 2.5),
//'   muT = 21, sigmaT = 2.3)
//'
//' @export
// [[Rcpp::export]]
DataFrame f_dose_ongoing_cpp(const StringVector usubjid,
                             const NumericVector V,
                             const NumericVector C,
                             const NumericVector D,
                             const std::string model_ki,
                             const NumericVector theta_ki,
                             const double muT,
                             const double sigmaT) {

  StringVector subject = Rcpp::StringVector(0);
  NumericVector day = Rcpp::NumericVector(0);

  int n = usubjid.size();
  double z, ki, Ti, Vi;

  for (int i=0; i<n; i++) { // loop by subject
    // initialize Vi
    Vi = V(i);

    // generate the number of skipped visits
    if (model_ki == "constant") {
      ki = theta_ki(0);
    } else if (model_ki == "poisson") {
      ki = R::rpois(theta_ki(0));
    } else if (model_ki == "zip") {
      z = R::rbinom(1, theta_ki(0));
      if (z == 1) {
        ki = 0; // extra zeros
      } else {
        ki = R::rpois(theta_ki(1));
      }
    } else if (model_ki == "zinb") {
      z = R::rbinom(1, theta_ki(0));
      if (z == 1) {
        ki = 0; // extra zeros
      } else {
        ki = R::rnbinom(theta_ki(1), theta_ki(2));
      }
    } else {
      ki = 0;
    }

    // time from Vi to the first drug dispensing visit after data cut
    Ti = rtnormcpp((ki+1)*muT, sigmaT, C(i) - Vi, R_PosInf);
    Ti = std::max(std::round(Ti), C(i) - Vi + 1.0);

    while (Vi + Ti <= D(i)) {
      // next dispensing visit
      Vi = Vi + Ti;

      // add the new dispensing visit information
      subject.push_back(usubjid(i));
      day.push_back(Vi + 1);

      // repeat for the next dispensing visit
      if (model_ki == "constant") {
        ki = theta_ki(0);
      } else if (model_ki == "poisson") {
        ki = R::rpois(theta_ki(0));
      } else if (model_ki == "zip") {
        z = R::rbinom(1, theta_ki(0));
        if (z == 1) {
          ki = 0; // extra zeros
        } else {
          ki = R::rpois(theta_ki(1));
        }
      } else if (model_ki == "zinb") {
        z = R::rbinom(1, theta_ki(0));
        if (z == 1) {
          ki = 0; // extra zeros
        } else {
          ki = R::rnbinom(theta_ki(1), theta_ki(2));
        }
      } else {
        ki = 0;
      }

      // time from Vi to the next drug dispensing visit
      Ti = rtnormcpp((ki+1)*muT, sigmaT, 0, R_PosInf);
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
//' @param model_k0 The model for the number of skipped visits between
//'   randomization and the first drug dispensing visit.
//' @param theta_k0 The model parameters for the number of skipped
//'   visits between randomization and the first drug dispensing visit.
//' @param model_t0 The model for the gap time between randomization
//'   and the first drug dispensing visit when there is no visit skipping.
//' @param theta_t0 The model parameters for the gap time between
//'   randomization and the first drug dispensing visit when there is
//'   no visit skipping.
//' @param mu0 The regression coefficient for the linear model for the gap
//'   time between randomization and the first drug dispensing visit
//'   when there is visit skipping.
//' @param sigma0 The residual standard deviation for the linear model
//'   for the gap time between randomization and the first drug
//'   dispensing visit when there is visit skipping.
//' @param model_ki The model for the number of skipped visits between
//'   two consecutive drug dispensing visits.
//' @param theta_ki The model parameters for the number of skipped
//'   visits between two consecutive drug dispensing visits.
//' @param muT The regression coefficient for the linear model for the gap
//'   time between two consecutive drug dispensing visits.
//' @param sigmaT The residual standard deviation the linear model for the
//'   gap time between two consecutive drug dispensing visits.
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
//'   model_k0 = "zip", theta_k0 = c(0.6, 1.1),
//'   model_t0 = "log-logistic", theta_t0 = c(-1.0, 0.7),
//'   mu0 = 21.5, sigma0 = 1.9,
//'   model_ki = "zip", theta_ki = c(0.1, 0.4),
//'   muT = 21, sigmaT = 2.3)
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
                         const double mu0,
                         const double sigma0,
                         const std::string model_ki,
                         const NumericVector theta_ki,
                         const double muT,
                         const double sigmaT) {

  StringVector subject = Rcpp::StringVector(0);
  NumericVector day = Rcpp::NumericVector(0);

  int n = usubjid.size();
  double z, ki, Ti, Vi, p;

  for (int i=0; i<n; i++) { // loop by subject
    // initialize Vi
    Vi = V(i);

    // generate the number of skipped visits between randomization and
    // the first drug dispensing visit
    // repeat for the next dispensing visit
    if (model_k0 == "constant") {
      ki = theta_k0(0);
    } else if (model_k0 == "poisson") {
      ki = R::rpois(theta_k0(0));
    } else if (model_k0 == "zip") {
      z = R::rbinom(1, theta_k0(0));
      if (z == 1) {
        ki = 0; // extra zeros
      } else {
        ki = R::rpois(theta_k0(1));
      }
    } else if (model_k0 == "zinb") {
      z = R::rbinom(1, theta_k0(0));
      if (z == 1) {
        ki = 0; // extra zeros
      } else {
        ki = R::rnbinom(theta_k0(1), theta_k0(2));
      }
    } else {
      ki = 0;
    }

    // generate the gap time between randomization and the first
    // drug dispensing visit
    if (ki == 0) {
      // simulate the gap time between randomization and the first
      // drug dispensing visit in the absence of skipping visits
      if (C(i) >= 0) {
        // draw from truncated distributions for ongoing patients
        if (model_t0 == "constant") {
          if (theta_t0(0) >= C(i) - Vi) {
            Ti = theta_t0(0);
          } else {
            Ti = C(i) - Vi;
            stop("error in drawing T0");
          }
        } else if (model_t0 == "exponential") {
          Ti = R::rexp(theta_t0(0)) + C(i) - Vi;
        } else if (model_t0 == "weibull") {
          p = R::pweibull(C(i) - Vi, theta_t0(0), theta_t0(1), 0, 0);
          Ti = R::qweibull(R::runif(0,1)*p, theta_t0(0), theta_t0(1),
                           0, 0);
        } else if (model_t0 == "log-logistic") {
          p = R::plogis(log(C(i) - Vi), theta_t0(0), theta_t0(1), 0, 0);
          Ti = exp(R::qlogis(R::runif(0,1)*p, theta_t0(0), theta_t0(1),
                             0, 0));
        } else if (model_t0 == "log-normal") {
          Ti = exp(rtnormcpp(theta_t0(0), theta_t0(1), log(C(i) - Vi),
                             R_PosInf));
        } else {
          Ti = C(i) - Vi;
          stop("incorrect model for T0");
        }
        Ti = std::max(std::floor(Ti), C(i) - Vi + 1.0);
      } else {
        // draw from regular distributions for new patients
        if (model_t0 == "constant") {
          Ti = theta_t0(0);
        } else if (model_t0 == "exponential") {
          Ti = R::rexp(theta_t0(0));
        } else if (model_t0 == "weibull") {
          Ti = R::rweibull(theta_t0(0), theta_t0(1));
        } else if (model_t0 == "log-logistic") {
          Ti = exp(R::rlogis(theta_t0(0), theta_t0(1)));
        } else if (model_t0 == "log-normal") {
          Ti = exp(R::rnorm(theta_t0(0), theta_t0(1)));
        } else {
          Ti = 0;
        }
        Ti = std::max(std::floor(Ti), 0.0);
      }
    } else {
      // simulate the gap time between randomization and the first
      // drug dispensing visit in the presence of skipping visits
      if (C(i) >= 0) { // draw from truncated normal for ongoing patients
        Ti = rtnormcpp(ki*mu0, sigma0, C(i) - Vi, R_PosInf);
        Ti = std::max(std::round(Ti), C(i) - Vi + 1.0);
      } else { // draw from regular normal for new patients
        Ti = R::rnorm(ki*mu0, sigma0);
        Ti = std::max(std::round(Ti), 1.0);
      }
    }

    // generate drug dispensing visits
    while (Vi + Ti <= D(i)) {
      // next dispensing visit
      Vi = Vi + Ti;

      // add the new dispensing visit information
      subject.push_back(usubjid(i));
      day.push_back(Vi + 1);

      // repeat for the next dispensing visit
      if (model_ki == "constant") {
        ki = theta_ki(0);
      } else if (model_ki == "poisson") {
        ki = R::rpois(theta_ki(0));
      } else if (model_ki == "zip") {
        z = R::rbinom(1, theta_ki(0));
        if (z == 1) {
          ki = 0; // extra zeros
        } else {
          ki = R::rpois(theta_ki(1));
        }
      } else if (model_ki == "zinb") {
        z = R::rbinom(1, theta_ki(0));
        if (z == 1) {
          ki = 0; // extra zeros
        } else {
          ki = R::rnbinom(theta_ki(1), theta_ki(2));
        }
      } else {
        ki = 0;
      }

      // time from Vi to the next drug dispensing visit
      Ti = rtnormcpp((ki+1)*muT, sigmaT, 0, R_PosInf);
      Ti = std::max(std::round(Ti), 1.0);
    }
  }

  DataFrame df = DataFrame::create(
    _["usubjid"] = subject,
    _["day"] = day);

  return df;
}


