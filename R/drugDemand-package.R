#' @name drugDemand-package
#' @aliases drugDemand-package
#' @docType package
#'
#' @title Drug Demand Forecasting
#'
#' @description Performs drug demand forecasting by modeling drug
#' dispensing data while taking into account predicted enrollment
#' and treatment discontinuation dates. The gap time between
#' randomization and the first drug dispensing visit is modeled
#' using interval-censored exponential, Weibull, log-logistic, or
#' log-normal distributions
#' (Anderson-Bergman (2017) <doi:10.18637/jss.v081.i12>).
#' The number of skipped visits is modeled using zero-inflated
#' Poisson or zero-inflated negative binomial distributions
#' (Zeileis, Kleiber & Jackman (2008) <doi:10.18637/jss.v027.i08>).
#' The gap time between two consecutive drug dispensing visis
#' is modeled using linear regression given the number of skipped
#' visits. The number of dispensed doses is modeled using a linear
#' mixed-effects model
#' (McCulloch & Searle (2001, ISBN:0-471-19364-X)).
#'
#' @details In clinical trials, patients do not always follow
#' protocol-specified visit and drug dispensing schedules.
#' Patients may encounter delays in their drug dispensing
#' appointments, skip visits altogether, or receive doses
#' different from the protocol-specified target.
#' Relying solely on protocol-based predictions tends to result
#' in an overestimation of drug demand. Consequently, we propose
#' a method that models observed drug dispensing data,
#' thereby accounting for these deviations.
#'
#' * \code{k0} The number of skipped visits between randomization
#' and the first drug dispensing visit.
#'
#' * \code{t0} The time elapsed between randomization and the first
#' drug dispensing visit when k0 equals 0.
#'
#' * \code{t1} The time elapsed between randomization and the first
#' drug dispensing visit when k0 is greater than 0.
#'
#' * \code{ki} The number of skipped visits between two consecutive
#' drug dispensing visits.
#'
#' * \code{ti} The time elapsed between two consecutive drug
#' dispensing visits.
#'
#' * \code{di} The number of kits dispensed at drug dispensing visits.
#'
#' For \code{k0} and \code{ki}, we explore several modeling options,
#' including constant, Poisson, zero-inflated Poisson (ZIP), and
#' zero-inflated negative binomial (ZINB) distributions.
#'
#' For \code{t0}, we consider various models such as constant,
#' exponential, Weibull, log-logistic, and log-normal.
#'
#' Linear regression models are applied to \code{t1} (given \code{k0})
#' and \code{ti} (given \code{ki}).
#'
#' For \code{di}, we evaluate constant, linear, and linear
#' mixed-effects models with subject random effects.
#'
#' Once the dosing models are fitted to the observed drug
#' dispensing data, we draw model parameters from their
#' approximate posterior distributions. Subsequently, we simulate
#' drug dispensing data after cutoff for both ongoing and new patients.
#'
#' Finally, we estimate the number of kits to dispense based on the
#' simulated data.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#'
#' Clifford Anderson-Bergman.
#' icenReg: Regression Models for Interval Censored Data in R.
#' J Stat Softw. 2017, Volume 81, Issue 12.
#'
#' Achim Zeileis, Christian Kleiber, and Simon Jackman.
#' Regression models for count data in R.
#' J Stat Softw. 2008, Volume 27, Issue 8.
#'
#' Charles E. McCulloch and Shayler R. Searle.
#' Generalized, Linear, and Mixed Models.
#' John Wiley & Sons: New York, 2001,
#' ISBN:0-471-19364-X
#'
#' @useDynLib drugDemand, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom mvtnorm pmvnorm rmvnorm
#' @importFrom dplyr %>% arrange as_tibble bind_cols bind_rows
#'   cross_join filter group_by inner_join lead left_join
#'   mutate n rename row_number select slice summarise tibble
#' @importFrom plotly add_trace layout plot_ly
#' @importFrom icenReg ic_par
#' @importFrom pscl zeroinfl
#' @importFrom nlme lme
#' @importFrom parallel detectCores makeCluster
#' @importFrom foreach %do% %dopar% foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG %dorng%
#' @importFrom survival Surv survfit
#' @importFrom stats AIC BIC dnbinom dpois glm lm logLik optimHess
#'   pexp plogis plnorm pnorm poisson pweibull quantile rchisq
#'   rnorm rstandard var vcov
#' @importFrom erify check_bool check_class check_content check_n
#' @importFrom rlang .data
#' @importFrom purrr map_dfr
#' @importFrom eventPred getPrediction
#' @importFrom tictoc tic toc
#'
NULL

