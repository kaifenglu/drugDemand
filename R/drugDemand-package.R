#' @name drugDemand-package
#' @aliases drugDemand-package
#' @docType package
#'
#' @title Drug Demand Forecast
#'
#' @description Performs drug demand forecast by modeling drug
#' dispensing data along with predicted enrollment and
#' treatment discontinuation dates.
#'
#' @details In clinical trials, patients do not always follow the
#' protocol-specified visit and drug dispensing schedules.
#' Patients may be late to the drug dispensing visit.
#' Patients may skip drug dispensing visits.
#' Patients may be dispensed doses different from the target dose.
#' Predictions based solely on protocols will likely overestimate
#' the drug demand. Therefore, we propose to model the observed
#' drug dispensing data to account for the various deviations.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
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
#'   pexp plogis pnorm poisson pweibull quantile rchisq rnorm
#'   rstandard var vcov
#' @importFrom erify check_bool check_class check_content check_n
#' @importFrom rlang .data
#' @importFrom purrr map_dfr
#' @importFrom eventPred getPrediction
#'
NULL

