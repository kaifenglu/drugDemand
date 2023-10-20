#' @title Cumulative Dose
#' @description Obtains the cumulative dose given treatment duration and
#' dosing schedule.
#'
#' @param x The treatment duration.
#' @param w The number of days per treatment cycle for the drug.
#' @param d The number of kits per treatment cycle for the drug.
#' @param N The maximum number of treatment cycles for the drug.
#'
#' @return The cumulative dose to dispense for the drug over a specified
#' treatment duration.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' f_cum_dose(c(28, 70), 21, 2, 10000)
#'
#' @export
f_cum_dose <- function(x, w, d, N) {
  m = length(w)
  u = c(0, cumsum(w*N))
  i = pmin(findInterval(x, u), m)
  z = pmin(x, u[m+1]) - u[i]
  n = pmin(floor(z/w[i]) + 1, N[i])
  v = c(0, cumsum(d*N))
  v[i] + n*d[i]
}


#' @title Drug Demand Per Protocol
#' @description Obtains drug demand prediction based on protocol-assumed
#' visit and dosing schedules.
#'
#' @param dosing_summary_t0 The cumulative doses dispensed
#'   before the cutoff date. It contains the following variables:
#'   \code{drug}, \code{drug_name}, \code{dose_unit},
#'   and \code{cum_dose_t0}.
#' @param newEvents A data frame containing the imputed event data
#'   for both ongoing and new patients, typically obtained from
#'   the output of the \code{eventPred::getPrediction} function.
#'   It contains the following variables:
#'   \code{draw}, \code{usubjid}, \code{arrivalTime}, \code{treatment},
#'   \code{treatment_description}, \code{time}, \code{event},
#'   \code{dropout}, and \code{totalTime}.
#' @param treatment_by_drug_df A data frame indicating the treatments
#'   associated with each drug, including the following variables:
#'   \code{treatment}, \code{drug}, \code{drug_name}, and
#'   \code{dose_unit}.
#' @param dosing_schedule_df A data frame providing dosing schedule
#'   information. It contains the following variables: \code{drug},
#'   \code{target_days}, \code{target_kits}, and \code{max_cycles}.
#' @param t0 The cutoff date relative to the trial start date.
#' @param t A vector of new time points for drug dispensing predictions.
#' @param pilevel The prediction interval level.
#'
#' @return A data frame for dosing summary by drug and time point per
#' protocol. It contains the following variables:
#' \code{drug}, \code{drug_name}, \code{dose_unit}, \code{t}, \code{n},
#' \code{pilevel}, \code{lower}, \code{upper}, \code{mean},
#' and \code{var}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' \dontrun{
#' # Design stage drug demand predictions per protocol.
#'
#' set.seed(312)
#' library(dplyr)
#'
#' dosing_summary_t0 = drug_description_df %>%
#'   dplyr::mutate(cum_dose_t0 = 0)
#'
#' pred <- eventPred::getPrediction(
#'   df = NULL,
#'   to_predict = "enrollment and event",
#'   target_n = 250,
#'   target_d = 250,
#'   enroll_prior = list(
#'     model = "piecewise poisson",
#'     theta = c(-0.74, -1.18),
#'     vtheta = matrix(c(0.0087, 0, 0, 0.0082), 2, 2),
#'     accrualTime = c(0, 240)),
#'   event_prior = list(
#'     list(model = "log-logistic",
#'          theta = c(5.9, -0.2),
#'          vtheta = matrix(c(0.022, 0.004, 0.004, 0.012), 2, 2)),
#'     list(model = "log-logistic",
#'          theta = c(5.6, 0.02),
#'          vtheta = matrix(c(0.032, 0.003, 0.003, 0.012), 2, 2)),
#'     list(model = "log-logistic",
#'          theta = c(5.7, -0.3),
#'          vtheta = matrix(c(0.071, 0.013, 0.013, 0.054), 2, 2))),
#'   dropout_prior = NULL,
#'   pilevel = 0.95,
#'   nyears = 3,
#'   nreps = 200,
#'   showsummary = FALSE,
#'   showplot = FALSE,
#'   by_treatment = TRUE,
#'   ngroups = 3,
#'   alloc = c(2, 2, 1),
#'   treatment_label = c("Drug A + Drug B",
#'                       "Drug C + Placebo",
#'                       "Drug A + Placebo"))
#'
#' newEvents <- pred$event_pred$newEvents
#'
#' drug_name = drug_description_df$drug_name
#' dose_unit = drug_description_df$dose_unit
#' treatment_by_drug_df <- f_treatment_by_drug_df(
#'   treatment_by_drug, drug_name, dose_unit)
#'
#' t0 = 1
#' nyears = 3
#' t1 = t0 + nyears*365
#' t = c(seq(t0, t1, 30), t1)
#' pilevel = 0.95
#'
#' dosing_pred_pp <- f_dose_pp(
#'   dosing_summary_t0, newEvents,
#'   treatment_by_drug_df, dosing_schedule_df,
#'   t0, t, pilevel)
#'
#' head(dosing_pred_pp)
#' }
#'
#' @export
f_dose_pp <- function(dosing_summary_t0, newEvents,
                      treatment_by_drug_df, dosing_schedule_df,
                      t0, t, pilevel) {

  dosing_pred_pp <- dplyr::tibble()
  for (h in 1:length(unique(treatment_by_drug_df$drug))) {
    w = dosing_schedule_df$target_days[dosing_schedule_df$drug == h]
    d = dosing_schedule_df$target_kits[dosing_schedule_df$drug == h]
    N = dosing_schedule_df$max_cycles[dosing_schedule_df$drug == h]

    # treatments associated with the drug
    drug1 <- treatment_by_drug_df %>% dplyr::filter(.data$drug == h)

    # predicted event dates for ongoing patients
    df_ongoing <- newEvents %>%
      dplyr::filter(.data$arrivalTime <= t0 & .data$totalTime > t0) %>%
      dplyr::inner_join(drug1, by = "treatment") %>%
      dplyr::select(-c(.data$treatment, .data$treatment_description,
                       .data$event, .data$dropout))

    # predicted enrollment and event dates for new patients
    df_new <- newEvents %>%
      dplyr::filter(.data$arrivalTime > t0) %>%
      dplyr::inner_join(drug1, by = "treatment") %>%
      dplyr::select(-c(.data$treatment, .data$treatment_description,
                       .data$event, .data$dropout))

    # predicted drugs to dispense for ongoing subjects
    dfa = dplyr::tibble(t = t) %>%
      dplyr::cross_join(df_ongoing) %>%
      dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                      .data$t, .data$draw) %>%
      dplyr::summarise(dose_a = sum(
        f_cum_dose(pmin(.data$totalTime, .data$t) - .data$arrivalTime,
                   w, d, N) -
          f_cum_dose(t0 - .data$arrivalTime, w, d, N)),
        .groups = "drop_last")

    # predicted drugs to dispense for new subjects
    dfb = dplyr::tibble(t = t) %>%
      dplyr::cross_join(df_new) %>%
      dplyr::filter(.data$arrivalTime <= t) %>%
      dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                      .data$t, .data$draw) %>%
      dplyr::summarise(dose_b = sum(
        f_cum_dose(pmin(.data$totalTime, .data$t) - .data$arrivalTime,
                   w, d, N)),
        .groups = "drop_last")

    # obtain the total doses by time point
    dosing_summary <- dfa %>%
      dplyr::full_join(dfb, by = c("drug", "drug_name", "dose_unit",
                                   "t", "draw")) %>%
      dplyr::left_join(dosing_summary_t0, by = c('drug', 'drug_name',
                                                 'dose_unit')) %>%
      dplyr::mutate(total_dose = .data$cum_dose_t0 +
                      ifelse(is.na(.data$dose_a), 0, .data$dose_a) +
                      ifelse(is.na(.data$dose_b), 0, .data$dose_b))

    # obtain summary statistics for predicted total doses
    dosing_pred_pp1 <- dosing_summary %>%
      dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                      .data$t) %>%
      dplyr::summarise(n = quantile(.data$total_dose, probs = 0.5),
                       pilevel = pilevel,
                       lower = quantile(.data$total_dose,
                                        probs = (1 - pilevel)/2),
                       upper = quantile(.data$total_dose,
                                        probs = (1 + pilevel)/2),
                       mean = mean(.data$total_dose),
                       var = var(.data$total_dose),
                       .groups = 'drop_last')

    # combine the information for all drugs
    dosing_pred_pp <- dosing_pred_pp %>%
      dplyr::bind_rows(dosing_pred_pp1)
  }

  dosing_pred_pp
}
