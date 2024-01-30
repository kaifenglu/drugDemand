#' @title Cumulative Dose
#' @description Obtains the cumulative dose given treatment duration and
#' dosing schedule.
#'
#' @param x Treatment duration.
#' @param w Number of days per treatment cycle.
#' @param d Dose per treatment cycle.
#' @param N Maximum number of treatment cycles.
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
#' @param dosing_summary_t0 A data frame for the cumulative doses
#'   dispensed before the cutoff date. It contains the following
#'   variables: \code{kit}, \code{kit_name}, \code{dose_unit},
#'   and \code{cum_dose_t0}.
#' @param vf_ongoing The observed drug dispensing data for ongoing
#'   patients with drug dispensing records. It includes the following
#'   variables: \code{draw}, \code{kit}, \code{kit_name}, \code{dose_unit},
#'   \code{usubjid}, \code{day}, \code{dose}, \code{arrivalTime},
#'   \code{treatment}, \code{treatment_description},
#'   \code{time}, and \code{totalTime}.
#' @param vf_new A data frame for the randomization date for new patients
#'   and ongoing patients with no drug dispensing records.
#'   It includes the following variables:
#'   \code{draw}, \code{kit}, \code{kit_name}, \code{dose_unit},
#'   \code{usubjid}, \code{arrivalTime}, \code{treatment},
#'   \code{treatment_description}, \code{time}, and \code{totalTime}.
#' @param dosing_schedule_df A data frame providing dosing schedule
#'   information. It contains the following variables: \code{kit},
#'   \code{target_days}, \code{target_dose}, and \code{max_cycles}.
#' @param t0 The cutoff date relative to the trial start date.
#' @param t A vector of new time points for drug dispensing prediction.
#' @param pilevel The prediction interval level.
#'
#' @return A data frame for dosing summary by drug and time point per
#' protocol. It contains the following variables:
#' \code{kit}, \code{kit_name}, \code{dose_unit}, \code{t}, \code{n},
#' \code{pilevel}, \code{lower}, \code{upper}, \code{mean},
#' and \code{var}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' \donttest{
#' # Design stage drug demand forecasting per protocol.
#'
#' set.seed(312)
#' library(dplyr)
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
#' dosing_summary_t0 = drug_description_df %>%
#'   mutate(cum_dose_t0 = 0) %>%
#'   select(-c("drug", "drug_name", "dose_strength"))
#'
#' treatment_by_drug_df = f_treatment_by_drug_df(treatment_by_drug)
#'
#' vf_ongoing_new <- f_ongoing_new(
#'   newEvents, drug_description_df, treatment_by_drug_df, NULL)
#'
#' t0 = 1
#' nyears = 3
#' t1 = t0 + nyears*365
#' t = c(seq(t0, t1, 30), t1)
#' pilevel = 0.95
#'
#' dosing_pred_pp <- f_dose_pp(
#'   dosing_summary_t0, vf_ongoing_new$vf_ongoing,
#'   vf_ongoing_new$vf_new, dosing_schedule_df, t0, t, pilevel)
#'
#' head(dosing_pred_pp)
#' }
#'
#' @export
f_dose_pp <- function(dosing_summary_t0, vf_ongoing, vf_new,
                      dosing_schedule_df, t0, t, pilevel) {

  if (!is.null(vf_ongoing)) {
    df_ongoing <- vf_ongoing %>%
      group_by(.data$kit, .data$kit_name, .data$dose_unit,
               .data$draw, .data$usubjid) %>%
      slice(1) %>%
      select(-c("day", "dose"))

    if (!is.null(vf_new)) {
      df_ongoing <- df_ongoing %>%
        bind_rows(vf_new %>% filter(.data$arrivalTime <= t0))
    }
  } else {
    df_ongoing <- NULL
  }

  if (!is.null(vf_new)) {
    df_new <- vf_new %>% filter(.data$arrivalTime > t0)
  } else {
    df_new <- NULL
  }

  purrr::map_dfr(
    1:length(unique(dosing_schedule_df$kit)), function(h) {
      w = dosing_schedule_df$target_days[dosing_schedule_df$kit == h]
      d = dosing_schedule_df$target_dose[dosing_schedule_df$kit == h]
      N = dosing_schedule_df$max_cycles[dosing_schedule_df$kit == h]

      # predicted drugs to dispense for ongoing subjects
      if (!is.null(df_ongoing)) {
        dfa <- tibble(t = t) %>%
          cross_join(df_ongoing %>% filter(.data$kit == h)) %>%
          group_by(.data$kit, .data$kit_name, .data$dose_unit,
                   .data$t, .data$draw) %>%
          summarise(inc_dose = sum(
            f_cum_dose(pmin(.data$totalTime, .data$t) - .data$arrivalTime,
                       w, d, N) -
              f_cum_dose(t0 - .data$arrivalTime, w, d, N)),
            .groups = "drop_last")
      } else {
        dfa <- NULL
      }

      # predicted drugs to dispense for new subjects
      if (!is.null(df_new)) {
        dfb <- tibble(t = t) %>%
          cross_join(df_new %>% filter(.data$kit == h)) %>%
          filter(.data$arrivalTime <= t) %>%
          group_by(.data$kit, .data$kit_name, .data$dose_unit,
                   .data$t, .data$draw) %>%
          summarise(inc_dose = sum(
            f_cum_dose(pmin(.data$totalTime, .data$t) - .data$arrivalTime,
                       w, d, N)),
            .groups = "drop_last")
      } else {
        dfb <- NULL
      }

      # obtain the total doses by time point
      dosing_summary <- bind_rows(dfa, dfb) %>%
        group_by(.data$kit, .data$kit_name, .data$dose_unit,
                 .data$t, .data$draw) %>%
        summarise(inc_dose = sum(.data$inc_dose),
                  .groups = "drop_last") %>%
        left_join(dosing_summary_t0,
                  by = c("kit", "kit_name", "dose_unit")) %>%
        mutate(total_dose = .data$cum_dose_t0 + .data$inc_dose) %>%
        group_by(.data$kit, .data$kit_name, .data$dose_unit, .data$t) %>%
        summarise(n = quantile(.data$total_dose, probs = 0.5),
                  pilevel = pilevel,
                  lower = quantile(.data$total_dose, probs = (1 - pilevel)/2),
                  upper = quantile(.data$total_dose, probs = (1 + pilevel)/2),
                  mean = mean(.data$total_dose),
                  var = var(.data$total_dose),
                  .groups = "drop_last")
    })
}
