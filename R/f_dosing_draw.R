#' @title Drug Dispensing Visit Dates for One Iteration
#' @description Obtains drug dispensing visit dates for one iteration.
#'
#' @param i The iteration number.
#' @param fit_k0 The model fit for the number of skipped
#'   visits between randomization and the first drug dispensing visit.
#' @param fit_t0 The model fit for the gap time between randomization
#'   and the first drug dispensing visit when there is no visit skipping.
#' @param fit_t1 The model fit for the gap time between randomization
#'   and the first drug dispensing visit when there is visit skipping.
#' @param fit_ki The model fit for the number of skipped
#'   visits between two consecutive drug dispensing visits.
#' @param fit_ti The model fit for the gap time between two
#'   consecutive drug dispensing visits.
#' @param vf_ongoing1 The last observed drug dispensing date for
#'   ongoing patients with drug dispensing records, with or without
#'   the associated drug information.
#' @param vf_new1 The randomization date for new patients and ongoing
#'   patients with no drug dispensing records, with or without the
#'   associated drug information.
#'
#' @return A data frame containing the drug dispensing visit dates
#' at the subject level.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' \donttest{
#' set.seed(431)
#' library(dplyr)
#'
#' df <- df2 %>%
#'   mutate(arrivalTime = as.numeric(randdt - trialsdt + 1))
#'
#' vf <- visitview2 %>%
#'   inner_join(df, by = "usubjid") %>%
#'   mutate(day = as.numeric(date - randdt + 1)) %>%
#'   select(drug, drug_name, dose_unit, usubjid, treatment,
#'          treatment_description, arrivalTime,
#'          time, event, dropout, day, dispensed_quantity) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid, treatment,
#'            treatment_description, arrivalTime,
#'            time, event, dropout, day) %>%
#'   summarise(dose = sum(dispensed_quantity),
#'             .groups = "drop_last") %>%
#'   mutate(cum_dose = cumsum(dose)) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid) %>%
#'   mutate(row_id = row_number())
#'
#' pred <- eventPred::getPrediction(
#'   df = df,
#'   to_predict = "event only",
#'   target_d = 250,
#'   event_model = "log-logistic",
#'   dropout_model = "none",
#'   pilevel = 0.95,
#'   nyears = 3,
#'   nreps = 200,
#'   showsummary = FALSE,
#'   showplot = FALSE,
#'   by_treatment = TRUE)
#' newEvents <- pred$event_pred$newEvents
#'
#' drug_name = drug_description_df$drug_name
#' dose_unit = drug_description_df$dose_unit
#' treatment_by_drug_df <- f_treatment_by_drug_df(
#'   treatment_by_drug, drug_name, dose_unit)
#'
#' fit <- f_dispensing_models(
#'   target_days = dosing_schedule_df$target_days, vf,
#'   model_k0 = "zip", model_t0 = "log-logistic",
#'   model_ki = "zip", model_di = "lme",
#'   nreps = 200, showplot = FALSE)
#'
#' trialsdt = df$trialsdt[1]
#' cutoffdt = df$cutoffdt[1]
#' t0 = as.numeric(cutoffdt - trialsdt + 1)
#' nyears = 3
#' t1 = t0 + nyears*365
#' t = c(seq(t0, t1, 30), t1)
#'
#' nreps = length(unique(newEvents$draw))
#' l = length(unique(treatment_by_drug_df$drug))
#'
#' # all ongoing subjects
#' df_unames1 <- df %>% filter(event == 0)
#' unames1 <- df_unames1$usubjid
#'
#' # ongoing subjects with dosing records
#' df_unames2 <- vf %>% filter(event == 0) %>%
#'   group_by(usubjid) %>% slice(n())
#' unames2 <- df_unames2$usubjid
#'
#' ### dosing data for ongoing patients ###
#' vf1 <- vf %>%
#'   filter(usubjid %in% unames2) %>%
#'   select(drug, drug_name, dose_unit, usubjid, day, dose)
#'
#' # replicate for nreps times
#' vf1_rep = tibble(draw = 1:nreps) %>% cross_join(vf1)
#'
#' df1 <- newEvents %>%
#'   filter(usubjid %in% unames1) %>% select(-c(event, dropout))
#'
#' vf_ongoing <- vf1_rep %>%
#'   inner_join(df1, by = c("draw", "usubjid"))
#'
#' ### new patients and ongoing patients with no dosing records ###
#' df_new <- newEvents %>%
#'   filter(arrivalTime > t0 | usubjid %in% setdiff(unames1, unames2))
#'
#' vf_new <- purrr::map_dfr(
#'   1:l, function(h) {
#'     df_new %>%
#'       inner_join(treatment_by_drug_df %>% filter(drug == h),
#'                  by = "treatment") %>%
#'       select(-c(event, dropout))
#'   })
#'
#' # only keep the last record for each patient in each draw
#' vf_ongoing1 <- vf_ongoing %>%
#'   group_by(draw, usubjid) %>% slice(n()) %>%
#'   mutate(V = day - 1,
#'          C = as.numeric(t0 - arrivalTime),
#'          D = pmin(time - 1, t1 - arrivalTime)) %>%
#'   select(-c(drug, drug_name, dose_unit, day, dose))
#'
#' ### new patients and ongoing patients with no dosing records ###
#' vf_new1 <- vf_new %>%
#'   group_by(draw, usubjid) %>% slice(n()) %>%
#'   mutate(V = 0,
#'          C = as.numeric(t0 - arrivalTime),
#'          D = pmin(time - 1, t1 - arrivalTime)) %>%
#'   select(-c(drug, drug_name, dose_unit))
#'
#' dosing_subject_new1 <- f_dosing_draw_t_1(
#'   1, fit$fit_k0, fit$fit_t0, fit$fit_t1,
#'   fit$fit_ki, fit$fit_ti, vf_ongoing1, vf_new1)
#'
#' head(dosing_subject_new1)
#' }
#'
#' @export
f_dosing_draw_t_1 <- function(
    i, fit_k0, fit_t0, fit_t1,
    fit_ki, fit_ti,
    vf_ongoing1, vf_new1) {

  model_k0 = fit_k0$fit$model
  if (model_k0 == "constant") {
    theta_k0 = fit_k0$theta[i]
  } else if (model_k0 == "poisson") {
    theta_k0 = exp(fit_k0$theta[i])
  } else if (model_k0 == "zip") {
    theta_k0 = c(pi = plogis(fit_k0$theta[i,2]),
                 lambda = exp(fit_k0$theta[i,1]))
  } else if (model_k0 == "zinb") {
    mu = exp(fit_k0$theta[i,1])
    size = exp(fit_k0$theta[i,2])
    prob = size/(size + mu)
    theta_k0 = c(pi = plogis(fit_k0$theta[i,3]),
                 size = size,
                 prob = prob)
  }

  model_t0 = fit_t0$fit$model
  if (model_t0 == "constant") {
    theta_t0 = fit_t0$theta[i]
  } else if (model_t0 == "exponential") {
    theta_t0 = exp(fit_t0$theta[i])
  } else if (model_t0 == "weibull") {
    theta_t0 = c(kappa = exp(-fit_t0$theta[i,2]),
                 lambda = exp(fit_t0$theta[i,1]))
  } else if (model_t0 == "log-logistic") {
    theta_t0 = c(mu = fit_t0$theta[i,1],
                 sigma = exp(fit_t0$theta[i,2]))
  } else if (model_t0 == "log-normal") {
    theta_t0 = c(mu = fit_t0$theta[i,1],
                 sigma = exp(fit_t0$theta[i,2]))
  }

  mu0 = fit_t1$theta[i,1]
  sigma0 = fit_t1$theta[i,2]

  model_ki = fit_ki$fit$model
  if (model_ki == "constant") {
    theta_ki = fit_ki$theta[i]
  } else if (model_ki == "poisson") {
    theta_ki = exp(fit_ki$theta[i])
  } else if (model_ki == "zip") {
    theta_ki = c(pi = plogis(fit_ki$theta[i,2]),
                 lambda = exp(fit_ki$theta[i,1]))
  } else if (model_ki == "zinb") {
    mu = exp(fit_ki$theta[i,1])
    size = exp(fit_ki$theta[i,2])
    prob = size/(size + mu)
    theta_ki = c(pi = plogis(fit_ki$theta[i,3]),
                 size = size,
                 prob = prob)
  }

  muT = fit_ti$theta[i,1]
  sigmaT = fit_ti$theta[i,2]

  # impute dosing for ongoing patients
  df_ongoing1 <- vf_ongoing1 %>% dplyr::filter(.data$draw == i)

  # impute dosing dates for these ongoing patients
  df_ongoingi <- f_dose_ongoing_cpp(
    df_ongoing1$usubjid, df_ongoing1$V, df_ongoing1$C, df_ongoing1$D,
    model_ki, theta_ki, muT, sigmaT)

  # get other variables and combine with observed drug dispensing data
  df_ongoingi <- df_ongoingi %>%
    dplyr::left_join(df_ongoing1 %>% dplyr::select(-c(
      .data$V, .data$C, .data$D)), by = "usubjid") %>%
    dplyr::mutate(status = "ongoing")

  # impute dosing for new patients
  if (!is.null(vf_new1)) {
    # dosing data for new patients in draw i
    df_new1 <- vf_new1 %>% dplyr::filter(.data$draw == i)

    # impute dosing data for new patients
    df_newi <- f_dose_new_cpp(
      df_new1$usubjid, df_new1$V, df_new1$C, df_new1$D,
      model_k0, theta_k0, model_t0, theta_t0, mu0, sigma0,
      model_ki, theta_ki, muT, sigmaT)

    # get other variables
    df_newi <- df_newi %>%
      dplyr::left_join(df_new1 %>% dplyr::select(
        -c(.data$V, .data$C, .data$D)), by = "usubjid") %>%
      dplyr::mutate(status = "new")
  }

  # combine drug dispensing dates from ongoing and new patients
  if (is.null(vf_new1)) {  # real-time after enrollment completion
    dosing_subject_newi <- df_ongoingi
  } else { # real-time before enrollment completion
    dosing_subject_newi <- dplyr::bind_rows(df_ongoingi, df_newi)
  }

  # output drug dispensing dates
  dosing_subject_newi
}



#' @title Drug Dispensing Data for One Iteration
#' @description Obtains drug dispensing data for one iteration.
#'
#' @param i The iteration number.
#' @param common_time_model A Boolean variable that indicates whether
#'   a common time model is used for drug dispensing visits.
#' @param fit_k0 The model fit for the number of skipped
#'   visits between randomization and the first drug dispensing visit.
#' @param fit_t0 The model fit for the gap time between randomization
#'   and the first drug dispensing visit when there is no visit skipping.
#' @param fit_t1 The model fit for the gap time between randomization
#'   and the first drug dispensing visit when there is visit skipping.
#' @param fit_ki The model fit for the number of skipped
#'   visits between two consecutive drug dispensing visits.
#' @param fit_ti The model fit for the gap time between two
#'   consecutive drug dispensing visits.
#' @param fit_di The model fit for the dispensed doses at drug
#'   dispensing visits.
#' @param vf_ongoing The observed drug dispensing data for ongoing
#'   patients with drug dispensing records.
#' @param vf_ongoing1 The last observed drug dispensing date for
#'   ongoing patients with drug dispensing records, with or without
#'   the associated drug information.
#' @param vf_new1 The randomization date for new patients and ongoing
#'   patients with no drug dispensing records, with or without the
#'   associated drug information.
#' @param treatment_by_drug_df A data frame indicating the treatments
#'   associated with each drug, including the following variables:
#'   \code{treatment}, \code{drug}, \code{drug_name}, and
#'   \code{dose_unit}.
#' @param t A vector of new time points for drug dispensing predictions.
#'
#' @return A list of two components:
#'
#' * \code{dosing_subject_newi}: A data frame for the drug dispensing
#' data at the subject level by date for the given iteration.
#'
#' * \code{dosing_summary_newi}: A data frame for the drug dispensing
#' summary data by drug, time, and simulation draw for the given iteration.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' \donttest{
#' set.seed(431)
#' library(dplyr)
#'
#' df <- df2 %>%
#'   mutate(arrivalTime = as.numeric(randdt - trialsdt + 1))
#'
#' vf <- visitview2 %>%
#'   inner_join(df, by = "usubjid") %>%
#'   mutate(day = as.numeric(date - randdt + 1)) %>%
#'   select(drug, drug_name, dose_unit, usubjid, treatment,
#'          treatment_description, arrivalTime,
#'          time, event, dropout, day, dispensed_quantity) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid, treatment,
#'            treatment_description, arrivalTime,
#'            time, event, dropout, day) %>%
#'   summarise(dose = sum(dispensed_quantity),
#'             .groups = "drop_last") %>%
#'   mutate(cum_dose = cumsum(dose)) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid) %>%
#'   mutate(row_id = row_number())
#'
#' pred <- eventPred::getPrediction(
#'   df = df,
#'   to_predict = "event only",
#'   target_d = 250,
#'   event_model = "log-logistic",
#'   dropout_model = "none",
#'   pilevel = 0.95,
#'   nyears = 3,
#'   nreps = 200,
#'   showsummary = FALSE,
#'   showplot = FALSE,
#'   by_treatment = TRUE)
#' newEvents <- pred$event_pred$newEvents
#'
#' drug_name = drug_description_df$drug_name
#' dose_unit = drug_description_df$dose_unit
#' treatment_by_drug_df <- f_treatment_by_drug_df(
#'   treatment_by_drug, drug_name, dose_unit)
#'
#' fit <- f_dispensing_models(
#'   target_days = dosing_schedule_df$target_days, vf,
#'   model_k0 = "zip", model_t0 = "log-logistic",
#'   model_ki = "zip", model_di = "lme",
#'   nreps = 200, showplot = FALSE)
#'
#' trialsdt = df$trialsdt[1]
#' cutoffdt = df$cutoffdt[1]
#' t0 = as.numeric(cutoffdt - trialsdt + 1)
#' nyears = 3
#' t1 = t0 + nyears*365
#' t = c(seq(t0, t1, 30), t1)
#'
#' nreps = length(unique(newEvents$draw))
#' l = length(unique(treatment_by_drug_df$drug))
#'
#' # all ongoing subjects
#' df_unames1 <- df %>% filter(event == 0)
#' unames1 <- df_unames1$usubjid
#'
#' # ongoing subjects with dosing records
#' df_unames2 <- vf %>% filter(event == 0) %>%
#'   group_by(usubjid) %>% slice(n())
#' unames2 <- df_unames2$usubjid
#'
#' ### dosing data for ongoing patients ###
#' vf1 <- vf %>%
#'   filter(usubjid %in% unames2) %>%
#'   select(drug, drug_name, dose_unit, usubjid, day, dose)
#'
#' # replicate for nreps times
#' vf1_rep = tibble(draw = 1:nreps) %>% cross_join(vf1)
#'
#' df1 <- newEvents %>%
#'   filter(usubjid %in% unames1) %>% select(-c(event, dropout))
#'
#' vf_ongoing <- vf1_rep %>%
#'   inner_join(df1, by = c("draw", "usubjid"))
#'
#' ### new patients and ongoing patients with no dosing records ###
#' df_new <- newEvents %>%
#'   filter(arrivalTime > t0 | usubjid %in% setdiff(unames1, unames2))
#'
#' vf_new <- purrr::map_dfr(
#'   1:l, function(h) {
#'     df_new %>%
#'       inner_join(treatment_by_drug_df %>% filter(drug == h),
#'                  by = "treatment") %>%
#'       select(-c(event, dropout))
#'   })
#'
#' # only keep the last record for each patient in each draw
#' vf_ongoing1 <- vf_ongoing %>%
#'   group_by(draw, usubjid) %>% slice(n()) %>%
#'   mutate(V = day - 1,
#'          C = as.numeric(t0 - arrivalTime),
#'          D = pmin(time - 1, t1 - arrivalTime)) %>%
#'   select(-c(drug, drug_name, dose_unit, day, dose))
#'
#' ### new patients and ongoing patients with no dosing records ###
#' vf_new1 <- vf_new %>%
#'   group_by(draw, usubjid) %>% slice(n()) %>%
#'   mutate(V = 0,
#'          C = as.numeric(t0 - arrivalTime),
#'          D = pmin(time - 1, t1 - arrivalTime)) %>%
#'   select(-c(drug, drug_name, dose_unit))
#'
#' # first iteration to extract subject and summary data
#' list1 <- f_dosing_draw_1(
#'   1, fit$common_time_model,
#'   fit$fit_k0, fit$fit_t0, fit$fit_t1,
#'   fit$fit_ki, fit$fit_ti, fit$fit_di,
#'   vf_ongoing, vf_ongoing1, vf_new1,
#'   treatment_by_drug_df, t)
#'
#' head(list1$dosing_subject_newi)
#' head(list1$dosing_summary_newi)
#' }
#'
#' @export
f_dosing_draw_1 <- function(
    i, common_time_model,
    fit_k0, fit_t0, fit_t1,
    fit_ki, fit_ti, fit_di,
    vf_ongoing, vf_ongoing1, vf_new1,
    treatment_by_drug_df, t) {

  l = length(unique(treatment_by_drug_df$drug))

  # impute drug dispensing visit dates
  if (common_time_model) {
    dosing_subject_new1 <- f_dosing_draw_t_1(
      i, fit_k0, fit_t0, fit_t1,
      fit_ki, fit_ti,
      vf_ongoing1, vf_new1)

    # add drug information for each subject
    dosing_subject_new2 <- purrr::map_dfr(
      1:l, function(h) {
        dosing_subject_new1 %>%
          dplyr::inner_join(treatment_by_drug_df %>%
                              dplyr::filter(.data$drug == h),
                            by = "treatment")
      })
  } else {
    dosing_subject_new2 <- purrr::map_dfr(
      1:l, function(h) {
        f_dosing_draw_t_1(
          i, fit_k0[[h]], fit_t0[[h]], fit_t1[[h]],
          fit_ki[[h]], fit_ti[[h]],
          vf_ongoing1 %>% dplyr::filter(.data$drug == h),
          vf_new1 %>% dplyr::filter(.data$drug == h))
      })
  }


  # impute doses to dispense
  dosing_subject_new3 <- purrr::map_dfr(
    1:l, function(h) {
      mud = fit_di[[h]]$theta$fixed[i,1]
      sigmab = fit_di[[h]]$theta$fixed[i,2]
      sigmae = fit_di[[h]]$theta$fixed[i,3]
      df_ran = dplyr::tibble(usubjid = fit_di[[h]]$theta$usubjid,
                             b1 = fit_di[[h]]$theta$random[i,])

      df_ongoing2 <- dosing_subject_new2 %>%
        dplyr::filter(.data$drug == h & .data$status == "ongoing") %>%
        dplyr::inner_join(df_ran, by = "usubjid")

      df_ongoing2$dose <- pmax(round(rnorm(nrow(df_ongoing2))*sigmae +
                                       mud + df_ongoing2$b1), 1.0)

      df_ongoing2 <- df_ongoing2 %>%
        dplyr::select(-c(.data$status, .data$b1))

      df_new2 <- dosing_subject_new2 %>%
        dplyr::filter(.data$drug == h & .data$status == "new")
      n_new = nrow(df_new2)

      if (n_new > 0) {
        b1 = rnorm(n_new)*sigmab
        df_new2$dose <- pmax(round(rnorm(n_new)*sigmae + mud + b1), 1.0)
        df_new2 <- df_new2 %>% dplyr::select(-.data$status)
      }

      if (n_new == 0) {
        df_all <- df_ongoing2
      } else {
        df_all <- df_ongoing2 %>% dplyr::bind_rows(df_new2)
      }

      df_all
    })


  # add observed drug dispensing data
  dosing_subject_newi <- vf_ongoing %>%
    dplyr::filter(.data$draw == i) %>%
    dplyr::bind_rows(dosing_subject_new3) %>%
    dplyr::arrange(.data$drug, .data$drug_name, .data$dose_unit,
                   .data$draw, .data$usubjid, .data$day)

  # drug dispensed for ongoing and new subjects by drug, t, and draw
  dosing_summary_newi <- dplyr::tibble(t = t) %>%
    dplyr::cross_join(dosing_subject_newi) %>%
    dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                    .data$t, .data$draw, .data$usubjid) %>%
    dplyr::filter(.data$arrivalTime + .data$day - 1 <= .data$t) %>%
    dplyr::mutate(cum_dose = cumsum(.data$dose)) %>%
    dplyr::slice(dplyr::n()) %>%
    dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                    .data$t, .data$draw) %>%
    dplyr::summarise(total_dose_b = sum(.data$cum_dose),
                     .groups = "drop_last")

  # output drug dispensing data at the subject and summary levels
  list(dosing_subject_newi = dosing_subject_newi,
       dosing_summary_newi = dosing_summary_newi)
}



#' @title Drug Dispensing Data Simulation
#' @description Simulates drug dispensing data after cutoff for
#' both ongoing and new patients.
#'
#' @param df A data frame for subject-level enrollment and event data,
#'   including the following variables:
#'   \code{trialsdt}, \code{usubjid}, \code{randdt},
#'   \code{treatment}, \code{treatment_description},
#'   \code{time}, \code{event}, \code{dropout}, and \code{cutoffdt}.
#' @param vf A data frame for subject-level drug dispensing data,
#'   including the following variables:
#'   \code{drug}, \code{drug_name}, \code{dose_unit},
#'   \code{usubjid}, \code{treatment}, \code{treatment_description},
#'   \code{arrivalTime}, \code{time}, \code{event}, \code{dropout},
#'   \code{day}, \code{dose}, \code{cum_dose}, and \code{row_id}.
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
#' @param common_time_model A Boolean variable that indicates whether
#'   a common time model is used for drug dispensing visits.
#' @param fit_k0 The model fit for the number of skipped
#'   visits between randomization and the first drug dispensing visit.
#' @param fit_t0 The model fit for the gap time between randomization
#'   and the first drug dispensing visit when there is no visit skipping.
#' @param fit_t1 The model fit for the gap time between randomization
#'   and the first drug dispensing visit when there is visit skipping.
#' @param fit_ki The model fit for the number of skipped
#'   visits between two consecutive drug dispensing visits.
#' @param fit_ti The model fit for the gap time between two
#'   consecutive drug dispensing visits.
#' @param fit_di The model fit for the dispensed doses at drug
#'   dispensing visits.
#' @param t0 The cutoff date relative to the trial start date.
#' @param t A vector of new time points for drug dispensing predictions.
#' @param n.cores.max The maximum number of cores to use for parallel
#'   computing. The actual number of cores used will be the minimum of
#'   \code{n.cores.max} and half of the detected number of cores.
#'
#' @return A list with two components:
#'
#' * \code{dosing_subject_new} A data frame containing observed and
#' imputed subject-level dosing records for ongoing and new patients
#' for the first iteration. It contains the following variables:
#' \code{draw}, \code{drug}, \code{drug_name}, \code{dose_unit},
#' \code{usubjid}, \code{day}, \code{dose}, \code{arrivalTime},
#' \code{treatment}, \code{treatment_description}, \code{time},
#' and \code{totalTime}.
#'
#' * \code{dosing_summary_new} A data frame providing dosing summaries
#' by drug, future time point, and simulation draw for ongoing
#' and new patients. It contains the following variables:
#' \code{drug}, \code{drug_name}, \code{dose_unit}, \code{t},
#' \code{draw}, and \code{total_dose_b}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' \donttest{
#' set.seed(431)
#' library(dplyr)
#'
#' df <- df2 %>%
#'   mutate(arrivalTime = as.numeric(randdt - trialsdt + 1))
#'
#' vf <- visitview2 %>%
#'   inner_join(df, by = "usubjid") %>%
#'   mutate(day = as.numeric(date - randdt + 1)) %>%
#'   select(drug, drug_name, dose_unit, usubjid, treatment,
#'          treatment_description, arrivalTime,
#'          time, event, dropout, day, dispensed_quantity) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid, treatment,
#'            treatment_description, arrivalTime,
#'            time, event, dropout, day) %>%
#'   summarise(dose = sum(dispensed_quantity),
#'             .groups = "drop_last") %>%
#'   mutate(cum_dose = cumsum(dose)) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid) %>%
#'   mutate(row_id = row_number())
#'
#' pred <- eventPred::getPrediction(
#'   df = df,
#'   to_predict = "event only",
#'   target_d = 250,
#'   event_model = "log-logistic",
#'   dropout_model = "none",
#'   pilevel = 0.95,
#'   nyears = 3,
#'   nreps = 200,
#'   showsummary = FALSE,
#'   showplot = FALSE,
#'   by_treatment = TRUE)
#' newEvents <- pred$event_pred$newEvents
#'
#' drug_name = drug_description_df$drug_name
#' dose_unit = drug_description_df$dose_unit
#' treatment_by_drug_df <- f_treatment_by_drug_df(
#'   treatment_by_drug, drug_name, dose_unit)
#'
#' fit <- f_dispensing_models(
#'   target_days = dosing_schedule_df$target_days, vf,
#'   model_k0 = "zip", model_t0 = "log-logistic",
#'   model_ki = "zip", model_di = "lme",
#'   nreps = 200, showplot = FALSE)
#'
#' trialsdt = df$trialsdt[1]
#' cutoffdt = df$cutoffdt[1]
#' t0 = as.numeric(cutoffdt - trialsdt + 1)
#' nyears = 3
#' t1 = t0 + nyears*365
#' t = c(seq(t0, t1, 30), t1)
#'
#' a <- f_dosing_draw(
#'   df, vf, newEvents, treatment_by_drug_df,
#'   fit$common_time_model,
#'   fit$fit_k0, fit$fit_t0, fit$fit_t1,
#'   fit$fit_ki, fit$fit_ti, fit$fit_di,
#'   t0, t, n.cores.max = 2)
#'
#' head(a$dosing_subject_new)
#' head(a$dosing_summary_new)
#' }
#'
#' @export
f_dosing_draw <- function(
    df, vf, newEvents, treatment_by_drug_df,
    common_time_model,
    fit_k0, fit_t0, fit_t1,
    fit_ki, fit_ti, fit_di,
    t0, t, n.cores.max) {

  nreps = length(unique(newEvents$draw))
  l = length(unique(treatment_by_drug_df$drug))
  t1 = max(t)


  # all ongoing subjects
  df_unames1 <- df %>% dplyr::filter(.data$event == 0)
  unames1 <- df_unames1$usubjid

  # ongoing subjects with dosing records
  df_unames2 <- vf %>% dplyr::filter(.data$event == 0) %>%
    dplyr::group_by(.data$usubjid) %>%
    dplyr::slice(dplyr::n())
  unames2 <- df_unames2$usubjid

  ### dosing data for ongoing patients ###
  vf1 <- vf %>%
    dplyr::filter(.data$usubjid %in% unames2) %>%
    dplyr::select(.data$drug, .data$drug_name, .data$dose_unit,
                  .data$usubjid, .data$day, .data$dose)

  # replicate for nreps times
  vf1_rep = dplyr::tibble(draw = 1:nreps) %>% dplyr::cross_join(vf1)

  df1 <- newEvents %>%
    dplyr::filter(.data$usubjid %in% unames1) %>%
    dplyr::select(-c(.data$event, .data$dropout))

  vf_ongoing <- vf1_rep %>%
    dplyr::inner_join(df1, by = c("draw", "usubjid"))


  ### new patients and ongoing patients with no dosing records ###
  df_new <- newEvents %>%
    dplyr::filter(.data$arrivalTime > t0 |
                    .data$usubjid %in% setdiff(unames1, unames2))

  if (nrow(df_new) > 0) {
    vf_new <- purrr::map_dfr(
      1:l, function(h) {
        df_new %>%
          dplyr::inner_join(treatment_by_drug_df %>%
                              dplyr::filter(.data$drug == h),
                            by = "treatment") %>%
          dplyr::select(-c(.data$event, .data$dropout))
      })
  } else {
    vf_new <- NULL
  }


  # only keep the last record for each patient in each draw
  if (common_time_model) {
    vf_ongoing1 <- vf_ongoing %>%
      dplyr::group_by(.data$draw, .data$usubjid) %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::mutate(V = .data$day - 1,
                    C = as.numeric(t0 - .data$arrivalTime),
                    D = pmin(.data$time - 1, t1 - .data$arrivalTime)) %>%
      dplyr::select(-c(.data$drug, .data$drug_name, .data$dose_unit,
                       .data$day, .data$dose))

    ### new patients and ongoing patients with no dosing records ###
    if (!is.null(vf_new)) {
      vf_new1 <- vf_new %>%
        dplyr::group_by(.data$draw, .data$usubjid) %>%
        dplyr::slice(dplyr::n()) %>%
        dplyr::mutate(V = 0,
                      C = as.numeric(t0 - .data$arrivalTime),
                      D = pmin(.data$time - 1, t1 - .data$arrivalTime)) %>%
        dplyr::select(-c(.data$drug, .data$drug_name, .data$dose_unit))
    } else {
      vf_new1 <- NULL
    }
  } else {
    vf_ongoing1 <- vf_ongoing %>%
      dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                      .data$draw, .data$usubjid) %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::mutate(V = .data$day - 1,
                    C = as.numeric(t0 - .data$arrivalTime),
                    D = pmin(.data$time - 1, t1 - .data$arrivalTime)) %>%
      dplyr::select(-c(.data$day, .data$dose))

    if (!is.null(vf_new)) {
      vf_new1 <- vf_new %>%
        dplyr::mutate(V = 0,
                      C = as.numeric(t0 - .data$arrivalTime),
                      D = pmin(.data$time - 1, t1 - .data$arrivalTime))
    } else {
      vf_new1 <- NULL
    }
  }


  # first iteration to extract subject and summary data
  i = 1
  list1 <- f_dosing_draw_1(
    i, common_time_model,
    fit_k0, fit_t0, fit_t1,
    fit_ki, fit_ti, fit_di,
    vf_ongoing, vf_ongoing1, vf_new1,
    treatment_by_drug_df, t)

  dosing_subject_new <- list1$dosing_subject_newi


  # register parallel backend
  n.cores <- min(n.cores.max, parallel::detectCores()/2)
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)

  # subsequent iterations to extract summary data only
  dosing_summary_new <- foreach::foreach(
    i = 2:nreps, .combine = "bind_rows",
    .packages = c("dplyr", "mvtnorm")
  ) %dorng% {
    f_dosing_draw_1(
      i, common_time_model,
      fit_k0, fit_t0, fit_t1,
      fit_ki, fit_ti, fit_di,
      vf_ongoing, vf_ongoing1, vf_new1,
      treatment_by_drug_df, t)$dosing_summary_newi
  }

  # shut down the cluster of workers
  parallel::stopCluster(cl)


  # combine the summary data from all iterations
  dosing_summary_new <- list1$dosing_summary_newi %>%
    dplyr::bind_rows(dosing_summary_new)


  # output results for ongoing and new patients
  list(dosing_subject_new = dosing_subject_new,
       dosing_summary_new = dosing_summary_new)
}

