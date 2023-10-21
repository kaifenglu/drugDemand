#' @title Drug to Treatment Mapping
#' @description Obtains a data frame that indicates the treatment(s)
#' associated with each drug.
#'
#' @param treatment_by_drug The indicator matrix of treatment by drug
#'   combinations.
#' @param drug_name The name of the drug.
#' @param dose_unit The dose unit used for drug dispensing.
#'
#' @return A data frame indicating the treatments
#' associated with each drug, including the following variables:
#' \code{treatment}, \code{drug}, \code{drug_name}, and
#' \code{dose_unit}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' drug_name = drug_description_df$drug_name
#' dose_unit = drug_description_df$dose_unit
#' treatment_by_drug_df <- f_treatment_by_drug_df(
#'   treatment_by_drug, drug_name, dose_unit)
#' treatment_by_drug_df
#'
#' @export
f_treatment_by_drug_df <- function(
    treatment_by_drug, drug_name, dose_unit) {

  k = nrow(treatment_by_drug)
  l = ncol(treatment_by_drug)

  dplyr::tibble(treatment = rep(1:k, l),
                drug = rep(1:l, each = k),
                drug_name = rep(drug_name[1:l], each = k),
                dose_unit = rep(dose_unit[1:l], each = k),
                included = as.logical(treatment_by_drug)) %>%
    dplyr::filter(.data$included) %>%
    dplyr::select(.data$treatment, .data$drug, .data$drug_name,
                  .data$dose_unit)
}


#' @title Drug Demand Prediction
#' @description Obtains drug demand prediction via modeling and
#' simulation.
#'
#' @param df A data frame for subject-level enrollment and event data,
#'   including the following variables:
#'   \code{trialsdt}, \code{usubjid}, \code{randdt},
#'   \code{treatment}, \code{treatment_description},
#'   \code{time}, \code{event}, \code{dropout}, and \code{cutoffdt}.
#' @param newEvents A data frame containing the imputed event data
#'   for both ongoing and new patients, typically obtained from
#'   the output of the \code{eventPred::getPrediction} function.
#'   It contains the following variables:
#'   \code{draw}, \code{usubjid}, \code{arrivalTime}, \code{treatment},
#'   \code{treatment_description}, \code{time}, \code{event},
#'   \code{dropout}, and \code{totalTime}.
#' @param visitview A data frame containing the observed drug dispensing
#'   data, including the following variables:
#'   \code{usubjid}, \code{visit}, \code{date}, \code{drug},
#'   \code{drug_name}, \code{dose_unit}, \code{kit_number}, and
#'   \code{dispensed_quantity}.
#' @param drug_description_df The drug description data frame
#'   including \code{drug}, \code{drug_name}, and \code{dose_unit}.
#'   It must be specified at the design stage. It will be replaced with
#'   the observed information at the analysis stage.
#' @param treatment_by_drug The indicator matrix of treatment by drug
#'   combinations.
#' @param dosing_schedule_df A data frame providing dosing schedule
#'   information. It contains the following variables: \code{drug},
#'   \code{target_days}, \code{target_kits}, and \code{max_cycles}.
#' @param model_k0 The model for the number of skipped
#'   visits between randomization and the first drug dispensing visit.
#' @param model_t0 The model for the gap time between randomization
#'   and the first drug dispensing visit when there is no visit skipping.
#' @param model_ki The model for the number of skipped
#'   visits between two consecutive drug dispensing visits.
#' @param model_di The model for the dispensed doses at drug
#'   dispensing visits.
#' @param pilevel The prediction interval level.
#' @param nyears The number of years after the data cut for prediction.
#' @param n.cores.max The maximum number of cores to use for parallel
#'   computing. The actual number of cores used will be the minimum of
#'   \code{n.cores.max} and half of the detected number of cores.
#' @param showplot A Boolean variable that controls whether or not to
#'   show the drug dispensing model fit and drug demand prediction
#'   plots. It defaults to \code{TRUE}.
#'
#' @return A list with the following components:
#'
#' * \code{common_time_model} A Boolean variable that indicates whether
#' a common time model is used for drug dispensing visits.
#'
#' * \code{fit_k0} The model fit for the number of skipped
#' visits between randomization and the first drug dispensing visit.
#'
#' * \code{fit_t0} The model fit for the gap time between randomization
#' and the first drug dispensing visit when there is no visit skipping.
#'
#' * \code{fit_t1} The model fit for the gap time between randomization
#' and the first drug dispensing visit when there is visit skipping.
#'
#' * \code{fit_ki} The model fit for the number of skipped
#' visits between two consecutive drug dispensing visits.
#'
#' * \code{fit_ti} The model fit for the gap time between two
#' consecutive drug dispensing visits.
#'
#' * \code{fit_di} The model fit for the dispensed doses at drug
#' dispensing visits.
#'
#' * \code{dosing_subject} A data frame for the observed and imputed
#' subject-level dosing records.
#'
#' * \code{dosing_pred_df} A data frame for dosing summary by drug and
#' time point.
#'
#' * \code{dosing_pred_pp} A data frame for dosing summary by drug and
#' time point per protocol.
#'
#' * \code{dosing_pred_plot} A plot object for dosing prediction.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' \donttest{
#' set.seed(529)
#'
#' tictoc::tic("event prediction")
#'
#' pred <- eventPred::getPrediction(
#'   df = df2,
#'   to_predict = "event only",
#'   target_d = 250,
#'   event_model = "log-logistic",
#'   dropout_model = "none",
#'   pilevel = 0.95,
#'   nyears = 1,
#'   nreps = 200,
#'   showplot = FALSE,
#'   by_treatment = TRUE)
#'
#' tictoc::toc()
#'
#'
#' tictoc::tic("drug demand prediction")
#'
#' a <- f_drug_demand(
#'   df = df2,
#'   newEvents = pred$event_pred$newEvents,
#'   visitview = visitview2,
#'   treatment_by_drug = treatment_by_drug,
#'   dosing_schedule_df = dosing_schedule_df,
#'   model_k0 = "zip",
#'   model_t0 = "log-logistic",
#'   model_ki = "zip",
#'   model_di = "lme",
#'   pilevel = 0.95,
#'   nyears = 1,
#'   n.cores.max = 2,
#'   showplot = FALSE)
#'
#' tictoc::toc()
#'
#' a$dosing_pred_plot
#' }
#'
#' @export
f_drug_demand <- function(
    df = NULL,
    newEvents = NULL,
    visitview = NULL,
    drug_description_df = NULL,
    treatment_by_drug = NULL,
    dosing_schedule_df = NULL,
    model_k0 = "zip",
    model_t0 = "exponential",
    model_ki = "zip",
    model_di = "lme",
    pilevel = 0.9,
    nyears = 2,
    n.cores.max = 10,
    showplot = TRUE) {

  if (is.null(newEvents)) {
    stop("newEvents must be provided.")
  }

  nreps = length(unique(newEvents$draw))

  # trial start date and cutoff date
  if (!is.null(df)) {
    trialsdt = df$trialsdt[1]
    cutoffdt = df$cutoffdt[1]
    df <- df %>%
      dplyr::mutate(arrivalTime = as.numeric(
        .data$randdt - .data$trialsdt + 1))
  }


  # set up drug/subject/day drug dispensing data
  if (!is.null(visitview)) {
    vf <- visitview %>%
      dplyr::inner_join(df, by = "usubjid") %>%
      dplyr::mutate(day = as.numeric(.data$date - .data$randdt + 1)) %>%
      dplyr::select(.data$drug, .data$drug_name, .data$dose_unit,
                    .data$usubjid, .data$treatment,
                    .data$treatment_description, .data$arrivalTime,
                    .data$time, .data$event, .data$dropout,
                    .data$day, .data$dispensed_quantity) %>%
      dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                      .data$usubjid, .data$treatment,
                      .data$treatment_description, .data$arrivalTime,
                      .data$time, .data$event, .data$dropout,
                      .data$day) %>%
      dplyr::summarise(dose = sum(.data$dispensed_quantity),
                       .groups = "drop_last") %>%
      dplyr::mutate(cum_dose = cumsum(.data$dose)) %>%
      dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                      .data$usubjid) %>%
      dplyr::mutate(row_id = dplyr::row_number())


    # obtain the observed time points relative to trial start
    t_df <- vf %>% dplyr::mutate(t1 = .data$arrivalTime + .data$day - 1)
    t_obs <- sort(unique(t_df$t1))

    # obtain the subset of dosing records before each observed time point
    dosing_subject_t <- dplyr::tibble(t = t_obs) %>%
      dplyr::cross_join(vf) %>%
      dplyr::filter(.data$arrivalTime + .data$day - 1 <= .data$t) %>%
      dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                      .data$t, .data$usubjid) %>%
      dplyr::mutate(cum_dose = cumsum(.data$dose)) %>%
      dplyr::slice(dplyr::n())

    # tally the doses across patients
    dosing_summary_t <- dosing_subject_t %>%
      dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                      .data$t) %>%
      dplyr::summarise(n = sum(.data$cum_dose), .groups = "drop_last") %>%
      dplyr::mutate(pilevel = pilevel, lower = NA, upper = NA,
                    mean = .data$n, var = 0)

    # obtain the cumulative doses up to cutoff
    dosing_summary_t0 <- dosing_summary_t %>%
      dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit) %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::rename(cum_dose_t0 = .data$n) %>%
      dplyr::select(.data$drug, .data$drug_name, .data$dose_unit,
                    .data$cum_dose_t0)

    # extract drug description from observed data
    drug_description_df <- dosing_summary_t0 %>%
      dplyr::select(.data$drug, .data$drug_name, .data$dose_unit)
  } else {
    dosing_summary_t0 = drug_description_df %>%
      dplyr::mutate(cum_dose_t0 = 0)
  }

  # set up treatment by drug combinations
  l = nrow(drug_description_df)
  drug_name = drug_description_df$drug_name
  dose_unit = drug_description_df$dose_unit

  treatment_by_drug_df <- f_treatment_by_drug_df(
    treatment_by_drug, drug_name, dose_unit)

  # time points after cutoff to impute dosing data
  t0 = ifelse(is.null(df), 1, as.numeric(cutoffdt - trialsdt + 1))
  t1 = t0 + nyears*365
  t = c(seq(t0, t1, 30), t1)

  # dosing prediction per protocol
  dosing_pred_pp <- f_dose_pp(dosing_summary_t0, newEvents,
                              treatment_by_drug_df, dosing_schedule_df,
                              t0, t, pilevel)

  # dosing prediction based on modeling and simulation
  if (!is.null(visitview)) {
    # dosing summary for subjects who discontinued treatment before cutoff
    dosing_subject_stopped <- vf %>% dplyr::filter(.data$event == 1)

    dosing_summary_stopped <- dosing_subject_stopped %>%
      dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                      .data$usubjid) %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit) %>%
      dplyr::summarise(total_dose_a = sum(.data$cum_dose),
                       .groups = "drop_last")

    # model fit to the observed drug dispensing data
    fit <- f_dispensing_models(dosing_schedule_df$target_days, vf,
                               model_k0, model_t0, model_ki, model_di,
                               nreps, showplot)

    # impute drug dispensing data for ongoing and new patients
    a <- f_dosing_draw(df, vf, newEvents, treatment_by_drug_df,
                       fit$common_time_model,
                       fit$fit_k0, fit$fit_t0, fit$fit_t1,
                       fit$fit_ki, fit$fit_ti, fit$fit_di,
                       t0, t, n.cores.max)


    # subject level dosing data for the first simulation run
    dosing_subject <- dosing_subject_stopped %>%
      dplyr::select(-c(.data$event, .data$dropout)) %>%
      dplyr::bind_rows(
        a$dosing_subject_new %>%
          dplyr::select(-c(.data$draw, .data$totalTime)) %>%
          dplyr::group_by(.data$drug, .data$drug_name, .data$dose_unit,
                          .data$usubjid) %>%
          dplyr::mutate(cum_dose = cumsum(.data$dose),
                        row_id = dplyr::row_number())) %>%
      dplyr::mutate(
        subject_type = ifelse(
          .data$arrivalTime + .data$time - 1 <= t0, "discontinued",
          ifelse(.data$arrivalTime <= t0, "ongoing", "new")),
        imputed = ifelse(.data$arrivalTime + .data$day - 1 > t0, 1, 0)) %>%
      dplyr::arrange(.data$usubjid, .data$day,
                     .data$drug, .data$drug_name, .data$dose_unit) %>%
      dplyr::mutate(
        trialsdt = trialsdt, cutoffdt = cutoffdt,
        randdt = as.Date(.data$arrivalTime - 1, origin = trialsdt),
        adt = as.Date(.data$time - 1, origin = .data$randdt),
        date = as.Date(.data$day - 1, origin = .data$randdt))

    # dosing summary by drug, t, draw
    dosing_summary <- a$dosing_summary_new %>%
      dplyr::right_join(dosing_summary_stopped,
                        by = c("drug", "drug_name", "dose_unit")) %>%
      dplyr::mutate(total_dose = .data$total_dose_a +
                      ifelse(is.na(.data$total_dose_b), 0,
                             .data$total_dose_b))

    # dosing overview by drug, t
    dosing_overview <- dosing_summary %>%
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


    # initialize the dosing plot data set
    df0 <- dplyr::tibble(drug = 1:l, drug_name = drug_name,
                         dose_unit = dose_unit, t = 1, n = 0,
                         pilevel = pilevel, lower = NA, upper = NA,
                         mean = 0, var = 0)

    # combine with the prediction results
    dosing_pred_df <- df0 %>%
      dplyr::bind_rows(dosing_summary_t) %>%
      dplyr::bind_rows(dosing_overview) %>%
      dplyr::arrange(.data$drug, .data$drug_name, .data$dose_unit,
                     .data$t) %>%
      dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt))

    # add date information for prediction per protocol
    dosing_pred_pp <- dosing_pred_pp %>%
      dplyr::mutate(date = as.Date(.data$t - 1, origin = trialsdt))
  }


  # construct the plot
  fig <- list()
  if (is.null(df)) { # design stage
    for (j in 1:l) {
      dfb_pp <- dplyr::filter(dosing_pred_pp, .data$drug == j &
                                !is.na(.data$lower))

      fig[[j]] <- plotly::plot_ly() %>%
        plotly::add_lines(
          data = dfb_pp, x = ~t, y = ~n, name = "median prediction pp",
          line = list(width = 2)) %>%
        plotly::add_ribbons(
          data = dfb_pp, x = ~t, ymin = ~lower, ymax = ~upper,
          fill = "tonexty", line = list(width = 0),
          name = "prediction interval pp") %>%
        plotly::layout(
          xaxis = list(title = "Days since trial start", zeroline = FALSE),
          yaxis = list(title = paste0("Doses to dispense ",
                                      "(", dfb_pp$dose_unit[1], ")"),
                       zeroline = FALSE),
          legend = list(x = 0, y = 1.05, yanchor = "bottom",
                        orientation = 'h'),
          annotations = list(
            x = 0.5, y = 1,
            text = paste0("<b>", dfb_pp$drug_name[1], "</b>"),
            xanchor = "center", yanchor = "bottom",
            showarrow = FALSE, xref = 'paper', yref = 'paper'))
    }
  } else {
    for (j in 1:l) {
      dfa <- dplyr::filter(dosing_pred_df, .data$drug == j &
                             is.na(.data$lower))
      dfb <- dplyr::filter(dosing_pred_df, .data$drug == j &
                             !is.na(.data$lower))
      dfb_pp <- dplyr::filter(dosing_pred_pp, .data$drug == j &
                                !is.na(.data$lower))

      fig[[j]] <- plotly::plot_ly() %>%
        plotly::add_lines(
          data = dfa, x = ~date, y = ~n, name = "observed",
          line = list(shape ="hv", width = 2)) %>%
        plotly::add_lines(
          data = dfb, x = ~date, y = ~n, name = "median prediction",
          line = list(width = 2)) %>%
        plotly::add_ribbons(
          data = dfb, x = ~date, ymin = ~lower, ymax = ~upper,
          fill = "tonexty", line = list(width = 0),
          name = "prediction interval") %>%
        plotly::add_lines(
          data = dfb_pp, x = ~date, y = ~n, name = "median prediction pp",
          line = list(width = 2)) %>%
        plotly::add_ribbons(
          data = dfb_pp, x = ~date, ymin = ~lower, ymax = ~upper,
          fill = "tonexty", line = list(width = 0),
          name = "prediction interval pp") %>%
        plotly::add_lines(
          x = rep(cutoffdt, 2), y = c(min(dfa$n), max(dfb_pp$upper)),
          name = "cutoff", line = list(dash = "dash"),
          showlegend = FALSE) %>%
        plotly::layout(
          xaxis = list(title = "", zeroline = FALSE),
          yaxis = list(title = paste0("Doses to dispense ",
                                      "(", dfb_pp$dose_unit[1], ")"),
                       zeroline = FALSE),
          legend = list(x = 0, y = 1.05, yanchor = "bottom",
                        orientation = 'h'),
          annotations = list(
            x = 0.5, y = 1,
            text = paste0("<b>", dfb_pp$drug_name[1], "</b>"),
            xanchor = "center", yanchor = "bottom",
            showarrow = FALSE, xref = 'paper', yref = 'paper'))

      if (j==1) {
        fig[[j]] <- fig[[j]] %>%
          plotly::layout(
            annotations = list(
              x = cutoffdt, y = 0, text = 'cutoff',
              xanchor = "left", yanchor = "bottom",
              font = list(size = 12), showarrow = FALSE))
      }
    }
  }


  if (showplot) print(fig)

  # output results
  if (is.null(df)) {
    list(dosing_pred_pp = dosing_pred_pp,
         dosing_pred_plot = fig)
  } else {
    list(common_time_model = fit$common_time_model,
         fit_k0 = fit$fit_k0, fit_t0 = fit$fit_t0,
         fit_t1 = fit$fit_t1, fit_ki = fit$fit_ki,
         fit_ti = fit$fit_ti, fit_di = fit$fit_di,
         dosing_subject = dosing_subject,
         dosing_pred_df = dosing_pred_df,
         dosing_pred_pp = dosing_pred_pp,
         dosing_pred_plot = fig)
  }
}
