#' @title Drug to Treatment Mapping
#' @description Obtains a data frame that indicates the treatments
#' associated with each drug.
#'
#' @param treatment_by_drug The indicator matrix of treatment by drug
#'   combinations.
#'
#' @return A data frame indicating the treatments
#' associated with each drug, including the following variables:
#' \code{treatment} and \code{drug}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' f_treatment_by_drug_df(treatment_by_drug)
#'
#' @export
f_treatment_by_drug_df <- function(treatment_by_drug) {
  k = nrow(treatment_by_drug)
  l = ncol(treatment_by_drug)
  tibble(treatment = rep(1:k, l),
         drug = rep(1:l, each = k),
         included = as.logical(treatment_by_drug)) %>%
    filter(.data$included) %>%
    select(c("treatment", "drug"))
}


#' @title Random Number Generator for the Dirichlet Distribution
#' @description Generates cell probabilities from the Dirichlet distribution.
#'
#' @param n The number of observations.
#' @param alpha The shape parameters of the Dirichlet distribution.
#'
#' @return A matrix of n rows and k columns, where n is the number of
#' observations and k is the number of cells.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' rdirichlet(2, c(50, 20, 30))
#'
#' @export
rdirichlet <- function (n = 1, alpha) {
  Gam <- matrix(0, n, length(alpha))
  for (i in 1:length(alpha)) Gam[,i] <- rgamma(n, shape = alpha[i])
  Gam/rowSums(Gam)
}


#' @title Drug Demand Forecasting
#' @description Obtains drug demand forecasting via modeling and
#' simulation.
#'
#' @param df A data frame for subject-level enrollment and event data,
#'   including the following variables:
#'   \code{trialsdt}, \code{usubjid}, \code{randdt},
#'   \code{treatment}, \code{treatment_description},
#'   \code{time}, \code{event}, \code{dropout}, and \code{cutoffdt}.
#' @param newEvents A data frame containing the imputed event data
#'   for both ongoing and new patients, typically obtained from
#'   the output of the \code{getPrediction} function of the
#'   \code{eventPred} package. It contains the following variables:
#'   \code{draw}, \code{usubjid}, \code{arrivalTime}, \code{treatment},
#'   \code{treatment_description}, \code{time}, \code{event},
#'   \code{dropout}, and \code{totalTime}.
#' @param visitview A data frame containing the observed drug dispensing
#'   data, including the following variables:
#'   \code{usubjid}, \code{visit}, \code{date}, \code{drug},
#'   \code{drug_name}, \code{dose_strength}, \code{dose_unit},
#'   \code{kit_number}, and \code{dispensed_quantity}.
#' @param drug_description_df The drug description data frame
#'   including \code{drug}, \code{drug_name}, \code{dose_strength},
#'   \code{kit}, \code{kit_name}, and \code{dose_unit}.
#'   It must be specified at the design stage. It will be replaced with
#'   the observed information at the analysis stage.
#' @param treatment_by_drug The indicator matrix of treatment by drug
#'   combinations. It must be specified at the design stage. It will
#'   be replaced with the observed information at the analysis stage.
#' @param dosing_schedule_df A data frame providing dosing schedule
#'   information. It contains the following variables: \code{kit},
#'   \code{target_days}, \code{target_dose}, and \code{max_cycles}.
#' @param model_k0 The model for the number of skipped
#'   visits between randomization and the first drug dispensing visit.
#'   Options include "constant", "poisson", "zero-inflated poisson",
#'   and "negative binomial".
#' @param model_t0 The model for the gap time between randomization
#'   and the first drug dispensing visit when there is no visit skipping.
#'   Options include "constant", "exponential", "weibull",
#'   "log-logistic", and "log-normal".
#' @param model_t1 The model for the gap time between randomization
#'   and the first drug dispensing visit when there is visit skipping.
#'   Options include "least squares" and "least absolute deviations".
#' @param model_ki The model for the number of skipped
#'   visits between two consecutive drug dispensing visits.
#'   Options include "constant", "poisson", "zero-inflated poisson",
#'   and "negative binomial".
#' @param model_ti The model for the gap time between two consecutive
#'   drug dispensing visits. Options include "least squares"
#'   and "least absolute deviations".
#' @param model_di The model for the dispensed doses at drug
#'   dispensing visits. Options include "constant",
#'   "linear model", and "linear mixed-effects model".
#' @param pilevel The prediction interval level.
#' @param nyears The number of years after the data cut for prediction.
#' @param ncores_max The maximum number of cores to use for parallel
#'   computing. The actual number of cores used is the minimum of
#'   \code{ncores_max} and half of the detected number of cores.
#' @param pred_pp_only A Boolean variable that controls whether or not
#'   to make protocol-based predictions only.
#' @param showplot A Boolean variable that controls whether or not to
#'   show the drug dispensing model fit and drug demand prediction
#'   plots. It defaults to \code{TRUE}.
#'
#' @return For design-stage drug demand forecasting, a list with the
#' following components:
#'
#' * \code{dosing_pred_df}: A data frame for dosing summary by kit type
#'    and time point per protocol. It includes the following variables:
#'   \code{kit}, \code{kit_name}, \code{dose_unit},
#'   \code{t}, \code{n}, \code{pilevel},
#'   \code{lower}, \code{upper}, \code{mean}, \code{var},
#'   and \code{parameter}.
#'
#' * \code{dosing_pred_plot}: A plot object for dosing prediction.
#'
#' For analysis-stage drug demand forecasting, a list with the
#' following components:
#'
#' * \code{trialsdt}: The trial start date.
#'
#' * \code{cutoffdt}: The cutoff date.
#'
#' * \code{dosing_summary_t0}: A data frame for the cumulative doses
#'   dispensed before the cutoff date. It contains the following
#'   variables: \code{kit}, \code{kit_name}, \code{dose_unit}, and
#'   \code{cum_dose_t0}.
#'
#' * \code{cum_dispense_plot}: The step plot for the cumulative doses
#'   dispensed for each kit type.
#'
#' * \code{bar_t0_plot}: The bar chart for the time between
#'   randomization and the first drug dispensing visit.
#'
#' * \code{bar_ti_plot}: The bar chart for the gap time between two
#'   consecutive drug dispensing visits.
#'
#' * \code{bar_di_plot}: The bar chart for the doses dispensed at drug
#'   dispensing visits.
#'
#' * \code{common_time_model}: A Boolean variable that indicates
#'   whether a common time model is used for drug dispensing visits.
#'
#' * \code{k0_fit}: The model fit for the number of skipped
#'   visits between randomization and the first drug dispensing visit.
#'
#' * \code{t0_fit}: The model fit for the gap time between
#'   randomization and the first drug dispensing visit when there is
#'   no visit skipping.
#'
#' * \code{t1_fit}: The model fit for the gap time between
#'   randomization and the first drug dispensing visit when there is
#'   visit skipping.
#'
#' * \code{ki_fit}: The model fit for the number of skipped
#'   visits between two consecutive drug dispensing visits.
#'
#' * \code{ti_fit}: The model fit for the gap time between two
#'   consecutive drug dispensing visits.
#'
#' * \code{di_fit}: The model fit for the dispensed doses at drug
#'   dispensing visits.
#'
#' * \code{dosing_subject}: A data frame for the observed and imputed
#'   subject-level dosing records for the first iteration. It includes
#'   the following variables:
#'   \code{kit}, \code{kit_name}, \code{dose_unit},
#'   \code{usubjid}, \code{treatment}, \code{treatment_description},
#'   \code{arrivalTime}, \code{time}, \code{day}, \code{dose},
#'   \code{cum_dose}, \code{row_id}, \code{subject_type}, \code{imputed},
#'   \code{trialsdt}, \code{cutoffdt}, \code{randdt}, \code{adt}, and
#'   \code{date}.
#'
#' * \code{dosing_pred_df}: A data frame for dosing summary by drug and
#'   time point. It includes the following variables:
#'   \code{kit}, \code{kit_name}, \code{dose_unit},
#'   \code{t}, \code{n}, \code{pilevel},
#'   \code{lower}, \code{upper}, \code{mean}, \code{var}, \code{date},
#'   and \code{parameter}.
#'
#' * \code{dosing_pred_plot}: A plot object for dosing prediction.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @seealso \code{\link{f_fit_t0}}, \code{\link{f_fit_ki}},
#' \code{\link{f_fit_ti}}, \code{\link{f_fit_di}}
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
#' drug_demand <- f_drug_demand(
#'   df = df2,
#'   newEvents = pred$event_pred$newEvents,
#'   visitview = visitview2,
#'   dosing_schedule_df = dosing_schedule_df,
#'   model_k0 = "zero-inflated poisson",
#'   model_t0 = "log-logistic",
#'   model_t1 = "least squares",
#'   model_ki = "zero-inflated poisson",
#'   model_ti = "least squares",
#'   model_di = "linear mixed-effects model",
#'   pilevel = 0.95,
#'   nyears = 1,
#'   ncores_max = 2,
#'   showplot = FALSE)
#'
#' tictoc::toc()
#'
#' drug_demand$dosing_pred_plot
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
    model_k0 = "negative binomial",
    model_t0 = "log-logistic",
    model_t1 = "least squares",
    model_ki = "negative binomial",
    model_ti = "least absolute deviations",
    model_di = "linear mixed-effects model",
    pilevel = 0.95,
    nyears = 1,
    ncores_max = 10,
    pred_pp_only = FALSE,
    showplot = TRUE) {

  # check if df has the required columns
  if (!is.null(df)) {
    cols = colnames(df)
    req_cols = c("trialsdt", "usubjid", "randdt", "treatment",
                 "treatment_description", "time", "event",
                 "dropout", "cutoffdt")

    if (!all(req_cols %in% cols)) {
      stop(paste("The following columns are missing from df:",
                 paste(req_cols[!(req_cols %in% cols)], collapse = ", ")))
    }

    if (any(is.na(df[, req_cols]))) {
      stop(paste("The following columns of df have missing values:",
                 paste(req_cols[sapply(df, function(x) any(is.na(x)))],
                       collapse = ", ")))
    }
  }

  if (is.null(newEvents)) {
    stop("newEvents must be provided.")
  }

  if (!("treatment" %in% colnames(newEvents))) {
    newEvents$treatment = 1
    newEvents$treatment_description = "Treatment 1"
  }

  # check if visitview has the required columns
  if (!is.null(visitview)) {
    cols = colnames(visitview)
    req_cols = c("usubjid", "date", "drug", "drug_name", "dose_strength",
                 "dose_unit", "dispensed_quantity")

    if (!all(req_cols %in% cols)) {
      stop(paste("The following columns are missing from visitview:",
                 paste(req_cols[!(req_cols %in% cols)], collapse = ", ")))
    }

    if (any(is.na(visitview[, req_cols]))) {
      stop(paste("The following columns of visitview have missing values:",
                 paste(req_cols[sapply(visitview,
                                       function(x) any(is.na(x)))],
                       collapse = ", ")))
    }
  }

  if (is.null(df) + is.null(visitview) == 1) {
    stop("df and visitview must be provided at the same time.")
  }

  if (is.null(df) && is.null(drug_description_df)) {
    stop("drug_description_df must be provided if df is not provided.")
  }

  if (is.null(df) && is.null(treatment_by_drug)) {
    stop("treatment_by_drug must be provided if df is not provided.")
  }


  if (is.null(dosing_schedule_df)) {
    stop("dosing_schedule_df must be provided.")
  }

  # check if drug_description_df has the required columns
  if (!is.null(drug_description_df)) {
    cols = colnames(drug_description_df)
    req_cols = c("drug", "drug_name", "dose_strength", "dose_unit")

    if (is.null(df)) req_cols = c(req_cols, "p_kit")

    if (!all(req_cols %in% cols)) {
      stop(paste("The following columns are missing from",
                 "drug_description_df:",
                 paste(req_cols[!(req_cols %in% cols)], collapse = ", ")))
    }

    if (any(is.na(drug_description_df[, req_cols]))) {
      stop(paste("The following columns of drug_description_df",
                 "have missing values:",
                 paste(req_cols[sapply(drug_description_df,
                                       function(x) any(is.na(x)))],
                       collapse = ", ")))
    }

    # ensure that the kit probabilities add up to 1 for each drug
    if ("p_kit" %in% cols) {
      if (any(drug_description_df$p_kit <= 0)) {
        stop("p_kit must be positive")
      }

      df_total_p_kit <- drug_description_df %>%
        group_by(.data$drug) %>%
        summarise(total_p_kit = sum(.data$p_kit))

      drug_description_df <- drug_description_df %>%
        left_join(df_total_p_kit, by = "drug") %>%
        mutate(p_kit = .data$p_kit/.data$total_p_kit) %>%
        select(-.data$total_p_kit)
    }

    # derive the kit information
    drug_description_df <- drug_description_df %>%
      arrange(.data$drug, .data$drug_name, .data$dose_strength,
              .data$dose_unit) %>%
      mutate(kit = row_number(),
             kit_name = trimws(paste(.data$drug_name, .data$dose_strength)))
  }

  nreps = length(unique(newEvents$draw))

  # trial start date and cutoff date
  if (!is.null(df)) {
    trialsdt = df$trialsdt[1]
    cutoffdt = df$cutoffdt[1]
  }

  # set up drug/subject/day drug dispensing data
  if (!is.null(visitview)) {
    observed <- f_dose_observed(df, visitview, showplot = showplot)
    vf = observed$vf
    treatment_by_drug_df = observed$treatment_by_drug_df
    drug_description_df = observed$drug_description_df
    dosing_summary_t = observed$dosing_summary_t %>%
      mutate(pilevel = pilevel)
    dosing_summary_t0 = observed$dosing_summary_t0
  } else {
    vf <- NULL
    treatment_by_drug_df = f_treatment_by_drug_df(treatment_by_drug)
    dosing_summary_t0 = drug_description_df %>%
      mutate(cum_dose_t0 = 0) %>%
      select(-c("drug", "drug_name", "dose_strength"))
  }

  # prepare the dosing data sets to impute for ongoing and new subjects
  vf_ongoing_new <- f_ongoing_new(newEvents, drug_description_df,
                                  treatment_by_drug_df, vf)

  vf_ongoing <- vf_ongoing_new$vf_ongoing
  vf_new <- vf_ongoing_new$vf_new

  # number of kit types
  l = nrow(drug_description_df)
  kit_name = drug_description_df$kit_name
  dose_unit = drug_description_df$dose_unit

  # time points after cutoff to impute dosing data
  t0 = ifelse(is.null(df), 1, as.numeric(cutoffdt - trialsdt + 1))
  t1 = t0 + nyears*365
  t = c(seq(t0, t1, 30), t1)

  # dosing prediction per protocol
  dosing_pred_pp <- f_dose_pp(dosing_summary_t0, vf_ongoing, vf_new,
                              dosing_schedule_df, t0, t, pilevel) %>%
    mutate(parameter = "protocol based prediction")

  # dosing prediction based on modeling and simulation
  if (!is.null(visitview)) {
    # dosing summary for subjects who discontinued treatment before cutoff
    dosing_subject_stopped <- vf %>% filter(.data$event == 1)

    dosing_summary_stopped <- dosing_subject_stopped %>%
      group_by(.data$kit, .data$kit_name, .data$dose_unit,
               .data$usubjid) %>%
      slice(n()) %>%
      group_by(.data$kit, .data$kit_name, .data$dose_unit) %>%
      summarise(total_dose_a = sum(.data$cum_dose),
                .groups = "drop_last")

    # add arms for which all patients had discontinued treatment before cutoff
    dosing_pred_pp <- dosing_summary_stopped %>%
      cross_join(tibble(t = t, pilevel = pilevel,
                        parameter = "protocol based prediction")) %>%
      left_join(dosing_pred_pp,
                by = c("kit", "kit_name", "dose_unit", "t",
                       "pilevel", "parameter")) %>%
      mutate(miss = is.na(.data$n)) %>%
      mutate(n = ifelse(.data$miss, .data$total_dose_a, .data$n),
             lower = ifelse(.data$miss, .data$total_dose_a, .data$lower),
             upper = ifelse(.data$miss, .data$total_dose_a, .data$upper),
             mean = ifelse(.data$miss, .data$total_dose_a, .data$mean),
             var = ifelse(.data$miss, 0, .data$var)) %>%
      mutate(date = as.Date(.data$t - 1, origin = trialsdt)) %>%
      select(-c("total_dose_a", "miss"))

    # initialize the dosing plot data set
    df0 <- tibble(kit = 1:l, kit_name = kit_name,
                  dose_unit = dose_unit, t = 1, n = 0,
                  pilevel = pilevel, lower = NA, upper = NA,
                  mean = 0, var = 0)

    if (!pred_pp_only) {
      # model fit to the observed drug dispensing data
      fit <- f_dispensing_models(vf, dosing_schedule_df,
                                 model_k0, model_t0, model_t1,
                                 model_ki, model_ti, model_di,
                                 nreps, showplot)

      # impute drug dispensing data for ongoing and new patients
      a <- f_dose_draw(vf_ongoing, vf_new,
                       fit$common_time_model,
                       fit$k0_fit, fit$t0_fit, fit$t1_fit,
                       fit$ki_fit, fit$ti_fit, fit$di_fit,
                       t0, t, ncores_max)

      # subject level dosing data for the first simulation run
      dosing_subject <- dosing_subject_stopped %>%
        select(-c("event", "dropout")) %>%
        bind_rows(
          a$dosing_subject_new %>%
            select(-c("draw", "totalTime")) %>%
            group_by(.data$kit, .data$kit_name, .data$dose_unit,
                     .data$usubjid) %>%
            mutate(cum_dose = cumsum(.data$dose),
                   row_id = row_number())) %>%
        mutate(
          subject_type = ifelse(
            .data$arrivalTime + .data$time - 1 <= t0, "discontinued",
            ifelse(.data$arrivalTime <= t0, "ongoing", "new")),
          imputed = ifelse(.data$arrivalTime + .data$day - 1 > t0, 1, 0)) %>%
        arrange(.data$kit, .data$kit_name, .data$dose_unit,
                .data$usubjid, .data$day) %>%
        mutate(
          trialsdt = trialsdt,
          cutoffdt = cutoffdt,
          randdt = as.Date(.data$arrivalTime - 1, origin = trialsdt),
          adt = as.Date(.data$time - 1, origin = .data$randdt),
          date = as.Date(.data$day - 1, origin = .data$randdt))

      # set total_dose_b = 0 for treatment arms for which
      # all patients had discontinued treatment before cutoff
      dosing_summary_new <- drug_description_df %>%
        cross_join(tibble(draw = rep(1:nreps, each = length(t)),
                          t = rep(t, nreps))) %>%
        left_join(a$dosing_summary_new,
                  by = c("kit", "kit_name", "dose_unit",
                         "draw", "t")) %>%
        mutate(total_dose_b = ifelse(is.na(.data$total_dose_b), 0,
                                     .data$total_dose_b))

      # add summary from patients who had discontinued treatment before cutoff
      dosing_summary <- dosing_summary_new %>%
        right_join(dosing_summary_stopped,
                   by = c("kit", "kit_name", "dose_unit")) %>%
        mutate(total_dose = .data$total_dose_a +
                 ifelse(is.na(.data$total_dose_b), 0,
                        .data$total_dose_b))

      # dosing overview by drug, t
      dosing_overview <- dosing_summary %>%
        group_by(.data$kit, .data$kit_name, .data$dose_unit,
                 .data$t) %>%
        summarise(n = quantile(.data$total_dose, probs = 0.5),
                  pilevel = pilevel,
                  lower = quantile(.data$total_dose,
                                   probs = (1 - pilevel)/2),
                  upper = quantile(.data$total_dose,
                                   probs = (1 + pilevel)/2),
                  mean = mean(.data$total_dose),
                  var = var(.data$total_dose),
                  .groups = 'drop_last')

      # combine with the prediction results
      dosing_pred_df <- df0 %>%
        bind_rows(dosing_summary_t) %>%
        bind_rows(dosing_overview) %>%
        arrange(.data$kit, .data$kit_name, .data$dose_unit,
                .data$t) %>%
        mutate(date = as.Date(.data$t - 1, origin = trialsdt))
    } else {
      dosing_subject <- vf %>%
        mutate(
          subject_type = ifelse(.data$event == 1, "discontinued", "ongoing"),
          imputed = 0) %>%
        select(-c("event", "dropout")) %>%
        arrange(.data$kit, .data$kit_name, .data$dose_unit,
                .data$usubjid, .data$day) %>%
        mutate(
          trialsdt = trialsdt,
          cutoffdt = cutoffdt,
          randdt = as.Date(.data$arrivalTime - 1, origin = trialsdt),
          adt = as.Date(.data$time - 1, origin = .data$randdt),
          date = as.Date(.data$day - 1, origin = .data$randdt))

      dosing_pred_df <- df0 %>%
        bind_rows(dosing_summary_t) %>%
        arrange(.data$kit, .data$kit_name, .data$dose_unit,
                .data$t) %>%
        mutate(date = as.Date(.data$t - 1, origin = trialsdt))
    }

    dosing_pred_df <- dosing_pred_df %>%
      mutate(parameter = ifelse(is.na(.data$lower), "observed data",
                                "model based prediction")) %>%
      bind_rows(dosing_pred_pp) %>%
      arrange(.data$kit, .data$kit_name, .data$dose_unit, .data$t)
  } else {
    dosing_pred_df <- dosing_pred_pp
  }


  # construct the plot
  fig <- list()
  if (is.null(df)) { # design stage
    for (j in 1:l) {
      dfb_pp <- filter(dosing_pred_df, .data$kit == j &
                         .data$parameter ==
                         "protocol based prediction")

      fig[[j]] <- plotly::plot_ly() %>%
        plotly::add_lines(
          data = dfb_pp, x = ~t, y = ~n, name = "median prediction protocol",
          line = list(width = 2)) %>%
        plotly::add_ribbons(
          data = dfb_pp, x = ~t, ymin = ~lower, ymax = ~upper,
          fill = "tonexty", line = list(width = 0),
          name = "prediction interval protocol") %>%
        plotly::layout(
          xaxis = list(title = "Days since trial start", zeroline = FALSE),
          yaxis = list(title = paste0("Doses to dispense ",
                                      "(", dfb_pp$dose_unit[1], ")"),
                       zeroline = FALSE),
          annotations = list(
            x = 0.5, y = 1,
            text = paste0("<b>", dfb_pp$kit_name[1], "</b>"),
            xanchor = "center", yanchor = "bottom",
            showarrow = FALSE, xref = 'paper', yref = 'paper'))
    }
  } else {
    for (j in 1:l) {
      dfa <- filter(dosing_pred_df, .data$kit == j &
                      .data$parameter == "observed data")
      dfb <- filter(dosing_pred_df, .data$kit == j &
                      .data$parameter == "model based prediction")
      dfb_pp <- filter(dosing_pred_df, .data$kit == j &
                         .data$parameter ==
                         "protocol based prediction")

      fig[[j]] <- plotly::plot_ly() %>%
        plotly::add_lines(
          data = dfa, x = ~date, y = ~n, line = list(shape ="hv", width = 2),
          name = "observed") %>%
        plotly::add_lines(
          data = dfb, x = ~date, y = ~n, line = list(width = 2),
          name = "median prediction model") %>%
        plotly::add_ribbons(
          data = dfb, x = ~date, ymin = ~lower, ymax = ~upper,
          fill = "tonexty", line = list(width = 0),
          name = "prediction interval model") %>%
        plotly::add_lines(
          data = dfb_pp, x = ~date, y = ~n, line = list(width = 2),
          name = "median prediction protocol") %>%
        plotly::add_ribbons(
          data = dfb_pp, x = ~date, ymin = ~lower, ymax = ~upper,
          fill = "tonexty", line = list(width = 0),
          name = "prediction interval protocol") %>%
        plotly::add_lines(
          x = rep(cutoffdt, 2),
          y = c(min(dfa$n), max(dfb$upper, dfb_pp$upper, na.rm = TRUE)),
          line = list(dash = "dash"), showlegend = FALSE,
          name = "cutoff") %>%
        plotly::layout(
          xaxis = list(title = "", zeroline = FALSE),
          yaxis = list(title = paste0("Doses to dispense ",
                                      "(", dfa$dose_unit[1], ")"),
                       zeroline = FALSE),
          annotations = list(
            x = 0.5, y = 1,
            text = paste0("<b>", dfa$kit_name[1], "</b>"),
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
    list(dosing_pred_df = dosing_pred_df,
         dosing_pred_plot = fig)
  } else {
    if (!pred_pp_only) {
      list(trialsdt = trialsdt, cutoffdt = cutoffdt,
           dosing_summary_t0 = observed$dosing_summary_t0,
           cum_dispense_plot = observed$cum_dispense_plot,
           bar_t0_plot = observed$bar_t0_plot,
           bar_ti_plot = observed$bar_ti_plot,
           bar_di_plot = observed$bar_di_plot,
           common_time_model = fit$common_time_model,
           k0_fit = fit$k0_fit, t0_fit = fit$t0_fit,
           t1_fit = fit$t1_fit, ki_fit = fit$ki_fit,
           ti_fit = fit$ti_fit, di_fit = fit$di_fit,
           dosing_subject = dosing_subject,
           dosing_pred_df = dosing_pred_df,
           dosing_pred_plot = fig)
    } else {
      list(trialsdt = trialsdt, cutoffdt = cutoffdt,
           dosing_summary_t0 = observed$dosing_summary_t0,
           cum_dispense_plot = observed$cum_dispense_plot,
           bar_t0_plot = observed$bar_t0_plot,
           bar_ti_plot = observed$bar_ti_plot,
           bar_di_plot = observed$bar_di_plot,
           dosing_subject = dosing_subject,
           dosing_pred_df = dosing_pred_df,
           dosing_pred_plot = fig)
    }
  }
}
