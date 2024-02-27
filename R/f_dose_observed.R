f_bar_chart <- function(x) {
  u = table(x)
  count = as.numeric(u)
  y.obs = as.numeric(names(u))
  ymax = max(y.obs)
  y = 0:ymax
  n = rep(0, ymax+1)
  n[y.obs+1] = count
  tibble(y = y, n = n)
}


#' @title Observed Drug Dispensing Data Summary
#' @description Provides an overview of the observed drug dispensing data,
#' including the summary of cumulative dispensed doses, bar chart
#' of the gap time between randomization and the first drug dispensing
#' visit, the gap time between two consecutive drug dispensing visits,
#' and the dispensed doses at drug dispensing visits by drug.
#'
#' @param df A data frame for subject-level enrollment and event data,
#'   including the following variables:
#'   \code{trialsdt}, \code{usubjid}, \code{randdt},
#'   \code{treatment}, \code{treatment_description},
#'   \code{time}, \code{event}, \code{dropout}, and \code{cutoffdt}.
#' @param visitview A data frame containing the observed drug dispensing
#'   data, including the following variables:
#'   \code{usubjid}, \code{visit}, \code{date}, \code{drug},
#'   \code{drug_name}, \code{kit}, \code{kit_name},
#'   \code{kit_number}, and \code{dispensed_quantity}.
#' @param showplot A Boolean variable that controls whether or not to
#'   show the drug dispensing model fit and drug demand prediction
#'   plots. It defaults to \code{TRUE}.
#'
#' @return A list with the following components:
#'
#' * \code{trialsdt}: The trial start date.
#'
#' * \code{cutoffdt}: The cutoff date.
#'
#' * \code{vf}: A data frame for subject-level drug dispensing data,
#'   including the following variables:
#'   \code{drug}, \code{drug_name}, \code{kit}, \code{kit_name},
#'   \code{usubjid}, \code{treatment}, \code{treatment_description},
#'   \code{arrivalTime}, \code{time}, \code{event}, \code{dropout},
#'   \code{day}, \code{dose}, \code{cum_dose}, and \code{row_id}.
#'
#' * \code{treatment_by_drug_df}: A data frame indicating the treatments
#'   associated with each drug, including the following variables:
#'   \code{treatment} and \code{drug}.
#'
#' * \code{kit_description_df}: A data frame indicating the
#'   drug and kit descriptions, including the following variables:
#'   \code{drug}, \code{drug_name}, \code{kit}, and \code{kit_name}.
#'
#' * \code{dosing_summary_t}: A data frame for the cumulative doses
#'   dispensed by each observed time point. It contains the following
#'   variables:
#'   \code{kit}, \code{kit_name}, \code{t}, \code{n},
#'   \code{lower}, \code{upper}, \code{mean}, and \code{var}, where
#'   \code{lower} and \code{upper} have missing values,
#'   \code{mean = n}, and \code{var = 0}.
#'
#' * \code{dosing_summary_t0}: A data frame for the cumulative doses
#'   dispensed before the cutoff date. It contains the following
#'   variables:
#'   \code{kit}, \code{kit_name}, and \code{cum_dose_t0}.
#'
#' * \code{cum_dispense_plot}: The step plot for the cumulative doses
#'   dispensed for each kit type.
#'
#' * \code{bar_t0_plot}: The bar chart for the gap time between
#'   randomization and the first drug dispensing visit.
#'
#' * \code{bar_ti_plot}: The bar chart for the gap time between two
#'   consecutive drug dispensing visits.
#'
#' * \code{bar_di_plot}: The bar chart for the dispensed doses at drug
#'   dispensing visits.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' observed <- f_dose_observed(df = df2, visitview = visitview2)
#'
#' @export
f_dose_observed <- function(df = NULL, visitview = NULL, showplot = TRUE) {
  trialsdt = df$trialsdt[1]
  cutoffdt = df$cutoffdt[1]

  df <- df %>%
    mutate(arrivalTime = as.numeric(.data$randdt - .data$trialsdt + 1))

  # set up drug/kit/subject/day drug dispensing data
  vf <- visitview %>%
    inner_join(df, by = "usubjid") %>%
    mutate(day = as.numeric(.data$date - .data$randdt + 1)) %>%
    group_by(.data$drug, .data$drug_name, .data$kit, .data$kit_name,
             .data$usubjid, .data$treatment, .data$treatment_description,
             .data$arrivalTime, .data$time, .data$event, .data$dropout,
             .data$day) %>%
    summarise(dose = sum(.data$dispensed_quantity),
              .groups = "drop_last") %>%
    mutate(cum_dose = cumsum(.data$dose)) %>%
    group_by(.data$drug, .data$drug_name, .data$kit, .data$kit_name,
             .data$usubjid) %>%
    mutate(row_id = row_number()) %>%
    select(c("drug", "drug_name", "kit", "kit_name",
             "usubjid", "treatment", "treatment_description",
             "arrivalTime", "time", "event", "dropout",
             "day", "dose", "cum_dose", "row_id"))

  treatment_by_drug_df <- vf %>%
    group_by(.data$treatment, .data$drug) %>%
    slice(1) %>%
    select(c("treatment", "drug"))

  # extract drug description from observed data
  kit_description_df <- vf %>%
    group_by(.data$drug, .data$drug_name, .data$kit, .data$kit_name) %>%
    slice(1) %>%
    select(c("drug", "drug_name", "kit", "kit_name"))

  # obtain the observed time points relative to trial start
  t_df <- vf %>% mutate(t1 = .data$arrivalTime + .data$day - 1)
  t_obs <- sort(unique(t_df$t1))

  # tally the doses across patients
  dosing_summary_t <- tibble(t = t_obs) %>%
    cross_join(vf) %>%
    filter(.data$arrivalTime + .data$day - 1 <= .data$t) %>%
    group_by(.data$kit, .data$kit_name, .data$t) %>%
    summarise(n = sum(.data$dose), .groups = "drop_last") %>%
    mutate(lower = NA, upper = NA, mean = .data$n, var = 0)

  # obtain the cumulative doses up to cutoff
  dosing_summary_t0 <- dosing_summary_t %>%
    group_by(.data$kit, .data$kit_name) %>%
    slice(n()) %>%
    rename(cum_dose_t0 = .data$n) %>%
    select(c("kit", "kit_name", "cum_dose_t0"))

  # set up treatment by drug combinations
  l = nrow(kit_description_df)
  kit_name = kit_description_df$kit_name

  # construct the plot of cumulative drug dispensing data
  # initialize the dosing plot data set
  df0 <- tibble(kit = 1:l, kit_name = kit_name, t = 1, n = 0,
                lower = NA, upper = NA, mean = 0, var = 0)

  ad <- df0 %>%
    bind_rows(dosing_summary_t) %>%
    mutate(date = as.Date(.data$t - 1, origin = trialsdt))

  # convert kit_name to a factor to ensure the correct order
  ad$kit_name <- factor(ad$kit_name, levels = kit_name)

  cum_dispense_plot <- plotly::plot_ly(
    ad, x = ~date, y = ~n, color = ~kit_name, colors = "Set2") %>%
    plotly::add_lines(line = list(shape = "hv", width = 2)) %>%
    plotly::layout(
      xaxis = list(title = ""),
      yaxis = list(title = "Doses dispensed", zeroline = FALSE))

  bar_t0_df <- bar_ti_df <- bar_di_df <- tibble()
  for (h in 1:l) {
    # observed dosing data for the kit under consideration
    vf1 <- vf %>% filter(.data$kit == h)

    # time from randomization to the first drug dispensing visit
    df_t0 <- vf1 %>%
      filter(.data$row_id == 1) %>%
      mutate(time = .data$day - 1)

    # gap time and number of skipped visits between drug dispensing visits
    df_ti <- vf1 %>%
      mutate(time = lead(.data$day) - .data$day) %>%
      filter(.data$row_id < n())

    # construct the data sets for the bar charts
    bar_t0_df <- bar_t0_df %>%
      bind_rows(f_bar_chart(df_t0$time) %>% mutate(kit = h))

    bar_ti_df <- bar_ti_df %>%
      bind_rows(f_bar_chart(df_ti$time) %>% mutate(kit = h))

    bar_di_df <- bar_di_df %>%
      bind_rows(f_bar_chart(vf1$dose) %>% mutate(kit = h))
  }

  # construct the bar chart for t0
  bar_t0_df <- bar_t0_df %>%
    left_join(kit_description_df, by = "kit")

  # convert kit_name to a factor to ensure the correct order
  bar_t0_df$kit_name <- factor(bar_t0_df$kit_name, levels = kit_name)

  bar_t0_plot <- plotly::plot_ly(
    bar_t0_df, x = ~y, y = ~n, type = 'bar',
    color = ~kit_name, colors = "Set2") %>%
    plotly::layout(
      xaxis = list(title = paste('Days between randomization and',
                                 'the first drug dispensing visit')),
      yaxis = list(title = 'Count'),
      barmode = 'group')

  # construct the bar chart for ti
  bar_ti_df <- bar_ti_df %>%
    left_join(kit_description_df, by = "kit")

  # convert kit_name to a factor to ensure the correct order
  bar_ti_df$kit_name <- factor(bar_ti_df$kit_name, levels = kit_name)

  bar_ti_plot <- plotly::plot_ly(
    bar_ti_df, x = ~y, y = ~n, type = 'bar',
    color = ~kit_name, colors = "Set2") %>%
    plotly::layout(
      xaxis = list(title = paste('Days between consecutive',
                                 'drug dispensing visits')),
      yaxis = list(title = 'Count'),
      barmode = 'group')

  # construct the bar chart for di
  bar_di_df <- bar_di_df %>%
    left_join(kit_description_df, by = "kit")

  # convert kit_name to a factor to ensure the correct order
  bar_di_df$kit_name <- factor(bar_di_df$kit_name, levels = kit_name)

  bar_di_plot <- plotly::plot_ly(
    bar_di_df, x = ~y, y = ~n, type = 'bar',
    color = ~kit_name, colors = "Set2") %>%
    plotly::layout(
      xaxis = list(title = paste('Doses dispensed at',
                                 'drug dispensing visits')),
      yaxis = list(title = 'Count'),
      barmode = 'group')

  if (showplot) {
    print(cum_dispense_plot)
    print(bar_t0_plot)
    print(bar_ti_plot)
    print(bar_di_plot)
  }

  list(trialsdt = trialsdt,
       cutoffdt = cutoffdt,
       vf = vf,
       treatment_by_drug_df = treatment_by_drug_df,
       kit_description_df = kit_description_df,
       dosing_summary_t = dosing_summary_t,
       dosing_summary_t0 = dosing_summary_t0,
       cum_dispense_plot = cum_dispense_plot,
       bar_t0_plot = bar_t0_plot,
       bar_ti_plot = bar_ti_plot,
       bar_di_plot = bar_di_plot)
}

