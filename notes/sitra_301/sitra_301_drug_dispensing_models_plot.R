library(dplyr)
library(plotly)

tictoc::tic("plot")

# cutoffdt = as.Date("2023-03-01")
cutoffdt = as.Date("2022-10-01")

df <- readxl::read_excel(
  "notes/sitra_301/Sitra-301 upload R-Shiny - CN_22-Aug-2023 - no optional chemo.xlsx"
) %>%
  dplyr::mutate(trialsdt = as.Date(trialsdt),
                randdt = as.Date(randdt),
                cutoffdt = as.Date(cutoffdt))

trialsdt = df$trialsdt[1]
t0 = as.numeric(cutoffdt - trialsdt + 1)
t_half = as.numeric(df$cutoffdt[1] - trialsdt + 1)

drug_name = c("Tislelizumab", "Sitravatinib", "Docetaxel")
dose_unit = c("kit", "kit", "kit")
l = length(drug_name)

# visit view data supplemented with subject view data
visitview <- read.csv(
  "notes/sitra_301/BEIGE20210007_Blinded_VisitSummary_VisitView_22-Aug-2023 06_40_15.csv"
) %>%
  dplyr::rename(usubjid = Subject.Number) %>%
  dplyr::mutate(drug = ifelse(grepl('^1', Kit.Number) |
                                grepl('^2', Kit.Number), 2,
                              ifelse(grepl('^3', Kit.Number), 1,
                                     3))) %>%
  dplyr::mutate(drug_name = ifelse(drug == 1, "Tislelizumab",
                                   ifelse(drug == 2, "Sitravatinib",
                                          "Docetaxel"))) %>%
  dplyr::mutate(dose_unit = 'kit') %>%
  dplyr::inner_join(df %>% dplyr::mutate(
    adt = as.Date(time - 1, origin = randdt)),
    by = "usubjid") %>%
  dplyr::mutate(randdt = as.Date(randdt))


# calculate total number of doses per drug/subject/dispensing_date
vf <- visitview %>%
  dplyr::filter(!grepl("^\\s*$", Date.Kits.Dispensed)) %>%
  dplyr::mutate(date = as.Date(Date.Kits.Dispensed, "%d-%b-%Y")) %>%
  dplyr::select(usubjid, drug, drug_name, dose_unit, randdt, adt,
                event, dropout, date, Dispensed.Quantity) %>%
  dplyr::group_by(drug, drug_name, dose_unit, usubjid, randdt, adt,
                  event, dropout, date) %>%
  dplyr::summarise(dose = sum(Dispensed.Quantity),
                   .groups = "drop_last") %>%
  dplyr::mutate(cum_dose = cumsum(dose)) %>%
  dplyr::group_by(drug, drug_name, dose_unit, usubjid) %>%
  dplyr::mutate(row_id = dplyr::row_number())


# prediction based on protocol specified visit and drug dispensing schedules
dosing_overview_pp <- readxl::read_excel(
  paste0("notes/sitra_301/sitra_301_dosing_pred_pp_", cutoffdt, ".xlsx")) %>%
  dplyr::filter(t <= t_half + 30)

# prediction based on modeling of observed drug dispensing data
dosing_overview_mod <- readxl::read_excel(
  paste0("notes/sitra_301/sitra_301_dosing_pred_df_", cutoffdt, ".xlsx")) %>%
  dplyr::filter(t <= t_half + 30)


### plot
df0 <- dplyr::tibble(drug = 1:l,
                     drug_name = drug_name,
                     dose_unit = dose_unit,
                     t = 1, n = 0, lower = NA, upper = NA,
                     mean = 0, var = 0)


t_obs <- unique((vf %>% dplyr::mutate(
  t1 = as.numeric(date - trialsdt + 1)))$t1)

dosing_subject_t <- dplyr::tibble(t = t_obs) %>%
  dplyr::cross_join(vf) %>%
  dplyr::filter(date <= as.Date(t-1, origin = trialsdt)) %>%
  dplyr::group_by(drug, drug_name, dose_unit, t, usubjid) %>%
  dplyr::summarise(cum_dose = sum(dose), .groups = "drop_last")

dosing_summary_t <- dosing_subject_t %>%
  dplyr::group_by(drug, drug_name, dose_unit, t) %>%
  dplyr::summarise(n = sum(cum_dose), .groups = "drop_last") %>%
  dplyr::mutate(lower = NA, upper = NA, mean = n, var = 0)

dosing_pred_df <- df0 %>%
  dplyr::bind_rows(dosing_summary_t) %>%
  dplyr::mutate(date = as.Date(t-1, origin = trialsdt)) %>%
  dplyr::bind_rows(dosing_overview_mod) %>%
  dplyr::arrange(drug, drug_name, dose_unit, t)


g_dosing <- list()
for (j in 1:l) {
  dfa <- dplyr::filter(dosing_pred_df, drug == j & is.na(lower))
  dfb <- dplyr::filter(dosing_pred_df, drug == j & !is.na(lower))

  dfb_pp <- dplyr::filter(dosing_overview_pp, drug == j & !is.na(lower))

  g_dosing[[j]] <- plotly::plot_ly() %>%
    plotly::add_ribbons(
      data = dfb, x = ~date, ymin = ~lower, ymax = ~upper,
      fill = "tonexty", line = list(width=0),
      name = "prediction interval") %>%
    plotly::add_lines(
      data = dfb, x = ~date, y = ~n, name = "median prediction",
      line = list(width=2)) %>%
    plotly::add_lines(
      data = dfa, x = ~date, y = ~n, name = "observed",
      line = list(shape = "hv", width = 2)) %>%
    plotly::add_ribbons(
      data = dfb_pp, x = ~date, ymin = ~lower, ymax = ~upper,
      fill = "tonexty", line = list(width = 0),
      name = "prediction interval pp") %>%
    plotly::add_lines(
      data = dfb_pp, x = ~date, y = ~n, name = "median prediction pp",
      line = list(width = 2)) %>%
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


  if (j == 1) {
    g_dosing[[j]] <- g_dosing[[j]] %>%
      plotly::layout(
        annotations = list(
          x = cutoffdt, y = 0, text = 'cutoff',
          xanchor = "left", yanchor = "bottom",
          font = list(size = 12), showarrow = FALSE))

  }
}

g_dosing

tictoc::toc()
