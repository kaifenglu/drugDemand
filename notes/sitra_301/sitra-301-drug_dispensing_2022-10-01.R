library(dplyr)
library(eventPred)
library(drugDemand)

set.seed(2000)

df <- readxl::read_excel(
  "notes/sitra_301/sitra-301-datacut-2022-10-01.xlsx") %>%
  dplyr::mutate(trialsdt = as.Date(trialsdt),
                randdt = as.Date(randdt),
                cutoffdt = as.Date(cutoffdt))
cutoffdt = df$cutoffdt[1]


visitview <- read.csv(
  "notes/sitra_301/BEIGE20210007_Blinded_VisitSummary_VisitView_22-Aug-2023 06_40_15.csv"
) %>%
  dplyr::rename(usubjid = Subject.Number,
                visit = Visit.Description) %>%
  dplyr::mutate(drug = ifelse(grepl('^1', Kit.Number) |
                                grepl('^2', Kit.Number), 2,
                              ifelse(grepl('^3', Kit.Number), 1,
                                     3))) %>%
  dplyr::mutate(drug_name = ifelse(drug == 1, "Tislelizumab",
                                   ifelse(drug == 2, "Sitravatinib",
                                          "Docetaxel"))) %>%
  dplyr::mutate(dose_unit = 'kit') %>%
  dplyr::filter(!grepl("^\\s*$", Date.Kits.Dispensed)) %>%
  dplyr::mutate(date = as.Date(Date.Kits.Dispensed, "%d-%b-%Y"),
                dispensed_quantity = as.numeric(Dispensed.Quantity)) %>%
  dplyr::mutate(kit_number = Kit.Number) %>%
  dplyr::select(usubjid, visit, date, drug, drug_name, dose_unit,
                kit_number, dispensed_quantity) %>%
  dplyr::filter(date <= cutoffdt)


dosing_schedule_df = dplyr::tibble(
  drug = c(1, 2, 3),
  target_days = c(21, 21, 21),
  target_dose = c(2, 2, 7),
  max_cycles = c(10000, 10000, 10000))

tictoc::tic("event prediction")

pred <- eventPred::getPrediction(
  df = df,
  to_predict = "enrollment and event",
  target_n = 385,
  target_d = 385,
  enroll_model = "b-spline",
  nknots = 2,
  lags = 150,
  event_model = "log-logistic",
  dropout_model = "none",
  pilevel = 0.95,
  nyears = 2,
  nreps = 200,
  showplot = FALSE,
  by_treatment = TRUE,
  alloc = c(2, 2))

tictoc::toc()


tictoc::tic("drug demand prediction")

drug_demand <- f_drug_demand(
  df = df,
  newEvents = pred$event_pred$newEvents,
  visitview = visitview,
  dosing_schedule_df = dosing_schedule_df,
  model_k0 = "zip",
  model_t0 = "log-logistic",
  model_t1 = "lad",
  model_ki = "zip",
  model_ti = "lad",
  model_di = "lme",
  pilevel = 0.95,
  nyears = 2,
  showplot = FALSE)

tictoc::toc()


writexl::write_xlsx(drug_demand$dosing_subject, paste0(
  "notes/sitra_301/sitra_301_dosing_subject_", cutoffdt, ".xlsx"))

writexl::write_xlsx(drug_demand$dosing_pred_df, paste0(
  "notes/sitra_301/sitra_301_dosing_pred_df_", cutoffdt, ".xlsx"))

writexl::write_xlsx(drug_demand$dosing_pred_pp, paste0(
  "notes/sitra_301/sitra_301_dosing_pred_pp_", cutoffdt, ".xlsx"))

