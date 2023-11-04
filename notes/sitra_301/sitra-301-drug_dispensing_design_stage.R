library(dplyr)
library(eventPred)
library(drugDemand)

set.seed(2000)

drug_description_df = dplyr::tibble(
  drug = c(1, 2, 3),
  drug_name = c("Tislelizumab", "Sitravatinib", "Docetaxel"),
  dose_unit = c("kit", "kit", "kit"))

treatment_by_drug = matrix(c(1, 1, 0, 0, 0, 1),
                           nrow = 2, ncol = 3, byrow = TRUE)

dosing_schedule_df = dplyr::tibble(
  drug = c(1, 2, 3),
  target_days = c(21, 21, 21),
  target_dose = c(2, 2, 7),
  max_cycles = c(10000, 10000, 10000))

tictoc::tic("event prediction")

pred <- eventPred::getPrediction(
  df = NULL,
  to_predict = "enrollment and event",
  target_n = 385,
  target_d = 385,
  enroll_prior = list(
    model = "poisson",
    theta = -0.77,
    vtheta = 0.0030),
  event_prior = list(
    list(model = "log-logistic",
         theta = c(5.0, -0.54),
         vtheta = matrix(c(0.00658, 0.00037, 0.00037, 0.00515), 2, 2)),
    list(model = "log-logistic",
         theta = c(4.3, -0.35),
         vtheta = matrix(c(0.00860, 0.00011, 0.00011, 0.00513), 2, 2))),
  dropout_prior = NULL,
  pilevel = 0.95,
  nyears = 3,
  nreps = 200,
  showplot = FALSE,
  by_treatment = TRUE,
  ngroups = 2,
  alloc = c(2, 2),
  treatment_label = c("Tislelizumab + Sitravatinib",
                      "Docetaxel monotherapy"))

tictoc::toc()


tictoc::tic("drug demand prediction")

drug_demand <- f_drug_demand(
  df = NULL,
  newEvents = pred$event_pred$newEvents,
  drug_description_df = drug_description_df,
  treatment_by_drug = treatment_by_drug,
  dosing_schedule_df = dosing_schedule_df,
  pilevel = 0.95,
  nyears = 3,
  showplot = FALSE)

tictoc::toc()


writexl::write_xlsx(drug_demand$dosing_pred_pp, paste0(
  "notes/sitra_301/sitra_301_dosing_pred_design_stage.xlsx"))

