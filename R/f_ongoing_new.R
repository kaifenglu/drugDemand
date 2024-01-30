#' @title Observed Dosing for Ongoing and New Subjects
#' @description Prepares the dosing data sets to impute for ongoing and
#' new subjects.
#'
#' @param newEvents A data frame containing the imputed event data
#'   for both ongoing and new patients, typically obtained from
#'   the output of the \code{getPrediction} function of the
#'   \code{eventPred} package. It contains the following variables:
#'   \code{draw}, \code{usubjid}, \code{arrivalTime}, \code{treatment},
#'   \code{treatment_description}, \code{time}, \code{event},
#'   \code{dropout}, and \code{totalTime}.
#' @param drug_description_df The drug description data frame
#'   including \code{drug}, \code{drug_name}, \code{dose_strength},
#'   \code{kit}, \code{kit_name}, and \code{dose_unit}.
#'   It must be specified at the design stage. It will be replaced with
#'   the observed information at the analysis stage.
#' @param treatment_by_drug_df The data frame indicating the treatments
#'   associated with each drug, including the following variables:
#'   \code{treatment} and \code{drug}.
#' @param vf A data frame for subject-level drug dispensing data,
#'   including the following variables: \code{drug}, \code{drug_name},
#'   \code{dose_strength}, \code{kit}, \code{kit_name},
#'   \code{dose_unit}, \code{usubjid}, \code{treatment},
#'   \code{treatment_description}, \code{arrivalTime}, \code{time},
#'   \code{event}, \code{dropout}, \code{day}, \code{dose},
#'   \code{cum_dose}, and \code{row_id}.
#'
#' @return A list with the following components:
#'
#' * \code{vf_ongoing}: The observed drug dispensing data for ongoing
#'   patients with drug dispensing records. It includes the following
#'   variables: \code{draw}, \code{kit}, \code{kit_name}, \code{dose_unit},
#'   \code{usubjid}, \code{day}, \code{dose}, \code{arrivalTime},
#'   \code{treatment}, \code{treatment_description},
#'   \code{time}, and \code{totalTime}.
#'
#' * \code{vf_new}: A data frame for the randomization date for new patients
#'   and ongoing patients with no drug dispensing records.
#'   It includes the following variables:
#'   \code{draw}, \code{kit}, \code{kit_name}, \code{dose_unit},
#'   \code{usubjid}, \code{arrivalTime}, \code{treatment},
#'   \code{treatment_description}, \code{time}, and \code{totalTime}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' \donttest{
#' set.seed(2000)
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
#' observed <- f_dose_observed(df = df2, visitview = visitview2)
#'
#' vf_ongoing_new <- f_ongoing_new(
#'   newEvents = pred$event_pred$newEvents,
#'   drug_description_df = observed$drug_description_df,
#'   treatment_by_drug_df = observed$treatment_by_drug_df,
#'   vf = observed$vf)
#'
#' head(vf_ongoing_new$vf_ongoing)
#' head(vf_ongoing_new$vf_new)
#' }
#'
#' @export
f_ongoing_new <- function(
    newEvents, drug_description_df, treatment_by_drug_df, vf) {
  nreps = length(unique(newEvents$draw))

  if (!is.null(vf)) {
    ### dosing data for ongoing patients ###
    vf1 <- vf %>%
      filter(.data$event == 0) %>%
      ungroup() %>%
      select(c("kit", "kit_name", "dose_unit",
               "usubjid", "day", "dose"))

    # ongoing subjects with dosing records
    unames <- unique(vf1$usubjid)

    # replicate nreps times
    vf1_rep = tibble(draw = 1:nreps) %>% cross_join(vf1)

    df1_ongoing <- newEvents %>%
      filter(.data$usubjid %in% unames) %>%
      select(-c("event", "dropout"))

    vf_ongoing <- vf1_rep %>%
      inner_join(df1_ongoing, by = c("draw", "usubjid"))

    ### new patients and ongoing patients with no dosing records ###
    df_new <- newEvents %>% filter(!(.data$usubjid %in% unames))

    J = length(unique(drug_description_df$drug))
    if (nrow(df_new) > 0) {
      vf_new1 <- purrr::map_dfr(
        1:J, function(j) {
          df_new %>%
            inner_join(treatment_by_drug_df %>% filter(.data$drug == j),
                       by = "treatment")
        }) %>% select(-c("event", "dropout"))

      # draw kit probabilities from Dirichlet distribution
      vf2 <- vf %>%
        slice(1) %>%
        group_by(.data$drug, .data$drug_name, .data$dose_strength,
                 .data$kit, .data$kit_name, .data$dose_unit) %>%
        summarise(n = n(), .groups = "drop_last")

      p_kit <- purrr::map(1:J, function(j) {
        rdirichlet(nreps, vf2$n[vf2$drug == j] + 1)
      })

      # number of kits per drug
      nkits <- as.numeric(table(drug_description_df$drug))
      offset <- c(0, cumsum(nkits))[1:J]

      # loop over drugs
      vf_new <- purrr::map_dfr(1:J, function(j) {
        if (nkits[j] == 1) {
          vf_new1 %>% filter(.data$drug == j) %>%
            mutate(kit = offset[j] + 1) %>%
            inner_join(drug_description_df, by = c("drug", "kit"))
        } else {
          purrr::map_dfr(1:nreps, function(i) {
            df_new1 <- vf_new1 %>% filter(.data$draw == i & .data$drug == j)
            kit_ind <- t(rmultinom(nrow(df_new1), 1, p_kit[[j]][i,]))
            df_new1$kit <- as.numeric(kit_ind %*% seq(1,nkits[j])) + offset[j]
            df_new1 %>% inner_join(drug_description_df, by = c("drug", "kit"))
          })
        }
      }) %>% select(-c("drug", "drug_name", "dose_strength"))
    } else {
      vf_new <- NULL
    }
  } else {
    vf_ongoing <- NULL

    # drug and kit information for new subjects
    J = length(unique(drug_description_df$drug))
    vf_new1 <- purrr::map_dfr(
      1:J, function(j) {
        newEvents %>%
          inner_join(treatment_by_drug_df %>% filter(.data$drug == j),
                     by = "treatment")
      }) %>% select(-c("event", "dropout"))

    # add kit information for new subjects
    p_kit <- purrr::map(1:J, function(j) {
      drug_description_df$p_kit[drug_description_df$drug == j]
    })

    # number of kits per drug
    nkits <- as.numeric(table(drug_description_df$drug))
    offset <- c(0, cumsum(nkits))[1:J]

    # loop over drugs
    vf_new <- purrr::map_dfr(1:J, function(j) {
      if (nkits[j] == 1) {
        vf_new1 %>% filter(.data$drug == j) %>%
          mutate(kit = offset[j] + 1) %>%
          inner_join(drug_description_df, by = c("drug", "kit"))
      } else {
        purrr::map_dfr(1:nreps, function(i) {
          df_new1 <- vf_new1 %>% filter(.data$draw == i & .data$drug == j)
          kit_ind <- t(rmultinom(nrow(df_new1), 1, p_kit[[j]][i,]))
          df_new1$kit <- as.numeric(kit_ind %*% seq(1,nkits[j])) + offset[j]
          df_new1 %>% inner_join(drug_description_df, by = c("drug", "kit"))
        })
      }
    }) %>% select(-c("drug", "drug_name", "dose_strength"))
  }

  list(vf_ongoing = vf_ongoing, vf_new = vf_new)
}
