#' The kit description data frame.
#'
#' A data frame with the following columns:
#' \describe{
#'   \item{\code{drug}}{The numeric code of the drug.}
#'   \item{\code{drug_name}}{The name of the drug.}
#'   \item{\code{kit}}{The numeric code of the kit type.}
#'   \item{\code{kit_name}}{The name of the kit type.}
#'   \item{\code{p_kit}}{The prior probability of different kit types
#'   within a drug.}
#' }
"kit_description_df"


#' The data frame indicating the treatments associated with each drug.
#'
#' A data frame with the following columns:
#' \describe{
#'   \item{\code{treatment}}{The numeric code of the treatment group.}
#'   \item{\code{drug}}{The numeric code of the drug.}
#'  }
"treatment_by_drug_df"


#' The dosing schedule data frame.
#'
#' A data frame with the following columns:
#' \describe{
#'   \item{\code{kit}}{The numeric code of the kit type.}
#'   \item{\code{target_days}}{Number of days per treatment cycle.}
#'   \item{\code{target_dose}}{Dose per treatment cycle.}
#'   \item{\code{max_cycles}}{Maximum number of treatment cycles.}
#' }
"dosing_schedule_df"


#' The subject-level enrollment and event data before enrollment completion.
#'
#' A data frame with the following columns:
#' \describe{
#'   \item{\code{trialsdt}}{The trial start date.}
#'   \item{\code{usubjid}}{The unique subject ID.}
#'   \item{\code{randdt}}{The randomization date for each subject.}
#'   \item{\code{treatment}}{The treatment group.}
#'   \item{\code{treatment_description}}{Description of the treatment group.}
#'   \item{\code{time}}{The number of days elapsed since randomization.}
#'   \item{\code{event}}{The event indicator, with a value of 1 indicating
#'   the occurrence of an event, and 0 indicating no event.}
#'   \item{\code{dropout}}{The dropout indicator, where 1 corresponds to
#'   a dropout and 0 implies no dropout.}
#'   \item{\code{cutoffdt}}{The cutoff date.}
#' For drug demand forecasting, the event of interest is treatment
#' discontinuation. The dropout variable is set to 0 for all patients in
#' this context.
#' }
"df1"


#' The subject-level enrollment and event data after enrollment completion.
#'
#' A data frame with the following columns:
#' \describe{
#'   \item{\code{trialsdt}}{The trial start date.}
#'   \item{\code{usubjid}}{The unique subject ID.}
#'   \item{\code{randdt}}{The randomization date for each subject.}
#'   \item{\code{treatment}}{The treatment group.}
#'   \item{\code{treatment_description}}{Description of the treatment group.}
#'   \item{\code{time}}{The number of days elapsed since randomization.}
#'   \item{\code{event}}{The event indicator, with a value of 1 indicating
#'   the occurrence of an event, and 0 indicating no event.}
#'   \item{\code{dropout}}{The dropout indicator, where 1 corresponds to
#'   a dropout and 0 implies no dropout.}
#'   \item{\code{cutoffdt}}{The cutoff date.}
#' For drug demand forecasting, the event of interest is treatment
#' discontinuation. The dropout variable is set to 0 for all patients in
#' this context.
#' }
"df2"


#' The observed subject drug dispensing data before enrollment completion.
#'
#' A data frame with the following columns:
#' \describe{
#'   \item{\code{usubjid}}{The unique subject ID.}
#'   \item{\code{visit}}{The drug dispensing visit, e.g., "Cycle 1 Day 1".}
#'   \item{\code{date}}{The date of the drug dispensing visit.}
#'   \item{\code{drug}}{The numeric code of the drug.}
#'   \item{\code{drug_name}}{The name of the drug.}
#'   \item{\code{kit}}{The numeric code of the kit type.}
#'   \item{\code{kit_name}}{The name of the kit type.}
#'   \item{\code{kit_number}}{The kit number for drug dispensing.}
#'   \item{\code{dispensed_quantity}}{The number of kits dispensed
#'   at the visit.}
#' }
"visitview1"


#' The observed subject drug dispensing data after enrollment completion.
#'
#' A data frame with the following columns:
#' \describe{
#'   \item{\code{usubjid}}{The unique subject ID.}
#'   \item{\code{visit}}{The drug dispensing visit, e.g., "Cycle 1 Day 1".}
#'   \item{\code{date}}{The date of the drug dispensing visit.}
#'   \item{\code{drug}}{The numeric code of the drug.}
#'   \item{\code{drug_name}}{The name of the drug.}
#'   \item{\code{kit}}{The numeric code of the kit type.}
#'   \item{\code{kit_name}}{The name of the kit type.}
#'   \item{\code{kit_number}}{The kit number for drug dispensing.}
#'   \item{\code{dispensed_quantity}}{The number of kits dispensed
#'   at the visit.}
#' }
"visitview2"

