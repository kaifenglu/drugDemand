#' The drug description data frame.
#'
#' A data frame with 4 rows and 3 columns:
#' \describe{
#'   \item{\code{drug}}{The numeric code of the drug.}
#'   \item{\code{drug_name}}{The name of the drug.}
#'   \item{\code{dose_unit}}{The dose unit for drug dispensing.}
#' }
#' For drug demand forecasting, the default dose unit is "kit" for all drugs.
"drug_description_df"


#' The indicator matrix of treatment by drug combinations.
#'
#' A matrix with dimensions k x l, where k equals 3 representing the
#' number of treatment groups, and l equals 4 representing the number
#' of drugs. In this matrix, a value of 1 signifies the presence of
#' the drug within a treatment group, while a value of 0 indicates
#' the absence of the drug in that particular treatment group.
"treatment_by_drug"


#' The dosing schedule data frame.
#'
#' A data frame with 4 rows and 4 columns:
#' \describe{
#'   \item{\code{drug}}{The numeric code of the drug.}
#'   \item{\code{target_days}}{The target number of days per treatment cycle.}
#'   \item{\code{target_kits}}{The target number of kits per treatment cycle.}
#'   \item{\code{max_cycles}}{The maximum number of treatment cycles.}
#' }
"dosing_schedule_df"


#' The subject-level enrollment and event data before enrollment completion.
#'
#' A data frame with 175 rows and 9 columns:
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
#' A data frame with 250 rows and 9 columns:
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
#' A data frame with 2290 rows and 8 columns:
#' \describe{
#'   \item{\code{usubjid}}{The unique subject ID.}
#'   \item{\code{visit}}{The drug dispensing visit, e.g. "Cycle 1 Day 1".}
#'   \item{\code{date}}{The date of the drug dispensing visit.}
#'   \item{\code{drug}}{The numeric code of the drug.}
#'   \item{\code{drug_name}}{The name of the drug.}
#'   \item{\code{dose_unit}}{The dose unit for drug dispensing.}
#'   \item{\code{kit_number}}{The kit number for drug dispensing.}
#'   \item{\code{dispensed_quantity}}{The number of kits dispensed
#'   at the visit.}
#' }
"visitview1"


#' The observed subject drug dispensing data after enrollment completion.
#'
#' A data frame with 5006 rows and 8 columns:
#' \describe{
#'   \item{\code{usubjid}}{The unique subject ID.}
#'   \item{\code{visit}}{The drug dispensing visit, e.g. "Cycle 1 Day 1".}
#'   \item{\code{date}}{The date of the drug dispensing visit.}
#'   \item{\code{drug}}{The numeric code of the drug.}
#'   \item{\code{drug_name}}{The name of the drug.}
#'   \item{\code{dose_unit}}{The dose unit for drug dispensing.}
#'   \item{\code{kit_number}}{The kit number for drug dispensing.}
#'   \item{\code{dispensed_quantity}}{The number of kits dispensed
#'   at the visit.}
#' }
"visitview2"

