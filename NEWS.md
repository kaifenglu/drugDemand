# drugDemand 0.1.3

- remove the dependence on pscl and icenReg packages
- add stringr to imported packages
- add survival::survreg and stats::optim to imported functions
- use the Brent method for optimization to fit an interval-censored exponential distribution to the time from randomization to the first drug dispensing visit
- check the number of unique combinations of time and k1 in f_fit_ti
- use the VarCorr function from the nlme package to obtain the variance of random effects from the lme function call
- use data from all drugs to inform the common time model
- convert drug_name to a factor to ensure the correct order when creating a plot
- add the pred_pp_only parameter to f_drug_demand to make protocol-based predictions only
- ensure that df and visitview have all the required columns and none of the required columns have missing values if the data sets are provided in the f_drug_demand function call
- combine dosing_pred_df and dosing_pred_pp in the f_drug_demand output
- handle cases where all patients in a treatment arm discontinued before the cutoff
- replace round with formatC to retain the zeros after the decimal point
- rename fit_xx to xx_fit, where xx = k0, t0, t1, ki, ti, di to be consistent with enroll_fit and event_fit naming convention
- add drug_name as the sub plot title in f_dispensing_models.R

# drugDemand 0.1.2

- add least absolute deviations (LAD) regression as an option for modelling the gap time between drug dispensing visits and rename the original linear model as least squares (LS)
- use more descriptive names for drug dispensing models

# drugDemand 0.1.1

- add a reference for parametric analysis of interval-censored survival data
- only keep one record per subject and drug dispensing day when using a common time model
- rename f_dosing_draw, f_dosing_draw_1, and f_dosing_draw_t_1 to f_dose_draw, f_dose_draw_1, and f_dose_draw_t_1, respectively
- remove the vf_new parameter of the f_dose_draw_1 function
- remove the nreps parameter from the f_drug_demand function
- add the f_dose_observed function and incorporate it in the f_drug_demand function
- modify the condition "Vi + Ti <= D(i)" to "Vi + Ti < D(i)" in the f_dose_ongoing_cpp and f_dose_new_cpp functions
- change "as.numeric(exp(attr(a$apVar, "Pars")))" to "exp(as.numeric(attr(a$apVar, "Pars")))" in the f_fit_di function to avoid the error for non-numeric argument to mathematical function
- simplify the condition for common_time_model to "length(unique(target_days)) == 1"
- add dosing_summary_t0 to the output of the f_drug_demand function
- replace mutate and slice(n()) with summarise in the f_dose_observed and f_dose_draw functions to improve efficiency
- plot gap time t0 instead of t0 + 1 in the f_fit_t0 function
- plot the rounded value of di based on probability calculations alongside the observed value of di in the f_fit_di function
- replace the residual plot with the fitted gap time bar chart in the f_fit_ti function
- redefine row_id for vf1 if common_time_model is true in the f_dispensing_models function
- use df and visitview to derive treatment_by_drug_df for real-time drug demand forecasting
- update the examples of the f_fit_t0, f_fit_ki, and f_fit_ti functions
- add trialsdt and cutoffdt to the output of the f_dose_observed and f_drug_demand functions
- move the arrange operation of dosing_subject_newi out of the f_dose_draw_1 function into the f_dose_draw function to improve efficiency
- combine the two summarise operation of dosing_summary_newi in the f_dose_draw_1 function to improve efficiency
- replace the zero-inflated negative binomial distribution with the negative binomial distribution in the f_fit_ki function to avoid convergence issues
- print cum_dispense_plot if showplot is TRUE in the f_dose_observed function
- add colors = "Set2" to cum_dispense_plot, bar_t0_plot, bar_ti_plot, and bar_di_plot in the f_dose_observed function
- remove the custom legend of cum_dispense_plot in the f_dose_observed and f_drug_demand functions
- add parameter l to the f_dose_draw_1 function to improve efficiency
- add structure and more details to the function parameters and output descriptions
- combine the dosing_subject_t and dosing_summary_t steps in the f_dose_observed function to improve efficiency
- replace target_days with dosing_schedule_df in the argument for the f_dispensing_models function
- drop the creation of the status variable and use table instead of survfit for observed data summary in the f_fit_t0 function
- add the handling of one observation case in the f_fit_ti function
- remove the creation of the unames1 variable in the f_dose_draw function

# drugDemand 0.1.0

- initial release
