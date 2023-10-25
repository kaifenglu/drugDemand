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

# drugDemand 0.1.0

- initial release
