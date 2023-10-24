# drugDemand 0.1.1

- add a reference for parametric analysis of interval-censored survival data
- only keep one record per subject and drug dispensing day when using a common time model
- remove the vf_new parameter of the f_dosing_draw_1 function
- remove the nreps parameter from the f_drug_demand function
- add the f_dose_observed function and use it in f_drug_demand
- change the condition of Vi + Ti <= D(i) to Vi + Ti < D(i) in f_dose_ongoing_cpp and f_dose_new_cpp
- change "as.numeric(exp(attr(a$apVar, "Pars")))"" to "exp(as.numeric(attr(a$apVar, "Pars")))"" in the f_fit_di function to avoid the error for non-numeric argument to mathematical function
- simplify the condition for common_time_model to "length(unique(target_days)) == 1"

# drugDemand 0.1.0

- Initial release
