library(dplyr)
library(purrr)

set.seed(123)

n = 250
trialsdt = as.Date("2021-07-28")
cutoffdt = as.Date("2023-08-24")

# piecewise accrual
accrualTime = c(0, 240)
accrualIntensity = c(0.48, 0.31)

# generate arrival times from piecewise poisson process
f_enrollment <- function(n, accrualTime, accrualIntensity) {
  J = length(accrualTime)
  if (J>1) {
    psum = c(0, cumsum(accrualIntensity[1:(J-1)] * diff(accrualTime)))
  } else {
    psum = 0
  }
  rhs =  cumsum(rexp(n))
  j1 = findInterval(rhs, psum)
  accrualTime[j1] + (rhs - psum[j1])/accrualIntensity[j1]
}

arrivalTime <- ceiling(f_enrollment(n, accrualTime, accrualIntensity))

# treatment assignment based on given randomization ratio
f_randomization <- function(n, ngroups, alloc) {
  blocksize = sum(alloc)
  nblocks = ceiling(n/blocksize)
  treats = rep(1:ngroups, alloc)
  treatment = c(replicate(nblocks, sample(treats)))[1:n]
}

treatment <- f_randomization(n, 3, c(2, 2, 1))

# time to treatment discontinuation
muTd = c(5.8, 5.6, 5.7)
sigmaTd = c(0.85, 1.02, 0.73)
time = ceiling(exp(rlogis(n)*sigmaTd[treatment] + muTd[treatment]))

# generate subject-level data and apply cutoff
df2 <- tibble(trialsdt = trialsdt,
             cutoffdt = cutoffdt,
             usubjid = paste0("A-", 100000 + (1:n)),
             arrivalTime = arrivalTime,
             treatment = treatment,
             time = time) %>%
  mutate(treatment_description = ifelse(
    treatment == 1, "Drug A + Drug B",
    ifelse(treatment == 2, "Drug C + Placebo",
           "Drug A + Placebo"))) %>%
  mutate(randdt = as.Date(arrivalTime - 1, origin = trialsdt)) %>%
  mutate(adt = as.Date(time - 1, origin = randdt)) %>%
  filter(randdt <= cutoffdt) %>%
  mutate(event = ifelse(adt <= cutoffdt, 1, 0),
         dropout = 0) %>%
  mutate(time = pmin(time, as.numeric(cutoffdt - randdt + 1))) %>%
  select(trialsdt, usubjid, randdt, treatment,
         treatment_description, time, event, dropout, cutoffdt)


# drug description
drug_description_df = tibble(
  drug = c(1, 2, 3, 4),
  drug_name = c("Drug A", "Drug B", "Drug C", "Placebo"),
  dose_strength = rep("", 4),
  kit = c(1, 2, 3, 4),
  kit_name = c("Drug A", "Drug B", "Drug C", "Placebo"),
  dose_unit = c("kit", "kit", "kit", "kit"),
  p_kit = 1)

l = nrow(drug_description_df)

# treatment by drug combinations
treatment_by_drug = matrix(c(1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1),
                           nrow = 3, ncol = 4, byrow = TRUE)

f_treatment_by_drug_df <- function(
    treatment_by_drug, drug_name) {

  k = nrow(treatment_by_drug)
  l = ncol(treatment_by_drug)

  tibble(treatment = rep(1:k, l),
         drug = rep(1:l, each=k),
         drug_name = rep(drug_name[1:l], each=k),
         included = as.logical(treatment_by_drug)) %>%
    filter(included) %>%
    select(treatment, drug, drug_name)
}

l = nrow(drug_description_df)
drug_name = drug_description_df$drug_name

treatment_by_drug_df <- f_treatment_by_drug_df(
  treatment_by_drug, drug_name)



# function to impute dosing records

pi0 = 0.61
lambda0 = 2.01

muT0 = -0.43
sigmaT0 = 0.27

mu0 = 21.3
sigma0 = 2.60

pi = 0.39
lambda = 2.47

muT = 21.0
sigmaT = 2.84

mud = c(2.0, 3.0, 2.0, 1.0)
sigmab = c(0, 0, 0, 0)
sigmae = c(0.038, 0.042, 0.177, 0.080)


f_draw_t <- function(
    pi0, lambda0, muT0, sigmaT0, mu0, sigma0,
    pi, lambda, muT, sigmaT, vf_new1) {

  df <- tibble()

  for (i in 1:nrow(vf_new1)) {
    Vi = 0
    D = vf_new1$time[i] - 1
    cycle = 0

    z = rbinom(1, 1, pi0)
    if (z == 1) {
      ki = 0
    } else {
      ki = rpois(1, lambda0)
    }

    if (ki == 0) {
      Ti = exp(rlogis(1, muT0, sigmaT0))
    } else {
      Ti = rnorm(1, ki*mu0, sigma0)
    }
    Ti = max(round(Ti), 0)

    while (Vi + Ti <= D) {
      Vi = Vi + Ti

      cycle = cycle + ki + 1

      df = df %>%
        bind_rows(tibble(usubjid = vf_new1$usubjid[i],
                         visit = paste("Cycle", cycle, "Day 1"),
                         day = Vi + 1))

      z = rbinom(1, 1, pi)
      if (z == 1) {
        ki = 0
      } else {
        ki = rpois(1, lambda)
      }

      Ti = rnorm(1, (ki+1)*muT, sigmaT)
      Ti = max(round(Ti), 0)
    }
  }

  df
}


f_draw <- function(
    pi0, lambda0, muT0, sigmaT0, mu0, sigma0,
    pi, lambda, muT, sigmaT, mud, sigmab, sigmae,
    vf_new, vf_new1, treatment_by_drug_df) {

  l = length(unique(treatment_by_drug_df$drug))

  dosing_subject_new1 <- f_draw_t(
    pi0, lambda0, muT0, sigmaT0, mu0, sigma0,
    pi, lambda, muT, sigmaT, vf_new1)

  dosing_subject_new1 <- dosing_subject_new1 %>%
    left_join(vf_new1, by = "usubjid")

  # add drug information for each subject
  dosing_subject_new2 <- purrr::map_dfr(
    1:l, function(h) {
      dosing_subject_new1 %>%
        inner_join(treatment_by_drug_df %>%
                            filter(drug == h),
                          by = "treatment")
    })

  dosing_subject_new <- purrr::map_dfr(
    1:l, function(h) {
      mud1 = mud[[h]]
      sigmab1 = sigmab[[h]]
      sigmae1 = sigmae[[h]]

      df_new <- dosing_subject_new2 %>%
        filter(drug == h)
      n_new = nrow(df_new)

      b1 = rnorm(n_new)*sigmab1
      df_new$dose <- pmax(round(rnorm(n_new)*sigmae1 + mud1 + b1), 1.0)

      df_new
    })

  dosing_subject_new
}


# function to create kit number and dispensed quantity
f_kit <- function(drug, usubjid, visit, dose) {
  if (drug %in% c(1, 2, 3)) {
    kit_number = as.character(round(runif(dose, 5000000, 5310000)))
    dispensed_quantity = 1
  } else {
    kit_number = "None"
    dispensed_quantity = dose
  }

  tibble(drug = drug, usubjid = usubjid, visit = visit,
         kit_number = kit_number,
         dispensed_quantity = dispensed_quantity)
}


# generate drug dispensing data
vf_new <- purrr::map_dfr(
  1:l, function(h) {
    df2 %>%
      inner_join(treatment_by_drug_df %>%
                   filter(drug == h),
                 by = "treatment") %>%
      select(-c(event, dropout))
  })

vf_new1 <- df2 %>% select(-c(event, dropout))


vf <- f_draw(
  pi0, lambda0, muT0, sigmaT0, mu0, sigma0,
  pi, lambda, muT, sigmaT, mud, sigmab, sigmae,
  vf_new, vf_new1, treatment_by_drug_df)

vf1 <- purrr::pmap_dfr(
  list(vf$drug, vf$usubjid, vf$visit, vf$dose), f_kit)


visitview2 <- vf1 %>%
  inner_join(vf, by = c("drug", "usubjid", "visit")) %>%
  mutate(date = as.Date(day - 1, origin = randdt),
         dose_strength = "", dose_unit = "kit") %>%
  select(usubjid, visit, date, drug, drug_name, dose_strength,
         dose_unit, kit_number, dispensed_quantity) %>%
  arrange(usubjid, date, kit_number)


cutoffdt = as.Date("2022-10-01")

df1 <- df2 %>%
  select(-cutoffdt) %>%
  filter(randdt <= cutoffdt) %>%
  mutate(cutoffdt = cutoffdt,
         arrivalTime = as.numeric(randdt - trialsdt + 1),
         followupTime = as.numeric(cutoffdt - randdt + 1),
         event = ifelse(time <= followupTime, event, 0),
         dropout = ifelse(time <= followupTime, dropout, 0),
         time = pmin(time, followupTime)) %>%
  select(-c(arrivalTime, followupTime))


visitview1 <- visitview2 %>%
  filter(date <= cutoffdt)

dosing_schedule_df = dplyr::tibble(
  kit = c(1, 2, 3, 4),
  target_days = c(21, 21, 21, 21),
  target_dose = c(2, 3, 2, 1),
  max_cycles = c(10000, 10000, 10000, 10000))

# save to data/ folder

usethis::use_data(drug_description_df, overwrite = TRUE)
usethis::use_data(treatment_by_drug, overwrite = TRUE)
usethis::use_data(dosing_schedule_df, overwrite = TRUE)

usethis::use_data(df1, overwrite = TRUE)
usethis::use_data(df2, overwrite = TRUE)
usethis::use_data(visitview1, overwrite = TRUE)
usethis::use_data(visitview2, overwrite = TRUE)
