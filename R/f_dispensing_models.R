#' @title Model Fitting for Dispensing Delay After Randomization
#' @description Fits a specified time-to-event model to the gap time
#' between randomization and the first drug dispensing visit when
#' there is no visit skipping.
#'
#' @param df The subject-level dosing data, including the following
#'   variables:
#'
#'   * \code{time}: The number of days between randomization and the
#'     first drug dispensing visit (first drug dispensing visit date -
#'     randomization date + 1).
#'
#'   * \code{left}: Equals \code{time - 1}, used to indicate the
#'     left endpoint of an interval for interval censoring.
#'
#'   * \code{right}: Equals \code{time}, used to indicate the
#'     right endpoint of an interval for interval censoring.
#'
#' @param model The event model used to analyze the gap time
#'   between randomization and the first drug dispensing visit when
#'   there is no visit skipping, with options including "constant",
#'   "exponential", "weibull", "log-logistic", and "log-normal".
#' @param nreps The number of simulations for drawing posterior model
#'   parameter values.
#' @param showplot A Boolean variable that controls whether or not to
#'   show the fitted time-to-event bar chart. It defaults to \code{TRUE}.
#'
#' @return A list with three components:
#' * \code{fit}: A list of results from the model fit that includes
#'
#'     + \code{model}: The specific model used in the analysis.
#'
#'     + \code{theta}: The estimated model parameters.
#'
#'     + \code{vtheta}: The estimated covariance matrix of \code{theta}.
#'
#'     + \code{aic}: The Akaike Information Criterion value.
#'
#'     + \code{bic}: The Bayesian Information Criterion value.
#'
#' * \code{fit_plot}: A fitted time-to-event bar chart.
#'
#' * \code{theta}: Posterior draws of model parameters.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' library(dplyr)
#'
#' df <- df2 %>%
#'   mutate(arrivalTime = as.numeric(randdt - trialsdt + 1))
#'
#' vf <- visitview2 %>%
#'   inner_join(df, by = "usubjid") %>%
#'   mutate(day = as.numeric(date - randdt + 1)) %>%
#'   select(drug, drug_name, dose_unit, usubjid, treatment,
#'          treatment_description, arrivalTime,
#'          time, event, dropout, day, dispensed_quantity) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid, treatment,
#'            treatment_description, arrivalTime,
#'            time, event, dropout, day) %>%
#'   summarise(dose = sum(dispensed_quantity), .groups = "drop_last") %>%
#'   mutate(cum_dose = cumsum(dose)) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid) %>%
#'   mutate(row_id = row_number())
#'
#' vf <- vf %>%
#'   left_join(dosing_schedule_df, by = "drug")
#'
#' vf1 <- vf %>%
#'   group_by(usubjid, day) %>%
#'   slice(n()) %>%
#'   group_by(usubjid) %>%
#'   mutate(row_id = row_number())
#'
#' # time from randomization to the first drug dispensing visit
#' df_k0 <- vf1 %>%
#'   filter(row_id == 1) %>%
#'   mutate(time = day,
#'          skipped = floor((time - target_days/2)/target_days) + 1)
#'
#' # no skipping
#' df_t0 <- df_k0 %>%
#'   filter(skipped == 0) %>%
#'   mutate(left = time - 1, right = time)
#'
#' fit_t0 <- f_fit_t0(df_t0, model = "log-logistic", nreps = 200)
#'
#' @export
f_fit_t0 <- function(df, model, nreps, showplot = TRUE) {
  erify::check_class(df, "data.frame")
  names(df) <- tolower(names(df))
  n = nrow(df)

  model = tolower(model)
  erify::check_content(model, c("constant", "exponential", "weibull",
                                "log-logistic", "log-normal"))

  erify::check_n(nreps)
  erify::check_bool(showplot)

  # observed probabilities
  x = table(df$time - 1)
  count = as.numeric(x)
  y.obs = as.numeric(names(x))
  y = min(y.obs):max(y.obs)
  p.obs = rep(0, length(y))
  p.obs[findInterval(y.obs, y)] = count/sum(count)

  # fit time-to-event models
  l = length(unique(df$time))
  if (l == 1) {
    fit <- list(model = "Constant",
                theta = df$time[1] - 1,
                vtheta = 0,
                aic = NA,
                bic = NA)

    p.fit = 1

    post = rep(fit$theta, nreps)
  } else {
    if (l == 2 & model != "exponential") {
      model = "exponential"
    }

    if (model == "exponential") {
      a <- icenReg::ic_par(cbind(left, right) ~ 1, data = df,
                           model = "ph", dist = "exponential")

      fit <- list(model = "Exponential",
                  theta = -as.numeric(a$coefficients),  # log(rate)
                  vtheta = as.numeric(a$var),
                  aic = -2*a$llk + 2,
                  bic = -2*a$llk + log(n))

      rate = exp(fit$theta)
      p.fit = pexp(y+1, rate) - pexp(y, rate)

      post = rnorm(nreps, mean = fit$theta, sd = sqrt(fit$vtheta))
    } else if (model == "weibull") {
      a <- icenReg::ic_par(cbind(left, right) ~ 1, data = df,
                           model = "ph", dist = "weibull")

      # convert to (log(lambda) = log_scale, -log(kappa) = -log_shape)
      const = matrix(c(0, 1, -1, 0), 2, 2, byrow = TRUE)
      fit <- list(model = "Weibull",
                  theta = c(const %*% a$coefficients),
                  vtheta = const %*% a$var %*% t(const),
                  aic = -2*a$llk + 4,
                  bic = -2*a$llk + 2*log(n))

      lambda = exp(fit$theta[1])
      k = exp(-fit$theta[2])
      p.fit = pweibull(y+1, k, lambda) - pweibull(y, k, lambda)

      post = mvtnorm::rmvnorm(nreps, mean = fit$theta, sigma = fit$vtheta)
    } else if (model == "log-logistic") {
      a <- icenReg::ic_par(cbind(left, right) ~ 1, data = df,
                           model = "ph", dist = "loglogistic")

      # convert to (mu = log_scale, log(sigma) = -log_shape)
      const = matrix(c(1, 0, 0, -1), 2, 2, byrow = TRUE)
      fit <- list(model = "Log-logistic",
                  theta = c(const %*% a$coefficients),
                  vtheta = const %*% a$var %*% t(const),
                  aic = -2*a$llk + 4,
                  bic = -2*a$llk + 2*log(n))

      mu = fit$theta[1]
      sigma = exp(fit$theta[2])
      p.fit = plogis(log(y+1), mu, sigma) - plogis(log(y), mu, sigma)

      post = mvtnorm::rmvnorm(nreps, mean = fit$theta, sigma = fit$vtheta)
    } else if (model == "log-normal") {
      a <- icenReg::ic_par(cbind(left, right) ~ 1, data = df,
                           model = "ph", dist = "lnorm")

      # parametrization: (mu = mean.log, log(sigma) = log(sd.log))
      fit <- list(model = "Log-normal",
                  theta = as.numeric(a$coefficients),
                  vtheta = a$var,
                  aic = -2*a$llk + 4,
                  bic = -2*a$llk + 2*log(n))

      mu = fit$theta[1]
      sigma = exp(fit$theta[2])
      p.fit = plnorm(y+1, mu, sigma) - plnorm(y, mu, sigma)

      post = mvtnorm::rmvnorm(nreps, mean = fit$theta, sigma = fit$vtheta)
    } else {
      stop("incorrect model for T0")
    }
  }


  # graphically assess the model fit
  modeltext = fit$model
  if (modeltext != "Constant") {
    aictext = paste("AIC:", round(fit$aic,2))
    bictext = paste("BIC:", round(fit$bic,2))
  } else {
    modeltext = ""
    aictext = ""
    bictext = ""
  }

  gf <- dplyr::tibble(y = y, p.obs = p.obs, p.fit = p.fit)

  fig <- plotly::plot_ly(gf, x = ~y, y = ~p.obs, type = 'bar',
                         name = 'Observed')
  fig <- fig %>% plotly::add_trace(y = ~p.fit, name = 'Fitted')
  fig <- fig %>% plotly::layout(
    xaxis = list(title = paste("Days between randomization and the",
                               "first drug dispensing visit")),
    yaxis = list(title = 'Proportion'),
    barmode = 'group')

  fig <- fig %>% plotly::layout(
    annotations = list(
      x = c(0.7, 0.7, 0.7), y = c(0.95, 0.80, 0.65),
      xref = "paper", yref = "paper",
      text = paste('<i>', c(modeltext, aictext, bictext), '</i>'),
      xanchor = "left", font = list(size = 14, color = "red"),
      showarrow = FALSE))

  if (showplot) print(fig)

  list(fit = fit, fit_plot = fig, theta = post)
}


#' @title Model Fitting for Number of Skipped Visits
#' @description Fits a count model to the number of skipped visits
#' between two consecutive drug dispensing visits.
#'
#' @param df The subject-level dosing data, including \code{skipped} to
#'   indicate the number of skipped visits.
#' @param model The count model used to analyze the number of
#'   skipped visits, with options including
#'   "constant", "poisson", "zero-inflated poisson", and
#'   "negative binomial".
#' @param nreps The number of simulations for drawing posterior model
#'   parameter values.
#' @param showplot A Boolean variable that controls whether or not to
#'   show the fitted count bar chart. It defaults to \code{TRUE}.
#'
#' @return A list with three components:
#'
#' * \code{fit}: A list of results from the model fit that includes
#'
#'     + \code{model}: The specific model used in the analysis.
#'
#'     + \code{theta}: The estimated model parameters.
#'
#'     + \code{vtheta}: The estimated covariance matrix of \code{theta}.
#'
#'     + \code{aic}: The Akaike Information Criterion value.
#'
#'     + \code{bic}: The Bayesian Information Criterion value.
#'
#' * \code{fit_plot}: A fitted count bar chart.
#'
#' * \code{theta}: Posterior draws of model parameters.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' library(dplyr)
#'
#' df <- df2 %>%
#' mutate(arrivalTime = as.numeric(randdt - trialsdt + 1))
#'
#' vf <- visitview2 %>%
#'   inner_join(df, by = "usubjid") %>%
#'   mutate(day = as.numeric(date - randdt + 1)) %>%
#'   select(drug, drug_name, dose_unit, usubjid, treatment,
#'          treatment_description, arrivalTime,
#'          time, event, dropout, day, dispensed_quantity) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid, treatment,
#'            treatment_description, arrivalTime,
#'            time, event, dropout, day) %>%
#'   summarise(dose = sum(dispensed_quantity), .groups = "drop_last") %>%
#'   mutate(cum_dose = cumsum(dose)) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid) %>%
#'   mutate(row_id = row_number())
#'
#' vf <- vf %>%
#'   left_join(dosing_schedule_df, by = "drug")
#'
#' vf1 <- vf %>%
#'   group_by(usubjid, day) %>%
#'   slice(n()) %>%
#'   group_by(usubjid) %>%
#'   mutate(row_id = row_number())
#'
#' df_ti <- vf1 %>%
#'   mutate(time = lead(day) - day,
#'          skipped = pmax(floor((time - target_days/2)/target_days), 0),
#'          k1 = skipped + 1) %>%
#'   filter(row_id < n())
#'
#' fit_ki <- f_fit_ki(df_ti, model = "zero-inflated poisson", nreps = 200)
#'
#' @export
f_fit_ki <- function(df, model, nreps, showplot = TRUE) {
  erify::check_class(df, "data.frame")
  names(df) <- tolower(names(df))
  n = nrow(df)

  model = tolower(model)
  erify::check_content(model, c("constant", "poisson",
                                "zero-inflated poisson",
                                "negative binomial"))

  erify::check_n(nreps)
  erify::check_bool(showplot)

  # observed probabilities
  x = table(df$skipped)
  count = as.numeric(x)
  y.obs = as.numeric(names(x))
  y = min(y.obs):max(y.obs)
  p.obs = rep(0, length(y))
  p.obs[findInterval(y.obs, y)] = count/sum(count)

  # fit the count model
  l = length(unique(df$skipped))
  if (l == 1) {
    fit <- list(model = "Constant",
                theta = df$skipped[1],
                vtheta = 0,
                aic = NA,
                bic = NA)

    p.fit = 1

    post = rep(fit$theta, nreps)
  } else {
    if (l == 2 & model != "poisson") {
      model = "poisson"
    }

    if (model == "poisson") {
      a <- glm(skipped ~ 1, family = poisson(link = "log"), data = df)

      fit <- list(model = "Poisson",
                  theta = as.numeric(a$coefficients),  # log(rate)
                  vtheta = as.numeric(vcov(a)),
                  aic = a$aic,
                  bic = a$aic - 2 + log(n))

      rate = exp(fit$theta)
      p.fit = dpois(y, rate)

      post = rnorm(nreps, mean = fit$theta, sd = sqrt(fit$vtheta))
    } else if (model == "zero-inflated poisson") {
      a <- pscl::zeroinfl(skipped ~ 1 | 1, data = df, dist = "poisson")

      fit <- list(model = "Zero-inflated Poisson",
                  theta = as.numeric(a$coefficients),
                  vtheta = vcov(a),
                  aic = -2*a$loglik + 4,
                  bic = -2*a$loglik + 2*log(n))

      lambda = exp(fit$theta[1])
      pi = plogis(fit$theta[2])

      p.fit <- c(pi + (1-pi)*dpois(0, lambda), (1-pi)*dpois(y[-1], lambda))

      post = mvtnorm::rmvnorm(nreps, mean = fit$theta, sigma = fit$vtheta)
    } else if (model == "negative binomial") {
      a <- MASS::glm.nb(skipped ~ 1, data = df)

      # parametrization: log(mean), log(size)
      fit <- list(model = "Negative binomial",
                  theta = c(as.numeric(a$coefficients), log(a$theta)),
                  vtheta = diag(c(as.numeric(vcov(a)),
                                  a$SE.theta^2/a$theta^2)),
                  aic = a$aic,
                  bic = -a$twologlik + 2*log(n))

      mu = exp(fit$theta[1])
      size = exp(fit$theta[2])
      prob = size/(size + mu)
      p.fit = dnbinom(y, size, prob)

      post = mvtnorm::rmvnorm(nreps, mean = fit$theta, sigma = fit$vtheta)
    } else {
      stop("incorrect model for ki")
    }
  }


  # graphically assess the model fit
  modeltext = fit$model
  if (modeltext != "Constant") {
    aictext = paste("AIC:", round(fit$aic,2))
    bictext = paste("BIC:", round(fit$bic,2))
  } else {
    modeltext = ""
    aictext = ""
    bictext = ""
  }

  gf <- dplyr::tibble(y = y, p.obs = p.obs, p.fit = p.fit)

  fig <- plotly::plot_ly(gf, x = ~y, y = ~p.obs, type = 'bar',
                         name = 'Observed')
  fig <- fig %>% plotly::add_trace(y = ~p.fit, name = 'Fitted')
  fig <- fig %>% plotly::layout(
    xaxis = list(title = 'Number of skipped visits'),
    yaxis = list(title = 'Proportion'),
    barmode = 'group')

  fig <- fig %>% plotly::layout(
    annotations = list(
      x = c(0.7, 0.7, 0.7), y = c(0.95, 0.80, 0.65),
      xref = "paper", yref = "paper",
      text = paste('<i>', c(modeltext, aictext, bictext), '</i>'),
      xanchor = "left", font = list(size = 14, color = "red"),
      showarrow = FALSE))

  if (showplot) print(fig)

  list(fit = fit, fit_plot = fig, theta = post)
}


#' @title Model Fitting for Gap Times
#' @description Fits a linear regression model to the gap time
#' between two consecutive drug dispensing visits.
#'
#' @param df The subject-level dosing data, including the following
#'   variables:
#'
#'   * \code{time}: The gap time to the next drug dispensing visit.
#'
#'   * \code{skipped}: The number of skipped visits.
#'
#'   * \code{k1}: The covariate for the linear regression. It equals
#'     \code{skipped} for the gap time between randomization and
#'     the first drug dispensing visit and \code{skipped + 1}
#'     for the gap time between two consecutive drug dispensing visits.
#' @param model The model used to analyze the gap time. Options include
#'   "least squares" and "least absolute deviations".
#' @param nreps The number of simulations for drawing posterior model
#'   parameter values.
#' @param showplot A Boolean variable that controls whether or not to
#'   show the fitted gap time bar chart. It defaults to \code{TRUE}.
#'
#' @return A list with three components:
#'
#' * \code{fit}: A list of results from the model fit that includes
#'
#'     + \code{model}: The specific model used in the analysis.
#'
#'     + \code{beta}: The estimated regression coefficient for the covariate.
#'
#'     + \code{vbeta}: The estimated variance of \code{beta}.
#'
#'     + \code{sigma}: The estimated residual standard deviation.
#'
#'     + \code{df}: The residual degrees-of-freedom.
#'
#'     + \code{aic}: The Akaike Information Criterion value.
#'
#'     + \code{bic}: The Bayesian Information Criterion value.
#'
#' * \code{fit_plot}: A fitted gap time bar chart.
#'
#' * \code{theta}: Posterior draws of model parameters.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' library(dplyr)
#'
#' df <- df2 %>%
#' mutate(arrivalTime = as.numeric(randdt - trialsdt + 1))
#'
#' vf <- visitview2 %>%
#'   inner_join(df, by = "usubjid") %>%
#'   mutate(day = as.numeric(date - randdt + 1)) %>%
#'   select(drug, drug_name, dose_unit, usubjid, treatment,
#'          treatment_description, arrivalTime,
#'          time, event, dropout, day, dispensed_quantity) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid, treatment,
#'            treatment_description, arrivalTime,
#'            time, event, dropout, day) %>%
#'   summarise(dose = sum(dispensed_quantity), .groups = "drop_last") %>%
#'   mutate(cum_dose = cumsum(dose)) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid) %>%
#'   mutate(row_id = row_number())
#'
#' vf <- vf %>%
#'   left_join(dosing_schedule_df, by = "drug")
#'
#' vf1 <- vf %>%
#'   group_by(usubjid, day) %>%
#'   slice(n()) %>%
#'   group_by(usubjid) %>%
#'   mutate(row_id = row_number())
#'
#' df_ti <- vf1 %>%
#'   mutate(time = lead(day) - day,
#'          skipped = pmax(floor((time - target_days/2)/target_days), 0),
#'          k1 = skipped + 1) %>%
#'   filter(row_id < n())
#'
#' fit_ti <- f_fit_ti(df_ti, model = "least squares", nreps = 200)
#'
#' @export
f_fit_ti <- function(df, model, nreps, showplot = TRUE) {
  erify::check_class(df, "data.frame")
  names(df) <- tolower(names(df))

  model = tolower(model)
  erify::check_content(model, c("least squares",
                                "least absolute deviations"))

  erify::check_n(nreps)
  erify::check_bool(showplot)

  x = table(df$time)
  count = as.numeric(x)
  y.obs = as.numeric(names(x))
  y = min(y.obs):max(y.obs)
  p.obs = rep(0, length(y))
  p.obs[findInterval(y.obs, y)] = count/sum(count)

  if (nrow(df) == 1) {
    fit <- list(model = "Constant",
                beta = df$time[1]/df$k1[1],
                vbeta = 0,
                sigma = 0,
                df = 0,
                aic = NA,
                bic = NA)

    p.fit = 1
  } else if (model == "least squares") {
    a <- lm(time ~ k1 - 1, data = df)

    fit <- list(model = "Least squares",
                beta = as.numeric(a$coefficients),
                vbeta = as.numeric(vcov(a)),
                sigma = summary(a)$sigma,
                df = a$df.residual,
                aic = as.numeric(AIC(a)),
                bic = as.numeric(BIC(a)))

    p.fit <- purrr::map_vec(y, function(y) {
      mean(pnorm((y + 0.5 - fit$beta*df$k1)/fit$sigma) -
             pnorm((y - 0.5 - fit$beta*df$k1)/fit$sigma))
    })
  } else if (model == "least absolute deviations") {
    a <- L1pack::lad(time ~ k1 - 1, data = df)

    fit <- list(model = "Least absolute deviations",
                beta = as.numeric(a$coefficients),
                vbeta = as.numeric(vcov(a)),
                sigma = a$scale,
                df = a$dims[1] - a$dims[2],
                aic = as.numeric(AIC(a)),
                bic = as.numeric(BIC(a)))

    p.fit <- purrr::map_vec(y, function(y) {
      d <- df %>%
        dplyr::mutate(x1 = y+0.5-fit$beta*.data$k1,
                      x2 = y-0.5-fit$beta*.data$k1,
                      t1 = exp(.data$x1/fit$sigma),
                      t2 = exp(.data$x2/fit$sigma)) %>%
        dplyr::mutate(prob = ifelse(
          .data$x1 < 0, 0.5*(.data$t1 - .data$t2),
          ifelse(.data$x2 > 0, 0.5*(-1/.data$t1 + 1/.data$t2),
                 1 - 0.5*(1/.data$t1 + .data$t2))))

      mean(d$prob)
    })
  }

  gf = dplyr::tibble(y = y, p.obs = p.obs, p.fit = p.fit)

  modeltext = fit$model
  if (modeltext != "Constant") {
    aictext = paste("AIC:", round(fit$aic,2))
    bictext = paste("BIC:", round(fit$bic,2))
  } else {
    modeltext = ""
    aictext = ""
    bictext = ""
  }

  fig <- plotly::plot_ly(gf, x = ~y, y = ~p.obs, type = 'bar',
                         name = 'Observed')
  fig <- fig %>% plotly::add_trace(y = ~p.fit, name = 'Fitted')
  fig <- fig %>% plotly::layout(
    xaxis = list(title = 'Days between consecutive drug dispensing visits'),
    yaxis = list(title = 'Proportion'),
    barmode = 'group')

  fig <- fig %>% plotly::layout(
    annotations = list(
      x = c(0.7, 0.7, 0.7), y = c(0.95, 0.80, 0.65),
      xref = "paper", yref = "paper",
      text = paste('<i>', c(modeltext, aictext, bictext), '</i>'),
      xanchor = "left", font = list(size = 14, color = "red"),
      showarrow = FALSE))

  if (showplot) print(fig)


  if (fit$df == 0) {
    post = matrix(c(rep(fit$beta, nreps), rep(0, nreps)), nreps, 2)
  } else {
    # draw sigma and then beta given sigma from posterior
    b2 = sqrt(fit$df * fit$sigma^2 / rchisq(nreps, fit$df))
    b1 = fit$beta + rnorm(nreps) * sqrt(fit$vbeta) / fit$sigma * b2
    post = matrix(c(b1, b2), nreps, 2)
  }

  list(fit = fit, fit_plot = fig, theta = post)
}


#' @title Model Fitting for Dispensed Doses
#' @description Fits a linear mixed-effects model to the dispensed doses
#' at drug dispensing visits.
#'
#' @param df The subject-level dosing data, including \code{usubjid},
#'   \code{day}, \code{drug}, and \code{dose}.
#' @param model The model used to analyze the dispensed doses, with
#'   options including "constant", "linear model", and
#'   "linear mixed-effects model".
#' @param nreps The number of simulations for drawing posterior model
#'   parameters.
#' @param showplot A Boolean variable that controls whether or not to
#'   show the fitted dose bar chart. It defaults to \code{TRUE}.
#'
#' @return A list with three components:
#'
#' * \code{fit}: A list of results from the model fit that includes
#'
#'     + \code{model}: The specific model used in the analysis.
#'
#'     + \code{mud}: The estimated mean dose.
#'
#'     + \code{vmud}: The estimated variance of \code{mud}.
#'
#'     + \code{sigmab}: The estimated between-subject standard deviation.
#'
#'     + \code{sigmae}: The estimated within-subject residual standard
#'       deviation.
#'
#'     + \code{aic}: The Akaike Information Criterion value.
#'
#'     + \code{bic}: The Bayesian Information Criterion value.
#'
#' * \code{fit_plot}: A fitted dose bar chart.
#'
#' * \code{theta}: Posterior draws of model parameters.
#'
#'     + \code{fixed}: Posterior draws of fixed model parameters:
#'       \code{mud}, \code{sigmab}, and \code{sigmae}.
#'
#'     + \code{random}: Posterior draws of subject random effects.
#'
#'     + \code{usubjid}: The unique subject ID associated with
#'       the subject random effects.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' library(dplyr)
#'
#' df <- df2 %>%
#'   mutate(arrivalTime = as.numeric(randdt - trialsdt + 1))
#'
#' vf <- visitview2 %>%
#'   inner_join(df, by = "usubjid") %>%
#'   mutate(day = as.numeric(date - randdt + 1)) %>%
#'   select(drug, drug_name, dose_unit, usubjid, treatment,
#'          treatment_description, arrivalTime,
#'          time, event, dropout, day, dispensed_quantity) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid, treatment,
#'            treatment_description, arrivalTime,
#'            time, event, dropout, day) %>%
#'   summarise(dose = sum(dispensed_quantity), .groups = "drop_last") %>%
#'   mutate(cum_dose = cumsum(dose)) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid) %>%
#'   mutate(row_id = row_number())
#'
#' vf1 <- vf %>% filter(drug == 3)
#' fit_di <- f_fit_di(vf1, model = "linear mixed-effects model", nreps = 200)
#'
#' @export
f_fit_di <- function(df, model, nreps, showplot = TRUE) {
  erify::check_class(df, "data.frame")
  names(df) <- tolower(names(df))
  n = nrow(df)                   # total number of observations
  N = length(unique(df$usubjid)) # number of unique subjects

  model = tolower(model)
  erify::check_content(model, c("constant", "linear model",
                                "linear mixed-effects model"))

  erify::check_n(nreps)
  erify::check_bool(showplot)

  # observed probabilities
  x = table(df$dose)
  count = as.numeric(x)
  y.obs = as.numeric(names(x))
  y = min(y.obs):max(y.obs)
  p.obs = rep(0, length(y))
  p.obs[findInterval(y.obs, y)] = count/sum(count)

  l = length(unique(df$dose))
  if (l == 1) {
    fit <- list(model = "Constant",
                mud = df$dose[1],
                vmud = 0,
                sigmab = 0,
                sigmae = 0,
                aic = NA,
                bic = NA)

    p.fit = 1

    gf = dplyr::tibble(y = y, p.obs = p.obs, p.fit = p.fit)

    theta_fix = matrix(c(rep(fit$mud, nreps), rep(0, 2*nreps)), nreps, 3)
    theta_ran = matrix(0, nreps, N)
  } else {
    if (l == 2 & model != "linear model") {
      model = "linear model"
    }

    if (model == "linear model") {
      a <- lm(dose ~ 1, data = df)

      fit <- list(model = "Linear model",
                  mud = as.numeric(a$coefficients),
                  vmud = as.numeric(vcov(a)),
                  sigmab = 0,
                  sigmae = summary(a)$sigma,
                  aic = as.numeric(AIC(a)),
                  bic = as.numeric(BIC(a)))

      p.fit = pnorm(y + 0.5, fit$mud, fit$sigmae) -
        pnorm(y - 0.5, fit$mud, fit$sigmae)

      gf = dplyr::tibble(y = y, p.obs = p.obs, p.fit = p.fit)

      # draw sigmae and then mud given sigmae from posterior
      b2 = sqrt(a$df.residual * fit$sigmae^2 / rchisq(nreps, a$df.residual))
      b1 = fit$mud + rnorm(nreps)*sqrt(fit$vmud)/fit$sigmae*b2
      theta_fix = matrix(c(b1, rep(0, nreps), b2), nreps, 3)
      theta_ran = matrix(0, nreps, N)
    } else if (model == "linear mixed-effects model") {
      a <- nlme::lme(dose ~ 1, random = ~ 1 | usubjid, data = df)

      fit <- list(model = "Linear mixed-effects model",
                  mud = as.numeric(a$coefficients$fixed),
                  vmud = as.numeric(vcov(a)),
                  sigmab = exp(as.numeric(attr(a$apVar, "Pars")))[1],
                  sigmae = exp(as.numeric(attr(a$apVar, "Pars")))[2],
                  aic = as.numeric(AIC(a)),
                  bic = as.numeric(BIC(a)))

      est = as.numeric(a$fitted[,"usubjid"])
      sigmae = fit$sigmae

      p.fit <- purrr::map_vec(y, function(y) {
        mean(pnorm((y + 0.5 - est)/sigmae) - pnorm((y - 0.5 - est)/sigmae))
      })

      gf = dplyr::tibble(y = y, p.obs = p.obs, p.fit = p.fit)

      # number of observations and mean dose by subject
      df1 <- df %>%
        dplyr::summarise(n = dplyr::n(), d = mean(.data$dose),
                         .groups = "drop_last")

      # fix the variance parameters to avoid drawing too extreme values
      # of sigma_b due to the large variance of log(sigma_b) in case that
      # the variance component is likely to be zero
      sigma_b = rep(fit$sigmab, nreps)
      sigma_e = rep(fit$sigmae, nreps)
      sigma_b2 = sigma_b^2
      sigma_e2 = sigma_e^2

      mud = fit$mud + sqrt(fit$vmud)*rnorm(nreps)

      # fixed effects model parameters
      theta_fix = matrix(c(mud, sigma_b, sigma_e), nreps, 3)

      # random effects model parameters
      theta_ran = matrix(0, nreps, N)
      for (i in 1:N) {
        si1 = df1$n[i]/sigma_e2*(df1$d[i] - mud)
        si2 = df1$n[i]/sigma_e2 + 1/sigma_b2
        theta_ran[,i] = si1/si2 + sqrt(1/si2)*rnorm(nreps)
      }
    } else {
      stop("incorrect model for Ti")
    }
  }

  # unique subject id corresponding to the random effects
  usubjid = unique(df$usubjid)

  post = list(fixed = theta_fix, random = theta_ran, usubjid = usubjid)

  # graphically assess the model fit
  modeltext = fit$model
  if (modeltext != "Constant") {
    aictext = paste("AIC:", round(fit$aic,2))
    bictext = paste("BIC:", round(fit$bic,2))
  } else {
    modeltext = ""
    aictext = ""
    bictext = ""
  }

  fig <- plotly::plot_ly(gf, x = ~y, y = ~p.obs, type = 'bar',
                         name = 'Observed')
  fig <- fig %>% plotly::add_trace(y = ~p.fit, name = 'Fitted')
  fig <- fig %>% plotly::layout(
    xaxis = list(title = 'Dose dispensed at drug dispensing visits'),
    yaxis = list(title = 'Proportion'),
    barmode = 'group')

  fig <- fig %>% plotly::layout(
    annotations = list(
      x = c(0.7, 0.7, 0.7), y = c(0.95, 0.80, 0.65),
      xref = "paper", yref = "paper",
      text = paste('<i>', c(modeltext, aictext, bictext), '</i>'),
      xanchor = "left", font = list(size = 14, color = "red"),
      showarrow = FALSE))

  if (showplot) print(fig)

  list(fit = fit, fit_plot = fig, theta = post)
}


#' @title Drug Dispensing Model Fitting
#' @description Fits drug dispensing models to the observed drug
#' dispensing data.
#'
#' @param vf A data frame for subject-level drug dispensing data,
#'   including the following variables:
#'   \code{drug}, \code{drug_name}, \code{dose_unit},
#'   \code{usubjid}, \code{treatment}, \code{treatment_description},
#'   \code{arrivalTime}, \code{time}, \code{event}, \code{dropout},
#'   \code{day}, \code{dose}, \code{cum_dose}, and \code{row_id}.
#' @param dosing_schedule_df A data frame providing dosing schedule
#'   information. It contains the following variables: \code{drug},
#'   \code{target_days}, \code{target_dose}, and \code{max_cycles}.
#' @param model_k0 The model for the number of skipped
#'   visits between randomization and the first drug dispensing visit.
#'   Options include "constant", "poisson", "zero-inflated poisson",
#'   and "negative binomial".
#' @param model_t0 The model for the gap time between randomization
#'   and the first drug dispensing visit when there is no visit skipping.
#'   Options include "constant", "exponential", "weibull",
#'   "log-logistic", and "log-normal".
#' @param model_t1 The model for the gap time between randomization
#'   and the first drug dispensing visit when there is visit skipping.
#'   Options include "least squares", and "least absolute deviations".
#' @param model_ki The model for the number of skipped
#'   visits between two consecutive drug dispensing visits.
#'   Options include "constant", "poisson", "zero-inflated poisson",
#'   and "negative binomial".
#' @param model_ti The model for the gap time between two consecutive
#'   drug dispensing visits. Options include "least squares"
#'   and "least absolute deviations".
#' @param model_di The model for the dispensed doses at drug
#'   dispensing visits. Options include "constant",
#'   "linear model", and "linear mixed-effects model".
#' @param nreps The number of simulations for drawing posterior model
#'   parameters.
#' @param showplot A Boolean variable that controls whether or not to
#'   show the model fit plot. It defaults to \code{TRUE}.
#'
#' @return A list with the following components:
#'
#' * \code{common_time_model}: A Boolean variable that indicates
#'   whether a common time model is used for drug dispensing visits.
#'
#' * \code{fit_k0}: The model fit for the number of skipped
#'   visits between randomization and the first drug dispensing visit.
#'
#' * \code{fit_t0}: The model fit for the gap time between
#'   randomization and the first drug dispensing visit when there is
#'   no visit skipping.
#'
#' * \code{fit_t1}: The model fit for the gap time between
#'   randomization and the first drug dispensing visit when there is
#'   visit skipping.
#'
#' * \code{fit_ki}: The model fit for the number of skipped
#'   visits between two consecutive drug dispensing visits.
#'
#' * \code{fit_ti}: The model fit for the gap time between two
#'   consecutive drug dispensing visits.
#'
#' * \code{fit_di}: The model fit for the dispensed doses at drug
#'   dispensing visits.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @seealso \code{\link{f_fit_t0}}, \code{\link{f_fit_ki}},
#' \code{\link{f_fit_ti}}, \code{\link{f_fit_di}}
#'
#' @examples
#' library(dplyr)
#'
#' df <- df2 %>%
#'   mutate(arrivalTime = as.numeric(randdt - trialsdt + 1))
#'
#' vf <- visitview2 %>%
#'   inner_join(df, by = "usubjid") %>%
#'   mutate(day = as.numeric(date - randdt + 1)) %>%
#'   select(drug, drug_name, dose_unit, usubjid, treatment,
#'          treatment_description, arrivalTime,
#'          time, event, dropout, day, dispensed_quantity) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid, treatment,
#'            treatment_description, arrivalTime,
#'            time, event, dropout, day) %>%
#'   summarise(dose = sum(dispensed_quantity), .groups = "drop_last") %>%
#'   mutate(cum_dose = cumsum(dose)) %>%
#'   group_by(drug, drug_name, dose_unit, usubjid) %>%
#'   mutate(row_id = row_number())
#'
#' dispensing_models <- f_dispensing_models(
#'   vf, dosing_schedule_df,
#'   model_k0 = "zero-inflated poisson",
#'   model_t0 = "log-logistic", model_t1 = "least squares",
#'   model_ki = "zero-inflated poisson", model_ti = "least squares",
#'   model_di = "linear mixed-effects model",
#'   nreps = 200, showplot = FALSE)
#'
#' dispensing_models$fit_ki$fit_plot
#
#' @export
f_dispensing_models <- function(
    vf, dosing_schedule_df,
    model_k0, model_t0, model_t1, model_ki, model_ti, model_di,
    nreps, showplot = TRUE) {

  erify::check_content(tolower(model_k0),
                       c("constant", "poisson", "zero-inflated poisson",
                         "negative binomial"))
  erify::check_content(tolower(model_t0),
                       c("constant", "exponential", "weibull",
                         "log-logistic", "log-normal"))
  erify::check_content(tolower(model_t1),
                       c("least squares", "least absolute deviations"))
  erify::check_content(tolower(model_ki),
                       c("constant", "poisson", "zero-inflated poisson",
                         "negative binomial"))
  erify::check_content(tolower(model_ti),
                       c("least squares", "least absolute deviations"))
  erify::check_content(tolower(model_di),
                       c("constant", "linear model",
                         "linear mixed-effects model"))

  vf = vf %>% left_join(dosing_schedule_df, by = "drug")
  l = nrow(dosing_schedule_df)

  if (length(unique(dosing_schedule_df$target_days)) == 1) {
    # fit a common model for k0, t0, t1, ki, and ti across drugs
    common_time_model = TRUE

    # only keep one record per subject and drug dispensing day
    # need to redefine row_id due to drug dispensing at unscheduled visits
    vf1 <- vf %>%
      dplyr::group_by(.data$usubjid, .data$day) %>%
      dplyr::slice(dplyr::n()) %>%
      dplyr::group_by(.data$usubjid) %>%
      dplyr::mutate(row_id = dplyr::row_number())

    # time from randomization to the first drug dispensing visit
    df_k0 <- vf1 %>%
      dplyr::filter(.data$row_id == 1) %>%
      dplyr::mutate(
        time = .data$day,
        skipped = floor((.data$time - .data$target_days/2)/
                          .data$target_days) + 1)

    fit_k0 <- f_fit_ki(df_k0, model_k0, nreps, showplot)

    # no skipping
    df_t0 <- df_k0 %>%
      dplyr::filter(.data$skipped == 0) %>%
      dplyr::mutate(left = .data$time - 1, right = .data$time)

    fit_t0 <- f_fit_t0(df_t0, model_t0, nreps, showplot)

    # skipping
    df_t1 <- df_k0 %>%
      dplyr::filter(.data$skipped > 0) %>%
      dplyr::mutate(k1 = .data$skipped)

    if (nrow(df_t1) == 0) {
      fit_t1 <- list(fit = list(model = "Linear model",
                                beta = 0,
                                vbeta = 0,
                                sigma = 0,
                                df = NA,
                                aic = NA,
                                bic = NA),
                     fit_plot = NA,
                     theta = matrix(0, nreps, 2))
    } else {
      fit_t1 <- f_fit_ti(df_t1, model_t1, nreps, showplot)
    }

    # gap time and number of skipped visits between drug dispensing visits
    df_ti <- vf1 %>%
      dplyr::mutate(
        time = dplyr::lead(.data$day) - .data$day,
        skipped = pmax(floor((.data$time - .data$target_days/2)/
                               .data$target_days), 0),
        k1 = .data$skipped + 1) %>%
      dplyr::filter(.data$row_id < dplyr::n())

    fit_ki <- f_fit_ki(df_ti, model_ki, nreps, showplot)
    fit_ti <- f_fit_ti(df_ti, model_ti, nreps, showplot)
  } else {
    # fit separate models for k0, t0, t1, ki, ti across drugs
    common_time_model = FALSE
    fit_k0 <- fit_t0 <- fit_t1 <- fit_ki <- fit_ti <- list()

    for (h in 1:l) {
      # observed dosing data for the drug under consideration
      vf1 <- vf %>% dplyr::filter(.data$drug == h)

      # time from randomization to the first drug dispensing visit
      df_k0 <- vf1 %>%
        dplyr::filter(.data$row_id == 1) %>%
        dplyr::mutate(
          time = .data$day,
          skipped = floor((.data$time - .data$target_days/2)/
                            .data$target_days) + 1)

      fit_k0[[h]] <- f_fit_ki(df_k0, model_k0, nreps, showplot)

      # no skipping
      df_t0 <- df_k0 %>%
        dplyr::filter(.data$skipped == 0) %>%
        dplyr::mutate(left = .data$time - 1, right = .data$time)

      fit_t0[[h]] <- f_fit_t0(df_t0, model_t0, nreps, showplot)

      # skipping
      df_t1 <- df_k0 %>%
        dplyr::filter(.data$skipped > 0) %>%
        dplyr::mutate(k1 = .data$skipped)

      if (nrow(df_t1) == 0) {
        if (tolower(model_t1) == "least squares") {
          model_t1_x = "Least squares"
        } else if (tolower(model_t1) == "least absolute deviations") {
          model_t1_x = "Least absolute deviations"
        }

        fit_t1[[h]] <- list(fit = list(model = model_t1_x,
                                       beta = 0,
                                       vbeta = 0,
                                       sigma = 0,
                                       df = NA,
                                       aic = NA,
                                       bic = NA),
                            fit_plot = NA,
                            theta = matrix(0, nreps, 2))
      } else {
        fit_t1[[h]] <- f_fit_ti(df_t1, model_t1, nreps, showplot)
      }

      # gap time and number of skipped visits between drug dispensing visits
      df_ti <- vf1 %>%
        dplyr::mutate(
          time = dplyr::lead(.data$day) - .data$day,
          skipped = pmax(floor((.data$time - .data$target_days/2)/
                                 .data$target_days), 0),
          k1 = .data$skipped + 1) %>%
        dplyr::filter(.data$row_id < dplyr::n())

      fit_ki[[h]] <- f_fit_ki(df_ti, model_ki, nreps, showplot)
      fit_ti[[h]] <- f_fit_ti(df_ti, model_ti, nreps, showplot)
    }
  }

  # fit separate models for di for different drugs
  fit_di <- purrr::map(1:l, function(h) {
    # observed dosing data for the drug under consideration
    vf1 <- vf %>% dplyr::filter(.data$drug == h)

    f_fit_di(vf1, model_di, nreps, showplot)
  })

  # output model fit results
  list(common_time_model = common_time_model,
       fit_k0 = fit_k0, fit_t0 = fit_t0, fit_t1 = fit_t1,
       fit_ki = fit_ki, fit_ti = fit_ti, fit_di = fit_di)
}

