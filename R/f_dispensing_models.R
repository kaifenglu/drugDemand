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
#'     - \code{model}: The specific model used in the analysis.
#'
#'     - \code{theta}: The estimated model parameters.
#'
#'     - \code{vtheta}: The estimated covariance matrix of \code{theta}.
#'
#'     - \code{aic}: The Akaike Information Criterion value.
#'
#'     - \code{bic}: The Bayesian Information Criterion value.
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
#' observed <- f_dose_observed(df2, visitview2, showplot = FALSE)
#' vf <- observed$vf
#'
#' vf <- vf %>% left_join(dosing_schedule_df, by = "kit")
#'
#' # time from randomization to the first drug dispensing visit
#' df_k0 <- vf %>%
#'   filter(row_id == 1) %>%
#'   mutate(time = day,
#'          skipped = floor((time - target_days/2)/target_days) + 1)
#'
#' # no skipping
#' df_t0 <- df_k0 %>%
#'   filter(skipped == 0) %>%
#'   mutate(left = time - 1, right = time)
#'
#' t0_fit <- f_fit_t0(df_t0, model = "log-logistic", nreps = 200)
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

    f_nloglik <- function(theta, df, model) {
      if (tolower(model) == "exponential") {
        -sum(log(pexp(df$right, exp(theta)) -
                   pexp(df$left, exp(theta))))
      } else if (tolower(model) == "weibull") {
        -sum(log(pweibull(df$right, exp(-theta[2]), exp(theta[1])) -
                   pweibull(df$left, exp(-theta[2]), exp(theta[1]))))
      } else if (tolower(model) == "log-logistic") {
        -sum(log(plogis(log(df$right), theta[1], exp(theta[2])) -
                   plogis(log(df$left), theta[1], exp(theta[2]))))
      } else if (tolower(model) == "log-normal") {
        -sum(log(pnorm(log(df$right), theta[1], exp(theta[2])) -
                   pnorm(log(df$left), theta[1], exp(theta[2]))))
      }
    }

    if (model == "exponential") {
      a <- survreg(Surv(right) ~ 1, data = df, dist = "exponential")
      theta = -as.numeric(a$coefficients)

      opt = optim(theta, f_nloglik, gr = NULL, df, model,
                  method = "Brent", lower = theta - 2, upper = theta + 2,
                  hessian = TRUE)

      fit <- list(model = "Exponential",
                  theta = opt$par,
                  vtheta = as.numeric(solve(opt$hessian)),
                  aic = 2*opt$value + 2,
                  bic = 2*opt$value + log(n))

      rate = exp(fit$theta)
      p.fit = pexp(y+1, rate) - pexp(y, rate)

      post = rnorm(nreps, mean = fit$theta, sd = sqrt(fit$vtheta))
    } else if (model == "weibull") {
      a <- survreg(Surv(right) ~ 1, data = df, dist = "weibull")
      theta = c(as.numeric(a$coefficients), log(a$scale))

      opt = optim(theta, f_nloglik, gr = NULL, df, model, hessian = TRUE)

      fit <- list(model = "Weibull",
                  theta = opt$par,
                  vtheta = matrix(solve(opt$hessian), 2, 2),
                  aic = 2*opt$value + 4,
                  bic = 2*opt$value + 2*log(n))

      lambda = exp(fit$theta[1])
      k = exp(-fit$theta[2])
      p.fit = pweibull(y+1, k, lambda) - pweibull(y, k, lambda)

      post = mvtnorm::rmvnorm(nreps, mean = fit$theta, sigma = fit$vtheta)
    } else if (model == "log-logistic") {
      a <- survreg(Surv(right) ~ 1, data = df, dist = "loglogistic")
      theta = c(as.numeric(a$coefficients), log(a$scale))

      opt = optim(theta, f_nloglik, gr = NULL, df, model, hessian = TRUE)

      fit <- list(model = "Log-logistic",
                  theta = opt$par,
                  vtheta = matrix(solve(opt$hessian), 2, 2),
                  aic = 2*opt$value + 4,
                  bic = 2*opt$value + 2*log(n))

      mu = fit$theta[1]
      sigma = exp(fit$theta[2])
      p.fit = plogis(log(y+1), mu, sigma) - plogis(log(y), mu, sigma)

      post = mvtnorm::rmvnorm(nreps, mean = fit$theta, sigma = fit$vtheta)
    } else if (model == "log-normal") {
      a <- survreg(Surv(right) ~ 1, data = df, dist = "lognormal")
      theta = c(as.numeric(a$coefficients), log(a$scale))

      opt = optim(theta, f_nloglik, gr = NULL, df, model, hessian = TRUE)

      fit <- list(model = "Log-normal",
                  theta = opt$par,
                  vtheta = matrix(solve(opt$hessian), 2, 2),
                  aic = 2*opt$value + 4,
                  bic = 2*opt$value + 2*log(n))

      mu = fit$theta[1]
      sigma = exp(fit$theta[2])
      p.fit = plnorm(y+1, mu, sigma) - plnorm(y, mu, sigma)

      post = mvtnorm::rmvnorm(nreps, mean = fit$theta, sigma = fit$vtheta)
    } else {
      stop("incorrect model for T0")
    }
  }

  # graphically assess the model fit
  if (fit$model != "Constant") {
    modeltext = fit$model
    aictext = paste("AIC:", formatC(fit$aic, format = "f", digits = 2))
    bictext = paste("BIC:", formatC(fit$bic, format = "f", digits = 2))
  } else {
    modeltext = ""
    aictext = ""
    bictext = ""
  }

  gf <- tibble(y = y, p.obs = p.obs, p.fit = p.fit)

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
#'     - \code{model}: The specific model used in the analysis.
#'
#'     - \code{theta}: The estimated model parameters.
#'
#'     - \code{vtheta}: The estimated covariance matrix of \code{theta}.
#'
#'     - \code{aic}: The Akaike Information Criterion value.
#'
#'     - \code{bic}: The Bayesian Information Criterion value.
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
#' observed <- f_dose_observed(df2, visitview2, showplot = FALSE)
#' vf <- observed$vf
#'
#' vf <- vf %>% left_join(dosing_schedule_df, by = "kit")
#'
#' df_ti <- vf %>%
#'   mutate(time = lead(day) - day,
#'          skipped = pmax(floor((time - target_days/2)/target_days), 0),
#'          k1 = skipped + 1) %>%
#'   filter(row_id < n())
#'
#' ki_fit <- f_fit_ki(df_ti, model = "zero-inflated poisson", nreps = 200)
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
                  theta = as.numeric(a$coefficients),
                  vtheta = as.numeric(vcov(a)),
                  aic = a$aic,
                  bic = a$aic - 2 + log(n))

      rate = exp(fit$theta)
      p.fit = dpois(y, rate)

      post = rnorm(nreps, mean = fit$theta, sd = sqrt(fit$vtheta))
    } else if (model == "zero-inflated poisson") {
      f_nloglik <- function(theta, df) {
        lambda = exp(theta[1])
        pi = plogis(theta[2])
        -sum(log(pi*(df$skipped == 0) + (1-pi)*dpois(df$skipped, lambda)))
      }

      theta = c(log(sum(df$skipped)/n), 0)
      opt <- optim(theta, f_nloglik, gr = NULL, df, hessian = TRUE)

      fit <- list(model = "Zero-inflated Poisson",
                  theta = opt$par,
                  vtheta = matrix(solve(opt$hessian), 2, 2),
                  aic = 2*opt$value + 4,
                  bic = 2*opt$value + 2*log(n))

      lambda = exp(fit$theta[1])
      pi = plogis(fit$theta[2])

      p.fit <- pi*(y == 0) + (1-pi)*dpois(y, lambda)

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
  if (fit$model != "Constant") {
    modeltext = fit$model
    aictext = paste("AIC:", formatC(fit$aic, format = "f", digits = 2))
    bictext = paste("BIC:", formatC(fit$bic, format = "f", digits = 2))
  } else {
    modeltext = ""
    aictext = ""
    bictext = ""
  }

  gf <- tibble(y = y, p.obs = p.obs, p.fit = p.fit)

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
#'     - \code{model}: The specific model used in the analysis.
#'
#'     - \code{beta}: The estimated regression coefficient for the covariate.
#'
#'     - \code{vbeta}: The estimated variance of \code{beta}.
#'
#'     - \code{sigma}: The estimated residual standard deviation.
#'
#'     - \code{df}: The residual degrees-of-freedom.
#'
#'     - \code{aic}: The Akaike Information Criterion value.
#'
#'     - \code{bic}: The Bayesian Information Criterion value.
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
#' observed <- f_dose_observed(df2, visitview2, showplot = FALSE)
#' vf <- observed$vf
#'
#' vf <- vf %>% left_join(dosing_schedule_df, by = "kit")
#'
#' df_ti <- vf %>%
#'   mutate(time = lead(day) - day,
#'          skipped = pmax(floor((time - target_days/2)/target_days), 0),
#'          k1 = skipped + 1) %>%
#'   filter(row_id < n())
#'
#' ti_fit <- f_fit_ti(df_ti, model = "least squares", nreps = 200)
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

  # number of unique combinations of response and covariates
  df_unique <- df %>%
    group_by(.data$time, .data$k1) %>%
    slice(1)

  if (nrow(df_unique) == 1) {
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

    mu = fit$beta*df$k1
    b = fit$sigma

    p.fit <- purrr::map_vec(y, function(y) {
      mean(pnorm(y + 0.5, mu, b) - pnorm(y - 0.5, mu, b))
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

    plaplace <- function(x, mu, b) {
      0.5*exp((x-mu)/b)*(x <= mu) + (1 - 0.5*exp(-(x-mu)/b))*(x > mu)
    }

    mu = fit$beta*df$k1
    b = fit$sigma

    p.fit <- purrr::map_vec(y, function(y) {
      mean(plaplace(y + 0.5, mu, b) - plaplace(y - 0.5, mu, b))
    })
  }

  gf = tibble(y = y, p.obs = p.obs, p.fit = p.fit)

  if (fit$model != "Constant") {
    modeltext = fit$model
    aictext = paste("AIC:", formatC(fit$aic, format = "f", digits = 2))
    bictext = paste("BIC:", formatC(fit$bic, format = "f", digits = 2))
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
#'   \code{day}, \code{kit}, and \code{dose}.
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
#'     - \code{model}: The specific model used in the analysis.
#'
#'     - \code{mud}: The estimated mean dose.
#'
#'     - \code{vmud}: The estimated variance of \code{mud}.
#'
#'     - \code{sigmab}: The estimated between-subject standard deviation.
#'
#'     - \code{sigmae}: The estimated within-subject residual standard
#'       deviation.
#'
#'     - \code{aic}: The Akaike Information Criterion value.
#'
#'     - \code{bic}: The Bayesian Information Criterion value.
#'
#' * \code{fit_plot}: A fitted dose bar chart.
#'
#' * \code{theta}: Posterior draws of model parameters.
#'
#'     - \code{fixed}: Posterior draws of fixed model parameters:
#'       \code{mud}, \code{sigmab}, and \code{sigmae}.
#'
#'     - \code{random}: Posterior draws of subject random effects.
#'
#'     - \code{usubjid}: The unique subject ID associated with
#'       the subject random effects.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' library(dplyr)
#'
#' observed <- f_dose_observed(df2, visitview2, showplot = FALSE)
#' vf <- observed$vf
#'
#' vf1 <- vf %>% filter(kit == 3)
#' di_fit <- f_fit_di(vf1, model = "linear mixed-effects model", nreps = 200)
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

    gf = tibble(y = y, p.obs = p.obs, p.fit = p.fit)

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

      gf = tibble(y = y, p.obs = p.obs, p.fit = p.fit)

      # draw sigmae and then mud given sigmae from posterior
      b2 = sqrt(a$df.residual * fit$sigmae^2 / rchisq(nreps, a$df.residual))
      b1 = fit$mud + rnorm(nreps) * sqrt(fit$vmud) / fit$sigmae * b2
      theta_fix = matrix(c(b1, rep(0, nreps), b2), nreps, 3)
      theta_ran = matrix(0, nreps, N)
    } else if (model == "linear mixed-effects model") {
      a <- nlme::lme(dose ~ 1, random = ~ 1 | usubjid, data = df)

      vc <- matrix(as.numeric(nlme::VarCorr(a)), 2, 2)

      fit <- list(model = "Linear mixed-effects model",
                  mud = as.numeric(a$coefficients$fixed),
                  vmud = as.numeric(vcov(a)),
                  sigmab = vc[1,2],
                  sigmae = vc[2,2],
                  aic = as.numeric(AIC(a)),
                  bic = as.numeric(BIC(a)))

      # mean incorporating subject random effects
      mu = as.numeric(a$fitted[,"usubjid"])
      sigmae = fit$sigmae

      p.fit <- purrr::map_vec(y, function(y) {
        mean(pnorm(y + 0.5, mu, sigmae) - pnorm(y - 0.5, mu, sigmae))
      })

      gf = tibble(y = y, p.obs = p.obs, p.fit = p.fit)

      # number of observations and mean dose by subject
      df1 <- df %>%
        summarise(n = n(), d = mean(.data$dose), .groups = "drop_last")

      # fix the variance parameters to avoid drawing extreme values
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
  if (fit$model != "Constant") {
    modeltext = fit$model
    aictext = paste("AIC:", formatC(fit$aic, format = "f", digits = 2))
    bictext = paste("BIC:", formatC(fit$bic, format = "f", digits = 2))
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
#'   \code{drug}, \code{drug_name}, \code{kit}, \code{kit_name},
#'   \code{usubjid}, \code{treatment}, \code{treatment_description},
#'   \code{arrivalTime}, \code{time}, \code{event}, \code{dropout},
#'   \code{day}, \code{dose}, \code{cum_dose}, and \code{row_id}.
#' @param dosing_schedule_df A data frame providing dosing schedule
#'   information. It contains the following variables:
#'   \code{kit}, \code{target_days}, \code{target_dose}, and
#'   \code{max_cycles}.
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
#' * \code{k0_fit}: The model fit for the number of skipped
#'   visits between randomization and the first drug dispensing visit.
#'
#' * \code{t0_fit}: The model fit for the gap time between
#'   randomization and the first drug dispensing visit when there is
#'   no visit skipping.
#'
#' * \code{t1_fit}: The model fit for the gap time between
#'   randomization and the first drug dispensing visit when there is
#'   visit skipping.
#'
#' * \code{ki_fit}: The model fit for the number of skipped
#'   visits between two consecutive drug dispensing visits.
#'
#' * \code{ti_fit}: The model fit for the gap time between two
#'   consecutive drug dispensing visits.
#'
#' * \code{di_fit}: The model fit for the dispensed doses at drug
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
#' observed <- f_dose_observed(df2, visitview2, showplot = FALSE)
#'
#' dispensing_models <- f_dispensing_models(
#'   observed$vf, dosing_schedule_df,
#'   model_k0 = "zero-inflated poisson",
#'   model_t0 = "log-logistic",
#'   model_t1 = "least squares",
#'   model_ki = "zero-inflated poisson",
#'   model_ti = "least squares",
#'   model_di = "linear mixed-effects model",
#'   nreps = 200, showplot = FALSE)
#'
#' dispensing_models$ki_fit$fit_plot
#
#' @export
f_dispensing_models <- function(
    vf, dosing_schedule_df, model_k0, model_t0, model_t1,
    model_ki, model_ti, model_di, nreps, showplot = TRUE) {

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

  vf = vf %>% left_join(dosing_schedule_df, by = "kit")
  l = nrow(dosing_schedule_df)

  vf0 <- vf %>%
    group_by(.data$kit, .data$kit_name) %>%
    slice(1)

  kit_name = vf0$kit_name

  f_kit_name <- function(fit_plot, kit_name, h) {
    fit_plot %>%
      plotly::layout(
        annotations = list(
          x = 0.5, y = 1,
          text = paste0("<b>", kit_name[h], "</b>"),
          xanchor = "center", yanchor = "middle",
          showarrow = FALSE, xref='paper', yref='paper'))
  }

  if (length(unique(dosing_schedule_df$target_days)) == 1) {
    # fit a common model for k0, t0, t1, ki, and ti across kit types
    common_time_model = TRUE

    # time from randomization to the first drug dispensing visit
    df_k0 <- vf %>%
      filter(.data$row_id == 1) %>%
      mutate(time = .data$day,
             skipped = floor((.data$time - .data$target_days/2)/
                               .data$target_days) + 1)

    k0_fit <- f_fit_ki(df_k0, model_k0, nreps, showplot)

    # no skipping
    df_t0 <- df_k0 %>%
      filter(.data$skipped == 0) %>%
      mutate(left = .data$time - 1, right = .data$time)

    t0_fit <- f_fit_t0(df_t0, model_t0, nreps, showplot)

    # skipping
    df_t1 <- df_k0 %>%
      filter(.data$skipped > 0) %>%
      mutate(k1 = .data$skipped)

    if (nrow(df_t1) == 0) {
      model_t1_x = stringr::str_to_title(model_t1)

      t1_fit <- list(fit = list(model = model_t1_x,
                                beta = 0,
                                vbeta = 0,
                                sigma = 0,
                                df = NA,
                                aic = NA,
                                bic = NA),
                     fit_plot = NULL,
                     theta = matrix(0, nreps, 2))
    } else {
      t1_fit <- f_fit_ti(df_t1, model_t1, nreps, showplot)
    }

    # gap time and number of skipped visits between drug dispensing visits
    df_ti <- vf %>%
      mutate(time = lead(.data$day) - .data$day,
             skipped = pmax(floor((.data$time - .data$target_days/2)/
                                    .data$target_days), 0),
             k1 = .data$skipped + 1) %>%
      filter(.data$row_id < n())

    ki_fit <- f_fit_ki(df_ti, model_ki, nreps, showplot)
    ti_fit <- f_fit_ti(df_ti, model_ti, nreps, showplot)
  } else {
    # fit separate models for k0, t0, t1, ki, ti across kit types
    common_time_model = FALSE
    k0_fit <- t0_fit <- t1_fit <- ki_fit <- ti_fit <- list()

    for (h in 1:l) {
      # observed dosing data for the kit type under consideration
      vf1 <- vf %>% filter(.data$kit == h)

      # time from randomization to the first drug dispensing visit
      df_k0 <- vf1 %>%
        filter(.data$row_id == 1) %>%
        mutate(time = .data$day,
               skipped = floor((.data$time - .data$target_days/2)/
                                 .data$target_days) + 1)

      k0_fit[[h]] <- f_fit_ki(df_k0, model_k0, nreps, showplot)
      k0_fit[[h]]$fit_plot <- f_kit_name(k0_fit[[h]]$fit_plot, kit_name, h)

      # no skipping
      df_t0 <- df_k0 %>%
        filter(.data$skipped == 0) %>%
        mutate(left = .data$time - 1, right = .data$time)

      t0_fit[[h]] <- f_fit_t0(df_t0, model_t0, nreps, showplot)
      t0_fit[[h]]$fit_plot <- f_kit_name(t0_fit[[h]]$fit_plot, kit_name, h)

      # skipping
      df_t1 <- df_k0 %>%
        filter(.data$skipped > 0) %>%
        mutate(k1 = .data$skipped)

      if (nrow(df_t1) == 0) {
        model_t1_x = stringr::str_to_title(model_t1)

        t1_fit[[h]] <- list(fit = list(model = model_t1_x,
                                       beta = 0,
                                       vbeta = 0,
                                       sigma = 0,
                                       df = NA,
                                       aic = NA,
                                       bic = NA),
                            fit_plot = NULL,
                            theta = matrix(0, nreps, 2))
      } else {
        t1_fit[[h]] <- f_fit_ti(df_t1, model_t1, nreps, showplot)
        t1_fit[[h]]$fit_plot <- f_kit_name(t1_fit[[h]]$fit_plot, kit_name, h)
      }

      # gap time and number of skipped visits between drug dispensing visits
      df_ti <- vf1 %>%
        mutate(time = lead(.data$day) - .data$day,
               skipped = pmax(floor((.data$time - .data$target_days/2)/
                                      .data$target_days), 0),
               k1 = .data$skipped + 1) %>%
        filter(.data$row_id < n())

      ki_fit[[h]] <- f_fit_ki(df_ti, model_ki, nreps, showplot)
      ki_fit[[h]]$fit_plot <- f_kit_name(ki_fit[[h]]$fit_plot, kit_name, h)

      ti_fit[[h]] <- f_fit_ti(df_ti, model_ti, nreps, showplot)
      ti_fit[[h]]$fit_plot <- f_kit_name(ti_fit[[h]]$fit_plot, kit_name, h)
    }
  }

  # fit separate models for di for different kit types
  di_fit <- list()

  for (h in 1:l) {
    # observed dosing data for the kit type under consideration
    vf1 <- vf %>% filter(.data$kit == h)

    di_fit[[h]] <- f_fit_di(vf1, model_di, nreps, showplot)
    di_fit[[h]]$fit_plot <- f_kit_name(di_fit[[h]]$fit_plot, kit_name, h)
  }

  # output model fit results
  list(common_time_model = common_time_model,
       k0_fit = k0_fit, t0_fit = t0_fit, t1_fit = t1_fit,
       ki_fit = ki_fit, ti_fit = ti_fit, di_fit = di_fit)
}

