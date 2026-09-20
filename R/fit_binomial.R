# -------------------------------------------------------------------------
# Maximum-likelihood estimation for the binomial distribution
# -------------------------------------------------------------------------


.check_BIN_sample <- function(x, size) {

  if (length(size) != 1L ||
      !is.numeric(size) ||
      !is.finite(size) ||
      size < 1 ||
      size != floor(size)) {
    stop(
      "'size' must be a strictly positive integer.",
      call. = FALSE
    )
  }

  if (!is.numeric(x) ||
      length(x) < 1L) {
    stop(
      "'x' must be a non-empty numeric vector.",
      call. = FALSE
    )
  }

  if (any(!is.finite(x))) {
    stop(
      "All observations in 'x' must be finite.",
      call. = FALSE
    )
  }

  if (any(
    x < 0 |
      x > size |
      x != floor(x)
  )) {
    stop(
      "All observations in 'x' must belong to {0, ..., size}.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Maximum-Likelihood Estimation for the Binomial Distribution
#'
#' Fits a binomial model to bounded count data with known group size.
#'
#' @param x Observed bounded counts.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param conf.level Confidence level for the Wald interval.
#'
#' @details
#' For independent observations
#' \eqn{X_i\sim Binomial(size,p)}, the maximum-likelihood estimator is
#' \deqn{\widehat p = \frac{\bar X}{size}.}
#'
#' For an interior estimate, the standard error is obtained from the
#' binomial Fisher information. A Wald confidence interval truncated to
#' \eqn{[0,1]} is reported.
#'
#' @return
#' A list containing parameter estimates, standard errors, confidence
#' limits, log-likelihood, information criteria, and model metadata.
#'
#' @examples
#' x <- c(0, 1, 0, 2, 1, 0)
#' fitBIN_mle(x, size = 4)
#'
#' @export
fitBIN_mle <- function(
    x,
    size,
    conf.level = 0.95) {

  x <- as.numeric(x)

  .check_BIN_sample(
    x = x,
    size = size
  )

  if (length(conf.level) != 1L ||
      !is.numeric(conf.level) ||
      !is.finite(conf.level) ||
      conf.level <= 0 ||
      conf.level >= 1) {
    stop(
      "'conf.level' must be strictly between 0 and 1.",
      call. = FALSE
    )
  }

  N <- length(x)
  k <- 1L

  prob <- mean(x) / size

  logLik <- sum(
    stats::dbinom(
      x = x,
      size = size,
      prob = prob,
      log = TRUE
    )
  )

  se <- if (
    prob > 0 &&
    prob < 1
  ) {
    sqrt(
      prob *
        (1 - prob) /
        (N * size)
    )
  } else {
    NA_real_
  }

  zcrit <- stats::qnorm(
    1 - (1 - conf.level) / 2
  )

  lower <- if (
    is.finite(se)
  ) {
    max(
      0,
      prob - zcrit * se
    )
  } else {
    NA_real_
  }

  upper <- if (
    is.finite(se)
  ) {
    min(
      1,
      prob + zcrit * se
    )
  } else {
    NA_real_
  }

  AIC <-
    -2 * logLik +
    2 * k

  AICc <- if (
    N > k + 1L
  ) {
    AIC +
      2 * k * (k + 1) /
      (N - k - 1)
  } else {
    NA_real_
  }

  BIC <-
    -2 * logLik +
    k * log(N)

  out <- list(
    call = match.call(),

    par = c(
      prob = prob
    ),

    se = c(
      prob = se
    ),

    coefficients = matrix(
      c(
        prob,
        se,
        lower,
        upper
      ),
      nrow = 1L,
      dimnames = list(
        "prob",
        c(
          "Estimate",
          "SE",
          "Lower",
          "Upper"
        )
      )
    ),

    logLik = logLik,
    AIC = AIC,
    AICc = AICc,
    BIC = BIC,

    convergence = 0L,
    success = TRUE,

    size = size,
    n = N,
    nobs = N,
    npar = k,

    conf.level = conf.level,
    data = x
  )

  out
}
