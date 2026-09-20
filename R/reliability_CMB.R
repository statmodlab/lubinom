# -------------------------------------------------------------------------
# Threshold-based reliability functions for the
# Conway-Maxwell-binomial distribution
# -------------------------------------------------------------------------


#' Threshold-Based Reliability Functions for the
#' Conway-Maxwell-Binomial Distribution
#'
#' Computes the discrete survival function, hazard function, and mean
#' residual count for a Conway-Maxwell-binomial random variable.
#'
#' For a Conway-Maxwell-binomial random variable \eqn{X}, define
#' \deqn{
#' S(x)=P(X\geq x),
#' }
#' \deqn{
#' h(x)=P(X=x\mid X\geq x),
#' }
#' and
#' \deqn{
#' m(x)=E(X-x\mid X\geq x).
#' }
#'
#' @param x Numeric vector of thresholds or support points.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param prob Probability parameter satisfying \eqn{0<prob<1}.
#' @param nu Positive dispersion parameter.
#' @param log.p Logical; if `TRUE`, `sCMB()` returns logarithmic survival
#'   probabilities.
#'
#' @return
#' `sCMB()` returns the survival probability,
#' `hCMB()` returns the discrete hazard probability, and
#' `mrlCMB()` returns the mean residual count.
#'
#' @examples
#' x <- 0:5
#' sCMB(x, size = 5, prob = 0.3, nu = 2)
#' hCMB(x, size = 5, prob = 0.3, nu = 2)
#' mrlCMB(x, size = 5, prob = 0.3, nu = 2)
#'
#' @name CMB-reliability
NULL


# -------------------------------------------------------------------------
# Survival function
# -------------------------------------------------------------------------

#' @rdname CMB-reliability
#' @export
sCMB <- function(
    x,
    size,
    prob,
    nu,
    log.p = FALSE) {

  .check_CMB(
    size = size,
    prob = prob,
    nu = nu
  )

  .check_logical_scalar(
    log.p,
    "log.p"
  )

  x <- as.numeric(x)

  probs <- .CMB_probs(
    size = size,
    prob = prob,
    nu = nu
  )

  if (is.null(probs)) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  surv <- rev(
    cumsum(
      rev(probs)
    )
  )

  surv[1L] <- 1

  vapply(
    x,
    function(xx) {

      if (is.na(xx)) {
        return(NA_real_)
      }

      if (xx == -Inf ||
          xx <= 0) {

        value <- 1

      } else if (
        xx == Inf ||
        xx > size
      ) {

        value <- 0

      } else if (
        !is.finite(xx)
      ) {

        return(NA_real_)

      } else {

        value <-
          surv[
            ceiling(xx) + 1L
          ]
      }

      if (log.p) {
        log(value)
      } else {
        value
      }
    },
    numeric(1)
  )
}


# -------------------------------------------------------------------------
# Discrete hazard function
# -------------------------------------------------------------------------

#' @rdname CMB-reliability
#' @export
hCMB <- function(
    x,
    size,
    prob,
    nu) {

  .check_CMB(
    size = size,
    prob = prob,
    nu = nu
  )

  x <- as.numeric(x)

  probs <- .CMB_probs(
    size = size,
    prob = prob,
    nu = nu
  )

  if (is.null(probs)) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  surv <- rev(
    cumsum(
      rev(probs)
    )
  )

  surv[1L] <- 1

  vapply(
    x,
    function(xx) {

      if (is.na(xx) ||
          !is.finite(xx) ||
          xx < 0 ||
          xx > size ||
          xx != floor(xx)) {
        return(NA_real_)
      }

      index <-
        as.integer(xx) + 1L

      sx <- surv[index]

      if (!is.finite(sx) ||
          sx <= 0) {
        return(NA_real_)
      }

      value <-
        probs[index] /
        sx

      if (xx == size) {
        value <- 1
      }

      value
    },
    numeric(1)
  )
}


# -------------------------------------------------------------------------
# Mean residual count
# -------------------------------------------------------------------------

#' @rdname CMB-reliability
#' @export
mrlCMB <- function(
    x,
    size,
    prob,
    nu) {

  .check_CMB(
    size = size,
    prob = prob,
    nu = nu
  )

  x <- as.numeric(x)

  probs <- .CMB_probs(
    size = size,
    prob = prob,
    nu = nu
  )

  if (is.null(probs)) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  support <- 0:size

  surv <- rev(
    cumsum(
      rev(probs)
    )
  )

  surv[1L] <- 1

  vapply(
    x,
    function(xx) {

      if (is.na(xx) ||
          !is.finite(xx) ||
          xx < 0 ||
          xx > size ||
          xx != floor(xx)) {
        return(NA_real_)
      }

      if (xx == size) {
        return(0)
      }

      index <-
        as.integer(xx) + 1L

      sx <- surv[index]

      if (!is.finite(sx) ||
          sx <= 0) {
        return(NA_real_)
      }

      sum(
        (
          support[
            (index + 1L):(size + 1L)
          ] -
            xx
        ) *
          probs[
            (index + 1L):(size + 1L)
          ]
      ) / sx
    },
    numeric(1)
  )
}
