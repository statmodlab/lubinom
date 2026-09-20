# -------------------------------------------------------------------------
# Threshold-based reliability functions for the beta-binomial distribution
# -------------------------------------------------------------------------


#' Threshold-Based Reliability Functions for the Beta-Binomial Distribution
#'
#' Computes the discrete survival function, hazard function, and mean
#' residual count for a beta-binomial random variable.
#'
#' For a beta-binomial random variable \eqn{X}, define
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
#' @param alpha Positive first beta shape parameter.
#' @param beta Positive second beta shape parameter.
#' @param log.p Logical; if `TRUE`, `sBB()` returns logarithmic survival
#'   probabilities.
#'
#' @return
#' `sBB()` returns the survival probability,
#' `hBB()` returns the discrete hazard probability, and
#' `mrlBB()` returns the mean residual count.
#'
#' @examples
#' x <- 0:5
#' sBB(x, size = 5, alpha = 1, beta = 3)
#' hBB(x, size = 5, alpha = 1, beta = 3)
#' mrlBB(x, size = 5, alpha = 1, beta = 3)
#'
#' @name BB-reliability
NULL


# -------------------------------------------------------------------------
# Survival function
# -------------------------------------------------------------------------

#' @rdname BB-reliability
#' @export
sBB <- function(
    x,
    size,
    alpha,
    beta,
    log.p = FALSE) {

  .check_BB_parameters(
    size = size,
    alpha = alpha,
    beta = beta
  )

  .check_logical_scalar(
    log.p,
    "log.p"
  )

  x <- as.numeric(x)

  probs <- .BB_probs(
    size = size,
    alpha = alpha,
    beta = beta
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

#' @rdname BB-reliability
#' @export
hBB <- function(
    x,
    size,
    alpha,
    beta) {

  .check_BB_parameters(
    size = size,
    alpha = alpha,
    beta = beta
  )

  x <- as.numeric(x)

  probs <- .BB_probs(
    size = size,
    alpha = alpha,
    beta = beta
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

#' @rdname BB-reliability
#' @export
mrlBB <- function(
    x,
    size,
    alpha,
    beta) {

  .check_BB_parameters(
    size = size,
    alpha = alpha,
    beta = beta
  )

  x <- as.numeric(x)

  probs <- .BB_probs(
    size = size,
    alpha = alpha,
    beta = beta
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
