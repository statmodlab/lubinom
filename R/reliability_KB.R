# -------------------------------------------------------------------------
# Threshold-based reliability functions for the
# Kumaraswamy-binomial distribution
# -------------------------------------------------------------------------


#' Threshold-Based Reliability Functions for the
#' Kumaraswamy-Binomial Distribution
#'
#' Computes the discrete survival function, hazard function, and mean
#' residual count for a Kumaraswamy-binomial random variable.
#'
#' For a Kumaraswamy-binomial random variable \eqn{X}, define
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
#' @param a Positive first Kumaraswamy shape parameter.
#' @param b Positive second Kumaraswamy shape parameter.
#' @param log.p Logical; if `TRUE`, `sKB()` returns logarithmic survival
#'   probabilities.
#'
#' @return
#' `sKB()` returns the survival probability,
#' `hKB()` returns the discrete hazard probability, and
#' `mrlKB()` returns the mean residual count.
#'
#' @examples
#' x <- 0:5
#' sKB(x, size = 5, a = 1, b = 3)
#' hKB(x, size = 5, a = 1, b = 3)
#' mrlKB(x, size = 5, a = 1, b = 3)
#'
#' @name KB-reliability
NULL


# -------------------------------------------------------------------------
# Survival function
# -------------------------------------------------------------------------

#' @rdname KB-reliability
#' @export
sKB <- function(
    x,
    size,
    a,
    b,
    log.p = FALSE) {

  .check_KB(
    size = size,
    a = a,
    b = b
  )

  .check_logical_scalar(
    log.p,
    "log.p"
  )

  x <- as.numeric(x)

  probs <- .KB_probs(
    size = size,
    a = a,
    b = b
  )

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

#' @rdname KB-reliability
#' @export
hKB <- function(
    x,
    size,
    a,
    b) {

  .check_KB(
    size = size,
    a = a,
    b = b
  )

  x <- as.numeric(x)

  probs <- .KB_probs(
    size = size,
    a = a,
    b = b
  )

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

#' @rdname KB-reliability
#' @export
mrlKB <- function(
    x,
    size,
    a,
    b) {

  .check_KB(
    size = size,
    a = a,
    b = b
  )

  x <- as.numeric(x)

  probs <- .KB_probs(
    size = size,
    a = a,
    b = b
  )

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
