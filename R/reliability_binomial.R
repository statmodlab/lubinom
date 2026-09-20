# -------------------------------------------------------------------------
# Threshold-based reliability functions for the binomial distribution
# -------------------------------------------------------------------------

#' Threshold-Based Reliability Functions for the Binomial Distribution
#'
#' Computes the discrete survival function, hazard function, and mean
#' residual count for a binomial random variable.
#'
#' For \eqn{X \sim Binomial(size,prob)}, define
#' \deqn{S(x)=P(X\geq x),}
#' \deqn{h(x)=P(X=x\mid X\geq x),}
#' and
#' \deqn{m(x)=E(X-x\mid X\geq x).}
#'
#' @param x Numeric vector of thresholds or support points.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param prob Binomial success probability in \eqn{[0,1]}.
#' @param log.p Logical; if `TRUE`, `sBIN()` returns logarithmic survival
#'   probabilities.
#'
#' @return
#' `sBIN()` returns the survival probability,
#' `hBIN()` returns the discrete hazard probability, and
#' `mrlBIN()` returns the mean residual count.
#'
#' @examples
#' x <- 0:5
#' sBIN(x, size = 5, prob = 0.3)
#' hBIN(x, size = 5, prob = 0.3)
#' mrlBIN(x, size = 5, prob = 0.3)
#'
#' @name BIN-reliability
NULL


#' @rdname BIN-reliability
#' @export
sBIN <- function(
    x,
    size,
    prob,
    log.p = FALSE) {

  .check_BIN(size, prob)
  .check_logical_scalar(log.p, "log.p")

  x <- as.numeric(x)

  vapply(
    x,
    function(xx) {

      if (is.na(xx)) {
        return(NA_real_)
      }

      if (xx == -Inf || xx <= 0) {
        return(
          if (log.p) 0 else 1
        )
      }

      if (xx == Inf || xx > size) {
        return(
          if (log.p) -Inf else 0
        )
      }

      if (!is.finite(xx)) {
        return(NA_real_)
      }

      stats::pbinom(
        q = ceiling(xx) - 1L,
        size = size,
        prob = prob,
        lower.tail = FALSE,
        log.p = log.p
      )
    },
    numeric(1)
  )
}


#' @rdname BIN-reliability
#' @export
hBIN <- function(
    x,
    size,
    prob) {

  .check_BIN(size, prob)

  x <- as.numeric(x)

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

      sx <- sBIN(
        x = xx,
        size = size,
        prob = prob
      )

      if (!is.finite(sx) ||
          sx <= 0) {
        return(NA_real_)
      }

      stats::dbinom(
        x = xx,
        size = size,
        prob = prob
      ) / sx
    },
    numeric(1)
  )
}


#' @rdname BIN-reliability
#' @export
mrlBIN <- function(
    x,
    size,
    prob) {

  .check_BIN(size, prob)

  x <- as.numeric(x)

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

      sx <- sBIN(
        x = xx,
        size = size,
        prob = prob
      )

      if (!is.finite(sx) ||
          sx <= 0) {
        return(NA_real_)
      }

      if (xx == size) {
        return(0)
      }

      jj <- (xx + 1L):size

      sum(
        sBIN(
          x = jj,
          size = size,
          prob = prob
        )
      ) / sx
    },
    numeric(1)
  )
}
