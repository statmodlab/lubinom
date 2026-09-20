# -------------------------------------------------------------------------
# Reliability functions for the Lambert-uniform binomial distribution
# -------------------------------------------------------------------------

#' Threshold-Based Reliability Functions for the LUB Distribution
#'
#' Computes the discrete survival function, hazard function, and mean
#' residual count for a Lambert-uniform binomial random variable.
#'
#' Let \eqn{X \sim LUB(size,\beta)} with support
#' \eqn{\{0,1,\ldots,size\}} and \eqn{\beta>0}. The functions are defined as
#' \deqn{S(x)=P(X\geq x),}
#' \deqn{h(x)=P(X=x\mid X\geq x)=\frac{P(X=x)}{S(x)},}
#' and
#' \deqn{m(x)=E(X-x\mid X\geq x).}
#'
#' @param x Numeric vector of thresholds or support points at which the
#'   function is evaluated.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param beta Positive Lambert-uniform parameter.
#' @param log.p Logical; if `TRUE`, `sLUB()` returns logarithmic survival
#'   probabilities.
#'
#' @details
#' The survival convention is \eqn{S(x)=P(X\geq x)}. Consequently,
#' \eqn{h(size)=1} and \eqn{m(size)=0}.
#'
#' For integer \eqn{x \in \{0,\ldots,size\}}, the mean residual count can
#' also be written as
#' \deqn{
#' m(x)=\frac{\sum_{j=x+1}^{size}S(j)}{S(x)}.
#' }
#'
#' When \eqn{\beta=1}, the LUB distribution reduces to the discrete uniform
#' distribution on \eqn{\{0,\ldots,size\}}, giving
#' \deqn{
#' S(x)=\frac{size-x+1}{size+1},\qquad
#' h(x)=\frac{1}{size-x+1},\qquad
#' m(x)=\frac{size-x}{2}.
#' }
#'
#' @return
#' `sLUB()` returns the discrete survival probability,
#' `hLUB()` returns the discrete hazard probability, and
#' `mrlLUB()` returns the mean residual count.
#'
#' @examples
#' x <- 0:5
#' sLUB(x, size = 5, beta = 2)
#' hLUB(x, size = 5, beta = 2)
#' mrlLUB(x, size = 5, beta = 2)
#'
#' @name LUB-reliability
NULL


# -------------------------------------------------------------------------
# Survival function
# S(x) = P(X >= x)
# -------------------------------------------------------------------------

#' @rdname LUB-reliability
#' @export
sLUB <- function(
    x,
    size,
    beta,
    log.p = FALSE) {

  .check_LUB(size, beta)
  .check_logical_scalar(log.p, "log.p")

  x <- as.numeric(x)

  probs <- .LUB_probs(
    size = size,
    beta = beta
  )

  if (is.null(probs)) {
    return(rep(NA_real_, length(x)))
  }

  # S(k) = P(X >= k), k = 0, ..., size
  surv <- rev(cumsum(rev(probs)))
  surv[1L] <- 1

  out <- vapply(
    x,
    function(xx) {

      if (is.na(xx)) {
        return(NA_real_)
      }

      if (xx == -Inf || xx <= 0) {
        value <- 1
      } else if (xx == Inf || xx > size) {
        value <- 0
      } else if (!is.finite(xx)) {
        return(NA_real_)
      } else {
        value <- surv[ceiling(xx) + 1L]
      }

      if (log.p) {
        log(value)
      } else {
        value
      }
    },
    numeric(1)
  )

  out
}


# -------------------------------------------------------------------------
# Discrete hazard function
# h(x) = P(X = x | X >= x)
# -------------------------------------------------------------------------

#' @rdname LUB-reliability
#' @export
hLUB <- function(
    x,
    size,
    beta) {

  .check_LUB(size, beta)

  x <- as.numeric(x)

  probs <- .LUB_probs(
    size = size,
    beta = beta
  )

  if (is.null(probs)) {
    return(rep(NA_real_, length(x)))
  }

  surv <- rev(cumsum(rev(probs)))
  surv[1L] <- 1

  out <- vapply(
    x,
    function(xx) {

      if (is.na(xx) ||
          !is.finite(xx) ||
          xx < 0 ||
          xx > size ||
          xx != floor(xx)) {
        return(NA_real_)
      }

      index <- as.integer(xx) + 1L
      sx <- surv[index]

      if (!is.finite(sx) || sx <= 0) {
        return(NA_real_)
      }

      value <- probs[index] / sx

      # Exact theoretical endpoint
      if (xx == size) {
        value <- 1
      }

      value
    },
    numeric(1)
  )

  out
}


# -------------------------------------------------------------------------
# Mean residual count
# m(x) = E(X - x | X >= x)
# -------------------------------------------------------------------------

#' @rdname LUB-reliability
#' @export
mrlLUB <- function(
    x,
    size,
    beta) {

  .check_LUB(size, beta)

  x <- as.numeric(x)

  probs <- .LUB_probs(
    size = size,
    beta = beta
  )

  if (is.null(probs)) {
    return(rep(NA_real_, length(x)))
  }

  support <- 0:size
  surv <- rev(cumsum(rev(probs)))
  surv[1L] <- 1

  out <- vapply(
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

      index <- as.integer(xx) + 1L
      sx <- surv[index]

      if (!is.finite(sx) || sx <= 0) {
        return(NA_real_)
      }

      sum(
        (support[(index + 1L):(size + 1L)] - xx) *
          probs[(index + 1L):(size + 1L)]
      ) / sx
    },
    numeric(1)
  )

  out
}
