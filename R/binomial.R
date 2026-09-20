# -------------------------------------------------------------------------
# Binomial distribution interface
# -------------------------------------------------------------------------


.check_BIN <- function(size, prob) {

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

  if (length(prob) != 1L ||
      !is.numeric(prob) ||
      !is.finite(prob) ||
      prob < 0 ||
      prob > 1) {
    stop(
      "'prob' must lie in [0, 1].",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Binomial Distribution Interface
#'
#' Convenience wrappers around the standard binomial distribution functions
#' provided by base R. These functions are included to provide a common
#' interface across the count models used by the package.
#'
#' @param x Numeric vector of quantiles.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Non-negative integer specifying the number of observations
#'   to generate.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param prob Binomial success probability in \eqn{[0,1]}.
#' @param log Logical; if `TRUE`, probabilities are returned on the
#'   logarithmic scale.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are
#'   \eqn{P(X \leq x)}; otherwise, \eqn{P(X > x)}.
#' @param log.p Logical; if `TRUE`, probabilities are given or returned
#'   on the logarithmic scale.
#'
#' @return
#' `dBIN()` returns the probability mass function,
#' `pBIN()` returns the distribution function,
#' `qBIN()` returns the quantile function, and
#' `rBIN()` returns a random sample.
#'
#' @examples
#' dBIN(0:4, size = 4, prob = 0.25)
#' pBIN(2, size = 4, prob = 0.25)
#' qBIN(0.5, size = 4, prob = 0.25)
#'
#' set.seed(123)
#' rBIN(10, size = 4, prob = 0.25)
#'
#' @name BIN
NULL


#' @rdname BIN
#' @export
dBIN <- function(
    x,
    size,
    prob,
    log = FALSE) {

  .check_BIN(size, prob)
  .check_logical_scalar(log, "log")

  stats::dbinom(
    x = x,
    size = size,
    prob = prob,
    log = log
  )
}


#' @rdname BIN
#' @export
pBIN <- function(
    q,
    size,
    prob,
    lower.tail = TRUE,
    log.p = FALSE) {

  .check_BIN(size, prob)
  .check_logical_scalar(lower.tail, "lower.tail")
  .check_logical_scalar(log.p, "log.p")

  stats::pbinom(
    q = q,
    size = size,
    prob = prob,
    lower.tail = lower.tail,
    log.p = log.p
  )
}


#' @rdname BIN
#' @export
qBIN <- function(
    p,
    size,
    prob,
    lower.tail = TRUE,
    log.p = FALSE) {

  .check_BIN(size, prob)
  .check_logical_scalar(lower.tail, "lower.tail")
  .check_logical_scalar(log.p, "log.p")

  stats::qbinom(
    p = p,
    size = size,
    prob = prob,
    lower.tail = lower.tail,
    log.p = log.p
  )
}


#' @rdname BIN
#' @export
rBIN <- function(
    n,
    size,
    prob) {

  .check_BIN(size, prob)

  if (length(n) != 1L ||
      !is.numeric(n) ||
      !is.finite(n) ||
      n < 0 ||
      n != floor(n)) {
    stop(
      "'n' must be a non-negative integer.",
      call. = FALSE
    )
  }

  stats::rbinom(
    n = n,
    size = size,
    prob = prob
  )
}
