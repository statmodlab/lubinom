# -------------------------------------------------------------------------
# Conway-Maxwell-binomial distribution
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Parameter validation
# -------------------------------------------------------------------------

.check_CMB <- function(
    size,
    prob,
    nu) {

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
      prob <= 0 ||
      prob >= 1) {
    stop(
      "'prob' must satisfy 0 < prob < 1.",
      call. = FALSE
    )
  }

  .check_positive_scalar(
    nu,
    "nu"
  )

  invisible(TRUE)
}


# -------------------------------------------------------------------------
# Unnormalized log-probability
# -------------------------------------------------------------------------

.logCMB_unnormalized <- function(
    x,
    size,
    prob,
    nu) {

  nu *
    lchoose(
      size,
      x
    ) +
    x * log(prob) +
    (size - x) *
      log1p(-prob)
}


# -------------------------------------------------------------------------
# Log normalizing constant
# -------------------------------------------------------------------------

.CMB_log_normalizer <- function(
    size,
    prob,
    nu) {

  .check_CMB(
    size = size,
    prob = prob,
    nu = nu
  )

  support <- 0:size

  logw <- .logCMB_unnormalized(
    x = support,
    size = size,
    prob = prob,
    nu = nu
  )

  m <- max(logw)

  if (!is.finite(m)) {
    return(NA_real_)
  }

  m +
    log(
      sum(
        exp(logw - m)
      )
    )
}


# -------------------------------------------------------------------------
# Internal normalized probability vector
# -------------------------------------------------------------------------

.CMB_probs <- function(
    size,
    prob,
    nu) {

  .check_CMB(
    size = size,
    prob = prob,
    nu = nu
  )

  support <- 0:size

  logw <- .logCMB_unnormalized(
    x = support,
    size = size,
    prob = prob,
    nu = nu
  )

  m <- max(logw)

  if (!is.finite(m)) {
    return(NULL)
  }

  weights <- exp(
    logw - m
  )

  total <- sum(weights)

  if (any(!is.finite(weights)) ||
      any(weights < 0) ||
      !is.finite(total) ||
      total <= 0) {
    return(NULL)
  }

  weights / total
}


# -------------------------------------------------------------------------
# Documentation
# -------------------------------------------------------------------------

#' Conway-Maxwell-Binomial Distribution
#'
#' Probability mass function, distribution function, quantile function,
#' and random generation for the Conway-Maxwell-binomial distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Non-negative integer specifying the number of observations
#'   to generate.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param prob Probability parameter satisfying \eqn{0<prob<1}.
#' @param nu Positive dispersion parameter.
#' @param log Logical; if `TRUE`, probabilities are returned on the
#'   logarithmic scale.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are
#'   \eqn{P(X \leq x)}; otherwise, \eqn{P(X > x)}.
#' @param log.p Logical; if `TRUE`, probabilities are given or returned
#'   on the logarithmic scale.
#'
#' @details
#' The Conway-Maxwell-binomial distribution has probability mass function
#' \deqn{
#' P(X=x)
#' =
#' \frac{
#' {size \choose x}^{\nu}
#' prob^x(1-prob)^{size-x}
#' }{
#' Z(size,prob,\nu)
#' },
#' \qquad x=0,\ldots,size,
#' }
#' where
#' \deqn{
#' Z(size,prob,\nu)
#' =
#' \sum_{j=0}^{size}
#' {size \choose j}^{\nu}
#' prob^j(1-prob)^{size-j}.
#' }
#'
#' The normalizing constant is evaluated using a log-sum-exp
#' representation for numerical stability.
#'
#' When \eqn{\nu=1}, the distribution reduces to the ordinary binomial
#' distribution with parameters `size` and `prob`.
#'
#' @return
#' `dCMB()` returns the probability mass function,
#' `pCMB()` returns the distribution function,
#' `qCMB()` returns the quantile function, and
#' `rCMB()` returns a random sample.
#'
#' @examples
#' dCMB(0:4, size = 4, prob = 0.25, nu = 2)
#' pCMB(2, size = 4, prob = 0.25, nu = 2)
#' qCMB(0.5, size = 4, prob = 0.25, nu = 2)
#'
#' set.seed(123)
#' rCMB(10, size = 4, prob = 0.25, nu = 2)
#'
#' @name CMB
NULL


# -------------------------------------------------------------------------
# Probability mass function
# -------------------------------------------------------------------------

#' @rdname CMB
#' @export
dCMB <- function(
    x,
    size,
    prob,
    nu,
    log = FALSE) {

  .check_CMB(
    size = size,
    prob = prob,
    nu = nu
  )

  .check_logical_scalar(
    log,
    "log"
  )

  x <- as.numeric(x)

  out <- rep(
    if (log) -Inf else 0,
    length(x)
  )

  ok <-
    is.finite(x) &
    x >= 0 &
    x <= size &
    x == floor(x)

  if (any(ok)) {

    logZ <- .CMB_log_normalizer(
      size = size,
      prob = prob,
      nu = nu
    )

    if (!is.finite(logZ)) {
      out[ok] <- NA_real_
    } else {

      lp <-
        .logCMB_unnormalized(
          x = x[ok],
          size = size,
          prob = prob,
          nu = nu
        ) -
        logZ

      out[ok] <- if (log) {
        lp
      } else {
        exp(lp)
      }
    }
  }

  out[is.na(x)] <- NA_real_

  out
}


# -------------------------------------------------------------------------
# Distribution function
# -------------------------------------------------------------------------

#' @rdname CMB
#' @export
pCMB <- function(
    q,
    size,
    prob,
    nu,
    lower.tail = TRUE,
    log.p = FALSE) {

  .check_CMB(
    size = size,
    prob = prob,
    nu = nu
  )

  .check_logical_scalar(
    lower.tail,
    "lower.tail"
  )

  .check_logical_scalar(
    log.p,
    "log.p"
  )

  q <- as.numeric(q)

  probs <- .CMB_probs(
    size = size,
    prob = prob,
    nu = nu
  )

  if (is.null(probs)) {
    return(
      rep(
        NA_real_,
        length(q)
      )
    )
  }

  cdf <- cumsum(probs)
  cdf[length(cdf)] <- 1

  out <- rep(
    NA_real_,
    length(q)
  )

  for (i in seq_along(q)) {

    qi <- q[i]

    if (is.na(qi)) {
      next
    }

    if (qi < 0) {

      value <- 0

    } else if (qi >= size) {

      value <- 1

    } else {

      value <-
        cdf[
          floor(qi) + 1L
        ]
    }

    if (!lower.tail) {
      value <- 1 - value
    }

    out[i] <- if (log.p) {
      log(value)
    } else {
      value
    }
  }

  out
}


# -------------------------------------------------------------------------
# Quantile function
# -------------------------------------------------------------------------

#' @rdname CMB
#' @export
qCMB <- function(
    p,
    size,
    prob,
    nu,
    lower.tail = TRUE,
    log.p = FALSE) {

  .check_CMB(
    size = size,
    prob = prob,
    nu = nu
  )

  .check_logical_scalar(
    lower.tail,
    "lower.tail"
  )

  .check_logical_scalar(
    log.p,
    "log.p"
  )

  p <- as.numeric(p)

  if (log.p) {
    p <- exp(p)
  }

  if (!lower.tail) {
    p <- 1 - p
  }

  out <- rep(
    NA_real_,
    length(p)
  )

  valid <-
    !is.na(p) &
    p >= 0 &
    p <= 1

  if (any(
    !valid &
      !is.na(p)
  )) {
    warning(
      "NaNs produced",
      call. = FALSE
    )
  }

  probs <- .CMB_probs(
    size = size,
    prob = prob,
    nu = nu
  )

  if (is.null(probs)) {
    return(out)
  }

  cdf <- cumsum(probs)
  cdf[length(cdf)] <- 1

  for (i in which(valid)) {

    if (p[i] == 0) {

      out[i] <- 0

    } else {

      out[i] <-
        which(
          cdf >= p[i]
        )[1L] - 1L
    }
  }

  out
}


# -------------------------------------------------------------------------
# Random generation
# -------------------------------------------------------------------------

#' @rdname CMB
#' @export
rCMB <- function(
    n,
    size,
    prob,
    nu) {

  .check_CMB(
    size = size,
    prob = prob,
    nu = nu
  )

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

  if (n == 0L) {
    return(integer(0))
  }

  probs <- .CMB_probs(
    size = size,
    prob = prob,
    nu = nu
  )

  if (is.null(probs)) {
    stop(
      "Unable to evaluate the CMB probability mass function.",
      call. = FALSE
    )
  }

  sample.int(
    size + 1L,
    size = n,
    replace = TRUE,
    prob = probs
  ) - 1L
}
