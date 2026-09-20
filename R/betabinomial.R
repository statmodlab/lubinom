# -------------------------------------------------------------------------
# Beta-binomial distribution
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Parameter validation
# -------------------------------------------------------------------------

.check_BB_parameters <- function(
    size,
    alpha,
    beta) {

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

  .check_positive_scalar(
    alpha,
    "alpha"
  )

  .check_positive_scalar(
    beta,
    "beta"
  )

  invisible(TRUE)
}


# -------------------------------------------------------------------------
# Internal normalized probability vector
# -------------------------------------------------------------------------

.BB_probs <- function(
    size,
    alpha,
    beta) {

  .check_BB_parameters(
    size = size,
    alpha = alpha,
    beta = beta
  )

  x <- 0:size

  logpmf <-
    lchoose(size, x) +
    lbeta(
      x + alpha,
      size - x + beta
    ) -
    lbeta(alpha, beta)

  probs <- exp(logpmf)

  if (any(!is.finite(probs)) ||
      any(probs < 0)) {
    return(NULL)
  }

  total <- sum(probs)

  if (!is.finite(total) ||
      total <= 0) {
    return(NULL)
  }

  probs / total
}


# -------------------------------------------------------------------------
# Documentation
# -------------------------------------------------------------------------

#' Beta-Binomial Distribution
#'
#' Probability mass function, distribution function, quantile function,
#' and random generation for the beta-binomial distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Non-negative integer specifying the number of observations
#'   to generate.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param alpha Positive first shape parameter of the beta mixing
#'   distribution.
#' @param beta Positive second shape parameter of the beta mixing
#'   distribution.
#' @param log Logical; if `TRUE`, probabilities are returned on the
#'   logarithmic scale.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are
#'   \eqn{P(X \leq x)}; otherwise, \eqn{P(X > x)}.
#' @param log.p Logical; if `TRUE`, probabilities are given or returned
#'   on the logarithmic scale.
#'
#' @details
#' Let
#' \deqn{
#' P \sim Beta(\alpha,\beta),
#' }
#' and, conditionally on \eqn{P=p},
#' \deqn{
#' X \mid P=p \sim Binomial(size,p).
#' }
#' Then \eqn{X} follows a beta-binomial distribution with probability
#' mass function
#' \deqn{
#' P(X=x)
#' =
#' {size \choose x}
#' \frac{
#' B(x+\alpha,size-x+\beta)
#' }{
#' B(\alpha,\beta)
#' },
#' \qquad x=0,\ldots,size.
#' }
#'
#' @return
#' `dBB()` returns the probability mass function,
#' `pBB()` returns the distribution function,
#' `qBB()` returns the quantile function, and
#' `rBB()` returns a random sample.
#'
#' @examples
#' dBB(0:4, size = 4, alpha = 1, beta = 3)
#' pBB(2, size = 4, alpha = 1, beta = 3)
#' qBB(0.5, size = 4, alpha = 1, beta = 3)
#'
#' set.seed(123)
#' rBB(10, size = 4, alpha = 1, beta = 3)
#'
#' @name BB
NULL


# -------------------------------------------------------------------------
# Probability mass function
# -------------------------------------------------------------------------

#' @rdname BB
#' @export
dBB <- function(
    x,
    size,
    alpha,
    beta,
    log = FALSE) {

  .check_BB_parameters(
    size = size,
    alpha = alpha,
    beta = beta
  )

  .check_logical_scalar(
    log,
    "log"
  )

  x <- as.numeric(x)

  ans <- rep(
    if (log) -Inf else 0,
    length(x)
  )

  ok <-
    is.finite(x) &
    x >= 0 &
    x <= size &
    x == floor(x)

  if (any(ok)) {

    logpmf <-
      lchoose(
        size,
        x[ok]
      ) +
      lbeta(
        x[ok] + alpha,
        size - x[ok] + beta
      ) -
      lbeta(
        alpha,
        beta
      )

    ans[ok] <- if (log) {
      logpmf
    } else {
      exp(logpmf)
    }
  }

  ans[is.na(x)] <- NA_real_

  ans
}


# -------------------------------------------------------------------------
# Distribution function
# -------------------------------------------------------------------------

#' @rdname BB
#' @export
pBB <- function(
    q,
    size,
    alpha,
    beta,
    lower.tail = TRUE,
    log.p = FALSE) {

  .check_BB_parameters(
    size = size,
    alpha = alpha,
    beta = beta
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

  probs <- .BB_probs(
    size = size,
    alpha = alpha,
    beta = beta
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

  ans <- rep(
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

    ans[i] <- if (log.p) {
      log(value)
    } else {
      value
    }
  }

  ans
}


# -------------------------------------------------------------------------
# Quantile function
# -------------------------------------------------------------------------

#' @rdname BB
#' @export
qBB <- function(
    p,
    size,
    alpha,
    beta,
    lower.tail = TRUE,
    log.p = FALSE) {

  .check_BB_parameters(
    size = size,
    alpha = alpha,
    beta = beta
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

  ans <- rep(
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

  probs <- .BB_probs(
    size = size,
    alpha = alpha,
    beta = beta
  )

  if (is.null(probs)) {
    return(ans)
  }

  cdf <- cumsum(probs)
  cdf[length(cdf)] <- 1

  for (i in which(valid)) {

    if (p[i] == 0) {

      ans[i] <- 0

    } else {

      ans[i] <-
        which(
          cdf >= p[i]
        )[1L] - 1L
    }
  }

  ans
}


# -------------------------------------------------------------------------
# Random generation
# -------------------------------------------------------------------------

#' @rdname BB
#' @export
rBB <- function(
    n,
    size,
    alpha,
    beta) {

  .check_BB_parameters(
    size = size,
    alpha = alpha,
    beta = beta
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

  prob <- stats::rbeta(
    n = n,
    shape1 = alpha,
    shape2 = beta
  )

  stats::rbinom(
    n = n,
    size = size,
    prob = prob
  )
}
