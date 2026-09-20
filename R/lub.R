#' Lambert-Uniform Binomial Distribution
#'
#' Probability mass function, distribution function, quantile function,
#' and random generation for the Lambert-uniform binomial (LUB)
#' distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Non-negative integer specifying the number of observations
#'   to generate.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param beta Positive Lambert-uniform parameter.
#' @param log Logical; if `TRUE`, probabilities are returned on the
#'   logarithmic scale.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are
#'   \eqn{P(X \leq x)}; otherwise, \eqn{P(X > x)}.
#' @param log.p Logical; if `TRUE`, probabilities are given or returned
#'   on the logarithmic scale.
#' @param tol Relative tolerance used by the internal evaluation of the
#'   confluent hypergeometric function.
#' @param maxit Maximum number of terms used in the hypergeometric series.
#'
#' @details
#' Let \eqn{P} have survival function
#' \deqn{
#' \bar F_P(p;\beta)
#' =
#' (1-p)\exp\{(1-\beta)p\},
#' \qquad 0<p<1,
#' }
#' with \eqn{\beta>0}. Its density is
#' \deqn{
#' f_P(p;\beta)
#' =
#' \exp\{(1-\beta)p\}
#' [p+\beta(1-p)].
#' }
#'
#' Conditionally on \eqn{P=p}, let
#' \deqn{
#' X \mid P=p \sim \mathrm{Binomial}(size,p).
#' }
#' Then \eqn{X} follows the Lambert-uniform binomial distribution.
#'
#' The probability mass function is
#' \deqn{
#' P(X=x)
#' =
#' \frac{1}{size+1}
#' \left[
#' \beta\,{}_1F_1(x+1;size+2;1-\beta)
#' +
#' (1-\beta)\frac{x+1}{size+2}
#' {}_1F_1(x+2;size+3;1-\beta)
#' \right].
#' }
#'
#' When \eqn{\beta=1}, the LUB distribution reduces to the discrete
#' uniform distribution on \eqn{\{0,\ldots,size\}}:
#' \deqn{
#' P(X=x)=\frac{1}{size+1}.
#' }
#'
#' @return
#' `dLUB()` returns the probability mass function,
#' `pLUB()` returns the distribution function,
#' `qLUB()` returns the quantile function, and
#' `rLUB()` returns a random sample.
#'
#' @examples
#' dLUB(0:4, size = 4, beta = 3)
#' pLUB(2, size = 4, beta = 3)
#' qLUB(0.5, size = 4, beta = 3)
#'
#' set.seed(123)
#' rLUB(10, size = 4, beta = 3)
#'
#' @name LUB
NULL


# -------------------------------------------------------------------------
# Numerical evaluation through the latent integral representation
# -------------------------------------------------------------------------

.dLUB_integral_scalar <- function(
    x,
    size,
    beta,
    rel.tol = 1e-10) {

  if (beta == 1) {
    return(1 / (size + 1))
  }

  integrand <- function(p) {

    log_kernel <-
      stats::dbinom(
        x,
        size = size,
        prob = p,
        log = TRUE
      ) +
      .logdLU(
        p,
        beta = beta
      )

    exp(log_kernel)
  }

  ans <- stats::integrate(
    integrand,
    lower = 0,
    upper = 1,
    rel.tol = rel.tol,
    subdivisions = 200L,
    stop.on.error = FALSE
  )

  if (!is.list(ans) ||
      !is.finite(ans$value) ||
      ans$value < 0) {
    return(NA_real_)
  }

  ans$value
}


.dLUB_integral <- function(
    x,
    size,
    beta,
    log = FALSE,
    rel.tol = 1e-10) {

  .check_LUB(size, beta)
  .check_logical_scalar(log, "log")
  .check_positive_scalar(rel.tol, "rel.tol")

  x <- as.numeric(x)

  out <- rep(
    if (log) -Inf else 0,
    length(x)
  )

  in_support <-
    is.finite(x) &
    x >= 0 &
    x <= size &
    x == floor(x)

  if (any(in_support)) {

    vals <- vapply(
      x[in_support],
      .dLUB_integral_scalar,
      numeric(1),
      size = size,
      beta = beta,
      rel.tol = rel.tol
    )

    bad <- !is.finite(vals) | vals < 0

    if (any(bad)) {
      vals[bad] <- NA_real_
    }

    out[in_support] <- if (log) {
      log(vals)
    } else {
      vals
    }
  }

  out[is.na(x)] <- NA_real_

  out
}


# -------------------------------------------------------------------------
# Numerical evaluation through the closed-form hypergeometric expression
# -------------------------------------------------------------------------

.dLUB_hyper_scalar <- function(
    x,
    size,
    beta,
    tol = 1e-14,
    maxit = 10000L) {

  if (beta == 1) {
    return(1 / (size + 1))
  }

  z <- 1 - beta

  M1 <- .hyp1f1_series(
    a = x + 1,
    b = size + 2,
    z = z,
    tol = tol,
    maxit = maxit
  )

  M2 <- .hyp1f1_series(
    a = x + 2,
    b = size + 3,
    z = z,
    tol = tol,
    maxit = maxit
  )

  value <- (
    beta * M1 +
      (1 - beta) *
        (x + 1) / (size + 2) *
        M2
  ) / (size + 1)

  if (!is.finite(value)) {
    return(NA_real_)
  }

  # Remove negligible negative values caused by floating-point error.
  if (value < 0 &&
      value > -100 * .Machine$double.eps) {
    value <- 0
  }

  value
}


.dLUB_hyper <- function(
    x,
    size,
    beta,
    log = FALSE,
    tol = 1e-14,
    maxit = 10000L) {

  .check_LUB(size, beta)
  .check_logical_scalar(log, "log")
  .check_LUB_control(tol, maxit)

  x <- as.numeric(x)

  out <- rep(
    if (log) -Inf else 0,
    length(x)
  )

  in_support <-
    is.finite(x) &
    x >= 0 &
    x <= size &
    x == floor(x)

  if (any(in_support)) {

    vals <- vapply(
      x[in_support],
      .dLUB_hyper_scalar,
      numeric(1),
      size = size,
      beta = beta,
      tol = tol,
      maxit = maxit
    )

    bad <- !is.finite(vals) | vals < 0

    if (any(bad)) {
      vals[bad] <- NA_real_
    }

    out[in_support] <- if (log) {
      log(vals)
    } else {
      vals
    }
  }

  out[is.na(x)] <- NA_real_

  out
}


# -------------------------------------------------------------------------
# Probability mass function
# -------------------------------------------------------------------------

#' @rdname LUB
#' @export
dLUB <- function(
    x,
    size,
    beta,
    log = FALSE,
    tol = 1e-14,
    maxit = 10000L) {

  .check_LUB(size, beta)
  .check_logical_scalar(log, "log")
  .check_LUB_control(tol, maxit)

  x <- as.numeric(x)

  # For strongly negative z = 1 - beta, evaluation of the
  # hypergeometric series becomes numerically delicate. In this region,
  # the latent integral representation is slower but more stable.
  use_integral <- beta > 5.5

  if (use_integral) {
    return(
      .dLUB_integral(
        x = x,
        size = size,
        beta = beta,
        log = log
      )
    )
  }

  out <- .dLUB_hyper(
    x = x,
    size = size,
    beta = beta,
    log = log,
    tol = tol,
    maxit = maxit
  )

  in_support <-
    is.finite(x) &
    x >= 0 &
    x <= size &
    x == floor(x)

  bad <- if (log) {

    in_support &
      (
        is.na(out) |
          !is.finite(out) |
          out > 100 * .Machine$double.eps
      )

  } else {

    in_support &
      (
        is.na(out) |
          !is.finite(out) |
          out < 0 |
          out > 1 + 100 * .Machine$double.eps
      )
  }

  if (any(bad)) {

    out[bad] <- .dLUB_integral(
      x = x[bad],
      size = size,
      beta = beta,
      log = log
    )
  }

  out
}


# -------------------------------------------------------------------------
# Internal normalized probability vector
# -------------------------------------------------------------------------

.LUB_probs <- function(
    size,
    beta,
    tol = 1e-14,
    maxit = 10000L) {

  .check_LUB(size, beta)
  .check_LUB_control(tol, maxit)

  probs <- dLUB(
    x = 0:size,
    size = size,
    beta = beta,
    tol = tol,
    maxit = maxit
  )

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
# Distribution function
# -------------------------------------------------------------------------

#' @rdname LUB
#' @export
pLUB <- function(
    q,
    size,
    beta,
    lower.tail = TRUE,
    log.p = FALSE,
    tol = 1e-14,
    maxit = 10000L) {

  .check_LUB(size, beta)
  .check_logical_scalar(lower.tail, "lower.tail")
  .check_logical_scalar(log.p, "log.p")
  .check_LUB_control(tol, maxit)

  q <- as.numeric(q)

  probs <- .LUB_probs(
    size = size,
    beta = beta,
    tol = tol,
    maxit = maxit
  )

  if (is.null(probs)) {
    return(rep(NA_real_, length(q)))
  }

  cdf <- cumsum(probs)
  cdf[length(cdf)] <- 1

  out <- rep(NA_real_, length(q))

  for (i in seq_along(q)) {

    qi <- q[i]

    if (is.na(qi)) {
      next
    }

    if (qi < 0) {

      Fq <- 0

    } else if (qi >= size) {

      Fq <- 1

    } else {

      Fq <- cdf[floor(qi) + 1L]
    }

    if (!lower.tail) {
      Fq <- 1 - Fq
    }

    out[i] <- if (log.p) {
      log(Fq)
    } else {
      Fq
    }
  }

  out
}


# -------------------------------------------------------------------------
# Quantile function
# -------------------------------------------------------------------------

#' @rdname LUB
#' @export
qLUB <- function(
    p,
    size,
    beta,
    lower.tail = TRUE,
    log.p = FALSE,
    tol = 1e-14,
    maxit = 10000L) {

  .check_LUB(size, beta)
  .check_logical_scalar(lower.tail, "lower.tail")
  .check_logical_scalar(log.p, "log.p")
  .check_LUB_control(tol, maxit)

  p <- as.numeric(p)

  if (log.p) {
    p <- exp(p)
  }

  if (!lower.tail) {
    p <- 1 - p
  }

  out <- rep(NA_real_, length(p))

  valid <- !is.na(p) & p >= 0 & p <= 1

  if (any(!valid & !is.na(p))) {
    warning(
      "NaNs produced",
      call. = FALSE
    )
  }

  probs <- .LUB_probs(
    size = size,
    beta = beta,
    tol = tol,
    maxit = maxit
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
        which(cdf >= p[i])[1L] - 1L
    }
  }

  out
}


# -------------------------------------------------------------------------
# Random generation
# -------------------------------------------------------------------------

#' @rdname LUB
#' @export
rLUB <- function(
    n,
    size,
    beta,
    tol = 1e-14,
    maxit = 10000L) {

  .check_LUB(size, beta)
  .check_LUB_control(tol, maxit)

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

  probs <- .LUB_probs(
    size = size,
    beta = beta,
    tol = tol,
    maxit = maxit
  )

  if (is.null(probs)) {
    stop(
      "Unable to evaluate the LUB probability mass function.",
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
