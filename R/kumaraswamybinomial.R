# -------------------------------------------------------------------------
# Kumaraswamy-binomial distribution
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Parameter validation
# -------------------------------------------------------------------------

.check_KB <- function(
    size,
    a,
    b) {

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
    a,
    "a"
  )

  .check_positive_scalar(
    b,
    "b"
  )

  invisible(TRUE)
}


# -------------------------------------------------------------------------
# Stable signed log-sum
# -------------------------------------------------------------------------

.KB_signed_sum <- function(
    logabs,
    sign) {

  pos <- sign > 0
  neg <- sign < 0

  logsumexp <- function(z) {

    if (length(z) == 0L) {
      return(-Inf)
    }

    m <- max(z)

    if (!is.finite(m)) {
      return(m)
    }

    m + log(
      sum(
        exp(z - m)
      )
    )
  }

  lp <- logsumexp(
    logabs[pos]
  )

  ln <- logsumexp(
    logabs[neg]
  )

  if (!is.finite(lp) &&
      !is.finite(ln)) {
    return(0)
  }

  if (!is.finite(ln)) {
    return(exp(lp))
  }

  if (!is.finite(lp)) {
    return(-exp(ln))
  }

  if (lp < ln) {
    return(NA_real_)
  }

  if (lp == ln) {
    return(0)
  }

  exp(lp) *
    (
      -expm1(
        ln - lp
      )
    )
}


# -------------------------------------------------------------------------
# Finite-sum PMF evaluation
# -------------------------------------------------------------------------

.dKB_one_sum <- function(
    x,
    size,
    a,
    b) {

  j <- 0:(size - x)

  logabs <-
    lchoose(
      size - x,
      j
    ) +
    lbeta(
      1 + (x + j) / a,
      b
    )

  signs <- ifelse(
    j %% 2L == 0L,
    1,
    -1
  )

  s <- .KB_signed_sum(
    logabs = logabs,
    sign = signs
  )

  if (!is.finite(s)) {
    return(NA_real_)
  }

  value <-
    exp(
      lchoose(
        size,
        x
      ) +
        log(b)
    ) *
    s

  value
}


# -------------------------------------------------------------------------
# Integral PMF evaluation
# -------------------------------------------------------------------------

.dKB_one_integral <- function(
    x,
    size,
    a,
    b,
    rel.tol = 1e-10) {

  integrand <- function(p) {

    log_kernel <-
      stats::dbinom(
        x = x,
        size = size,
        prob = p,
        log = TRUE
      ) +
      log(a) +
      log(b) +
      (a - 1) * log(p) +
      (b - 1) *
        log1p(
          -p^a
        )

    exp(log_kernel)
  }

  z <- stats::integrate(
    f = integrand,
    lower = 0,
    upper = 1,
    subdivisions = 500L,
    rel.tol = rel.tol,
    stop.on.error = FALSE
  )

  if (!is.list(z) ||
      !is.finite(z$value) ||
      z$value < 0) {
    return(NA_real_)
  }

  z$value
}


# -------------------------------------------------------------------------
# Internal normalized probability vector
# -------------------------------------------------------------------------

.KB_probs <- function(
    size,
    a,
    b) {

  .check_KB(
    size = size,
    a = a,
    b = b
  )

  support <- 0:size

  probs <- vapply(
    support,
    .dKB_one_sum,
    numeric(1),
    size = size,
    a = a,
    b = b
  )

  tiny.negative <-
    is.finite(probs) &
    probs < 0 &
    probs > -1e-12

  probs[tiny.negative] <- 0

  bad <-
    !is.finite(probs) |
    probs < 0

  if (any(bad)) {

    probs[bad] <- vapply(
      support[bad],
      .dKB_one_integral,
      numeric(1),
      size = size,
      a = a,
      b = b
    )
  }

  total <- sum(probs)

  # A substantial departure from one indicates numerical cancellation
  # in the finite expansion. In that case, evaluate the complete support
  # by adaptive quadrature.
  if (any(!is.finite(probs)) ||
      any(probs < 0) ||
      !is.finite(total) ||
      total <= 0 ||
      abs(total - 1) > 1e-7) {

    probs <- vapply(
      support,
      .dKB_one_integral,
      numeric(1),
      size = size,
      a = a,
      b = b
    )

    total <- sum(probs)
  }

  if (any(!is.finite(probs)) ||
      any(probs < 0) ||
      !is.finite(total) ||
      total <= 0) {
    stop(
      "Failed to evaluate the Kumaraswamy-binomial pmf.",
      call. = FALSE
    )
  }

  probs / total
}


# -------------------------------------------------------------------------
# Documentation
# -------------------------------------------------------------------------

#' Kumaraswamy-Binomial Distribution
#'
#' Probability mass function, distribution function, quantile function,
#' and random generation for the Kumaraswamy-binomial distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Non-negative integer specifying the number of observations
#'   to generate.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param a Positive first Kumaraswamy shape parameter.
#' @param b Positive second Kumaraswamy shape parameter.
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
#' P \sim Kumaraswamy(a,b),
#' }
#' with density
#' \deqn{
#' f_P(p)
#' =
#' abp^{a-1}(1-p^a)^{b-1},
#' \qquad 0<p<1.
#' }
#'
#' Conditionally on \eqn{P=p},
#' \deqn{
#' X \mid P=p
#' \sim Binomial(size,p).
#' }
#'
#' The marginal probability mass function admits the finite expansion
#' \deqn{
#' P(X=x)
#' =
#' {size \choose x}
#' b
#' \sum_{j=0}^{size-x}
#' (-1)^j
#' {size-x \choose j}
#' B\left(
#' 1+\frac{x+j}{a},
#' b
#' \right).
#' }
#'
#' The implementation evaluates this expression first and uses adaptive
#' numerical integration as a fallback when numerical cancellation is
#' detected.
#'
#' When \eqn{a=b=1}, the mixing distribution is uniform on \eqn{(0,1)}
#' and the marginal distribution of \eqn{X} is discrete uniform on
#' \eqn{\{0,\ldots,size\}}.
#'
#' @return
#' `dKB()` returns the probability mass function,
#' `pKB()` returns the distribution function,
#' `qKB()` returns the quantile function, and
#' `rKB()` returns a random sample.
#'
#' @examples
#' dKB(0:4, size = 4, a = 1, b = 3)
#' pKB(2, size = 4, a = 1, b = 3)
#' qKB(0.5, size = 4, a = 1, b = 3)
#'
#' set.seed(123)
#' rKB(10, size = 4, a = 1, b = 3)
#'
#' @name KB
NULL


# -------------------------------------------------------------------------
# Probability mass function
# -------------------------------------------------------------------------

#' @rdname KB
#' @export
dKB <- function(
    x,
    size,
    a,
    b,
    log = FALSE) {

  .check_KB(
    size = size,
    a = a,
    b = b
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

    probs <- .KB_probs(
      size = size,
      a = a,
      b = b
    )

    vals <-
      probs[
        as.integer(
          x[ok]
        ) + 1L
      ]

    out[ok] <- if (log) {
      log(vals)
    } else {
      vals
    }
  }

  out[is.na(x)] <- NA_real_

  out
}


# -------------------------------------------------------------------------
# Distribution function
# -------------------------------------------------------------------------

#' @rdname KB
#' @export
pKB <- function(
    q,
    size,
    a,
    b,
    lower.tail = TRUE,
    log.p = FALSE) {

  .check_KB(
    size = size,
    a = a,
    b = b
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

  probs <- .KB_probs(
    size = size,
    a = a,
    b = b
  )

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

#' @rdname KB
#' @export
qKB <- function(
    p,
    size,
    a,
    b,
    lower.tail = TRUE,
    log.p = FALSE) {

  .check_KB(
    size = size,
    a = a,
    b = b
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

  probs <- .KB_probs(
    size = size,
    a = a,
    b = b
  )

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

#' @rdname KB
#' @export
rKB <- function(
    n,
    size,
    a,
    b) {

  .check_KB(
    size = size,
    a = a,
    b = b
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

  u <- stats::runif(n)

  prob <- (
    1 -
      (1 - u)^(1 / b)
  )^(1 / a)

  stats::rbinom(
    n = n,
    size = size,
    prob = prob
  )
}
