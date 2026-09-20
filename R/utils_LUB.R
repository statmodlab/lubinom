# -------------------------------------------------------------------------
# Utility functions for the Lambert-uniform binomial (LUB) distribution
# -------------------------------------------------------------------------

.check_positive_scalar <- function(x, name) {

  if (length(x) != 1L ||
      !is.numeric(x) ||
      !is.finite(x) ||
      x <= 0) {
    stop(
      paste0("'", name, "' must be a strictly positive finite scalar."),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.check_logical_scalar <- function(x, name) {

  if (length(x) != 1L ||
      !is.logical(x) ||
      is.na(x)) {
    stop(
      paste0("'", name, "' must be TRUE or FALSE."),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.check_LUB <- function(size, beta) {

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

  .check_positive_scalar(beta, "beta")

  invisible(TRUE)
}


.check_LUB_control <- function(tol, maxit) {

  .check_positive_scalar(tol, "tol")

  if (length(maxit) != 1L ||
      !is.numeric(maxit) ||
      !is.finite(maxit) ||
      maxit < 1 ||
      maxit != floor(maxit)) {
    stop(
      "'maxit' must be a strictly positive integer.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# Latent Lambert-uniform log-density
.logdLU <- function(p, beta) {
  (1 - beta) * p + log(p + beta * (1 - p))
}


# Unconstrained parameterization: eta = log(beta)
.beta_from_eta <- function(eta) {

  if (length(eta) != 1L ||
      !is.numeric(eta) ||
      !is.finite(eta)) {
    return(NA_real_)
  }

  exp(eta)
}


# Confluent hypergeometric function 1F1(a; b; z)
.hyp1f1_series <- function(
    a,
    b,
    z,
    tol = 1e-14,
    maxit = 10000L) {

  if (length(a) != 1L ||
      !is.numeric(a) ||
      !is.finite(a) ||
      length(b) != 1L ||
      !is.numeric(b) ||
      !is.finite(b) ||
      b == 0 ||
      length(z) != 1L ||
      !is.numeric(z) ||
      !is.finite(z)) {
    stop(
      "'a', 'b', and 'z' must be finite numeric scalars, with b != 0.",
      call. = FALSE
    )
  }

  .check_LUB_control(tol, maxit)

  if (z == 0) {
    return(1)
  }

  term <- 1
  value <- 1

  for (k in seq_len(maxit)) {

    denom <- b + k - 1

    if (denom == 0) {
      return(NaN)
    }

    term <- term *
      (a + k - 1) / denom *
      z / k

    new_value <- value + term

    if (!is.finite(term) ||
        !is.finite(new_value)) {
      return(NaN)
    }

    if (abs(term) <= tol * max(1, abs(new_value))) {
      return(new_value)
    }

    value <- new_value
  }

  warning(
    "1F1 series did not converge within 'maxit' iterations.",
    call. = FALSE
  )

  value
}


# Finite penalty used by numerical optimizers
.LUB_penalty <- .Machine$double.xmax^0.25
