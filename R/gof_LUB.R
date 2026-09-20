# -------------------------------------------------------------------------
# Parametric-bootstrap goodness-of-fit for the LUB distribution
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Pearson statistic
# -------------------------------------------------------------------------

.LUB_pearson_stat <- function(obs, exp) {

  if (length(obs) != length(exp) ||
      any(!is.finite(obs)) ||
      any(!is.finite(exp)) ||
      any(obs < 0) ||
      any(exp <= 0)) {
    return(NA_real_)
  }

  sum(
    (obs - exp)^2 / exp
  )
}


# -------------------------------------------------------------------------
# Likelihood-ratio deviance
# -------------------------------------------------------------------------

.LUB_deviance_stat <- function(obs, exp) {

  if (length(obs) != length(exp) ||
      any(!is.finite(obs)) ||
      any(!is.finite(exp)) ||
      any(obs < 0) ||
      any(exp <= 0)) {
    return(NA_real_)
  }

  keep <- obs > 0

  if (!any(keep)) {
    return(0)
  }

  2 * sum(
    obs[keep] *
      log(
        obs[keep] /
          exp[keep]
      )
  )
}


# -------------------------------------------------------------------------
# Parametric-bootstrap goodness-of-fit
# -------------------------------------------------------------------------

#' Parametric-Bootstrap Goodness-of-Fit for the LUB Distribution
#'
#' Performs a parametric-bootstrap goodness-of-fit assessment for the
#' Lambert-uniform binomial distribution using Pearson and likelihood-ratio
#' deviance statistics.
#'
#' @param x Integer vector of observations from
#'   \eqn{\{0,1,\ldots,\texttt{size}\}}.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param B Positive integer giving the number of bootstrap replications.
#' @param seed Optional integer seed for reproducibility.
#' @param progress Logical; if `TRUE`, progress is printed every 100
#'   bootstrap replications.
#'
#' @details
#' The fitted null model is
#' \eqn{LUB(size,\widehat{\beta})}.
#'
#' For each bootstrap replication, a sample with the same number of
#' observations is generated from the fitted LUB model. The parameter
#' \eqn{\beta} is then re-estimated and the discrepancy statistics are
#' recomputed using the refitted probabilities.
#'
#' The Pearson statistic is
#' \deqn{
#' X^2 =
#' \sum_x
#' \frac{(O_x-E_x)^2}{E_x},
#' }
#' and the likelihood-ratio deviance is
#' \deqn{
#' G^2 =
#' 2
#' \sum_{x:O_x>0}
#' O_x
#' \log\left(\frac{O_x}{E_x}\right).
#' }
#'
#' For either discrepancy statistic \eqn{T}, the bootstrap p-value is
#' computed as
#' \deqn{
#' \widehat p =
#' \frac{
#' 1+\sum_{b=1}^{B_{\mathrm{valid}}}
#' I(T_b^* \geq T_{\mathrm{obs}})
#' }{
#' B_{\mathrm{valid}}+1
#' }.
#' }
#'
#' Bootstrap samples for which the LUB model cannot be refitted are
#' excluded from the valid bootstrap replications and are reported through
#' `B.valid` and `valid.rate`.
#'
#' @return
#' An object of class `"bootstrapGOF_LUB"` containing the fitted LUB model,
#' observed and expected frequencies, Pearson and deviance statistics,
#' bootstrap p-values, Monte Carlo standard errors, and bootstrap
#' diagnostics.
#'
#' @examples
#' \donttest{
#' x <- rep(
#'   0:4,
#'   times = c(69, 30, 9, 2, 0)
#' )
#'
#' gof <- bootstrapGOF_LUB(
#'   x = x,
#'   size = 4,
#'   B = 499,
#'   seed = 123,
#'   progress = FALSE
#' )
#'
#' print(gof)
#' }
#'
#' @export
bootstrapGOF_LUB <- function(
    x,
    size,
    B = 1000L,
    seed = NULL,
    progress = TRUE) {

  x <- as.numeric(x)

  .check_LUB_sample(
    x = x,
    size = size
  )

  if (length(B) != 1L ||
      !is.numeric(B) ||
      !is.finite(B) ||
      B < 1 ||
      B != floor(B)) {
    stop(
      "'B' must be a positive integer.",
      call. = FALSE
    )
  }

  .check_logical_scalar(
    progress,
    "progress"
  )

  if (!is.null(seed)) {

    if (length(seed) != 1L ||
        !is.numeric(seed) ||
        !is.finite(seed) ||
        seed != floor(seed)) {
      stop(
        "'seed' must be NULL or a finite integer.",
        call. = FALSE
      )
    }

    set.seed(seed)
  }

  x <- as.integer(x)

  N <- length(x)
  support <- 0:size


  # -----------------------------------------------------------------------
  # Fit observed sample
  # -----------------------------------------------------------------------

  fit.obs <- fitLUB_mle(
    x = x,
    size = size,
    diagnostic_bfgs = FALSE,
    profile.ci = FALSE
  )

  beta.hat <- unname(
    fit.obs$par["beta"]
  )

  probs.obs <- .LUB_probs(
    size = size,
    beta = beta.hat
  )

  if (is.null(probs.obs)) {
    stop(
      "Unable to evaluate fitted LUB probabilities.",
      call. = FALSE
    )
  }

  freq.obs <- tabulate(
    x + 1L,
    nbins = size + 1L
  )

  exp.obs <- N * probs.obs


  # -----------------------------------------------------------------------
  # Observed discrepancy statistics
  # -----------------------------------------------------------------------

  pearson.obs <- .LUB_pearson_stat(
    obs = freq.obs,
    exp = exp.obs
  )

  deviance.obs <- .LUB_deviance_stat(
    obs = freq.obs,
    exp = exp.obs
  )

  if (!is.finite(pearson.obs) ||
      !is.finite(deviance.obs)) {
    stop(
      "Unable to compute the observed goodness-of-fit statistics.",
      call. = FALSE
    )
  }


  # -----------------------------------------------------------------------
  # Bootstrap storage
  # -----------------------------------------------------------------------

  pearson.boot <- rep(
    NA_real_,
    B
  )

  deviance.boot <- rep(
    NA_real_,
    B
  )

  beta.boot <- rep(
    NA_real_,
    B
  )

  valid <- rep(
    FALSE,
    B
  )


  # -----------------------------------------------------------------------
  # Parametric bootstrap
  # -----------------------------------------------------------------------

  for (b in seq_len(B)) {

    xb <- rLUB(
      n = N,
      size = size,
      beta = beta.hat
    )

    fit.b <- try(
      fitLUB_mle(
        x = xb,
        size = size,
        diagnostic_bfgs = FALSE,
        profile.ci = FALSE
      ),
      silent = TRUE
    )

    if (!inherits(fit.b, "try-error")) {

      beta.b <- unname(
        fit.b$par["beta"]
      )

      if (is.finite(beta.b) &&
          beta.b > 0) {

        probs.b <- .LUB_probs(
          size = size,
          beta = beta.b
        )

        if (!is.null(probs.b)) {

          freq.b <- tabulate(
            xb + 1L,
            nbins = size + 1L
          )

          exp.b <- N * probs.b

          pearson.b <- .LUB_pearson_stat(
            obs = freq.b,
            exp = exp.b
          )

          deviance.b <- .LUB_deviance_stat(
            obs = freq.b,
            exp = exp.b
          )

          if (is.finite(pearson.b) &&
              is.finite(deviance.b)) {

            beta.boot[b] <- beta.b
            pearson.boot[b] <- pearson.b
            deviance.boot[b] <- deviance.b

            valid[b] <- TRUE
          }
        }
      }
    }

    if (progress &&
        (b %% 100L == 0L || b == B)) {

      cat(
        "LUB:",
        b,
        "of",
        B,
        "bootstrap samples\n"
      )
    }
  }


  # -----------------------------------------------------------------------
  # Valid bootstrap replications
  # -----------------------------------------------------------------------

  id <- which(valid)

  B.valid <- length(id)

  if (B.valid == 0L) {
    stop(
      "No valid bootstrap samples were obtained.",
      call. = FALSE
    )
  }


  # -----------------------------------------------------------------------
  # Bootstrap p-values
  # -----------------------------------------------------------------------

  p.pearson <- (
    1 +
      sum(
        pearson.boot[id] >=
          pearson.obs
      )
  ) / (
    B.valid + 1
  )

  p.deviance <- (
    1 +
      sum(
        deviance.boot[id] >=
          deviance.obs
      )
  ) / (
    B.valid + 1
  )


  # -----------------------------------------------------------------------
  # Monte Carlo standard errors
  # -----------------------------------------------------------------------

  mcse.pearson <- sqrt(
    p.pearson *
      (1 - p.pearson) /
      (B.valid + 1)
  )

  mcse.deviance <- sqrt(
    p.deviance *
      (1 - p.deviance) /
      (B.valid + 1)
  )


  # -----------------------------------------------------------------------
  # Output
  # -----------------------------------------------------------------------

  out <- list(
    call = match.call(),

    fit = fit.obs,

    parameters = c(
      beta = beta.hat
    ),

    observed = freq.obs,
    expected = exp.obs,

    pearson = c(
      statistic = pearson.obs,
      p.value = p.pearson,
      mcse = mcse.pearson
    ),

    deviance = c(
      statistic = deviance.obs,
      p.value = p.deviance,
      mcse = mcse.deviance
    ),

    B = B,
    B.valid = B.valid,
    B.invalid = B - B.valid,
    valid.rate = B.valid / B,

    beta.boot = beta.boot,
    pearson.boot = pearson.boot,
    deviance.boot = deviance.boot,
    valid = valid,

    size = size,
    n = N,
    nobs = N,
    support = support,

    seed = seed
  )

  class(out) <- "bootstrapGOF_LUB"

  out
}


# -------------------------------------------------------------------------
# Print method
# -------------------------------------------------------------------------

#' @export
print.bootstrapGOF_LUB <- function(
    x,
    digits = 6,
    ...) {

  cat(
    "Parametric-bootstrap goodness-of-fit for the LUB model\n",
    "------------------------------------------------------\n",
    sep = ""
  )

  cat(
    "Observations :",
    x$nobs,
    "\n"
  )

  cat(
    "Size         :",
    x$size,
    "\n"
  )

  cat(
    "Beta hat     :",
    format(
      x$parameters["beta"],
      digits = digits
    ),
    "\n"
  )

  cat(
    "Replications :",
    x$B,
    "\n"
  )

  cat(
    "Valid        :",
    x$B.valid,
    "\n"
  )

  cat(
    "Invalid      :",
    x$B.invalid,
    "\n"
  )

  cat(
    "Valid rate   :",
    format(
      x$valid.rate,
      digits = digits
    ),
    "\n\n"
  )

  cat(
    "Pearson X2   :",
    format(
      x$pearson["statistic"],
      digits = digits
    ),
    "\n"
  )

  cat(
    "Bootstrap p  :",
    format(
      x$pearson["p.value"],
      digits = digits
    ),
    "\n"
  )

  cat(
    "Monte Carlo SE:",
    format(
      x$pearson["mcse"],
      digits = digits
    ),
    "\n\n"
  )

  cat(
    "Deviance G2  :",
    format(
      x$deviance["statistic"],
      digits = digits
    ),
    "\n"
  )

  cat(
    "Bootstrap p  :",
    format(
      x$deviance["p.value"],
      digits = digits
    ),
    "\n"
  )

  cat(
    "Monte Carlo SE:",
    format(
      x$deviance["mcse"],
      digits = digits
    ),
    "\n"
  )

  invisible(x)
}
