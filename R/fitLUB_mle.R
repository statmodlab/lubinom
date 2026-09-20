# -------------------------------------------------------------------------
# Maximum-likelihood estimation for the Lambert-uniform binomial distribution
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Sample validation
# -------------------------------------------------------------------------

.check_LUB_sample <- function(x, size) {

  .check_LUB(size = size, beta = 1)

  if (!is.numeric(x) || length(x) < 1L) {
    stop(
      "'x' must be a non-empty numeric vector.",
      call. = FALSE
    )
  }

  if (any(!is.finite(x))) {
    stop(
      "All observations in 'x' must be finite.",
      call. = FALSE
    )
  }

  if (any(
    x < 0 |
      x > size |
      x != floor(x)
  )) {
    stop(
      "All observations in 'x' must belong to {0, ..., size}.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# -------------------------------------------------------------------------
# Frequency representation
# -------------------------------------------------------------------------

.LUB_frequency_data <- function(x, size) {

  .check_LUB_sample(x, size)

  freq <- tabulate(
    as.integer(x) + 1L,
    nbins = size + 1L
  )

  keep <- freq > 0L

  list(
    support = (0:size)[keep],
    freq = as.numeric(freq[keep]),
    full_freq = as.numeric(freq),
    N = length(x)
  )
}


# -------------------------------------------------------------------------
# Log-likelihood
# -------------------------------------------------------------------------

.loglik_LUB_beta <- function(
    beta,
    x,
    size,
    tol = 1e-14,
    maxit = 10000L) {

  .check_LUB_sample(x, size)
  .check_LUB_control(tol, maxit)

  if (length(beta) != 1L ||
      !is.numeric(beta) ||
      !is.finite(beta) ||
      beta <= 0) {
    return(-Inf)
  }

  dat <- .LUB_frequency_data(
    x = x,
    size = size
  )

  lp <- dLUB(
    x = dat$support,
    size = size,
    beta = beta,
    log = TRUE,
    tol = tol,
    maxit = maxit
  )

  if (any(!is.finite(lp)) ||
      any(lp > 100 * .Machine$double.eps)) {
    return(-Inf)
  }

  sum(dat$freq * lp)
}


# -------------------------------------------------------------------------
# Negative log-likelihood on eta = log(beta)
# -------------------------------------------------------------------------

.nll_LUB <- function(
    eta,
    x,
    size,
    tol = 1e-14,
    maxit = 10000L) {

  if (length(eta) != 1L ||
      !is.numeric(eta) ||
      !is.finite(eta)) {
    return(.LUB_penalty)
  }

  beta <- .beta_from_eta(eta)

  if (!is.finite(beta) ||
      beta <= 0) {
    return(.LUB_penalty)
  }

  ll <- .loglik_LUB_beta(
    beta = beta,
    x = x,
    size = size,
    tol = tol,
    maxit = maxit
  )

  if (!is.finite(ll)) {
    return(.LUB_penalty)
  }

  -ll
}


# -------------------------------------------------------------------------
# Adaptive Brent search
# -------------------------------------------------------------------------

.brent_LUB <- function(
    x,
    size,
    eta.interval = c(-8, 8),
    eta.max.abs = 40,
    tol = 1e-14,
    maxit = 10000L) {

  .check_LUB_sample(x, size)
  .check_LUB_control(tol, maxit)

  if (!is.numeric(eta.interval) ||
      length(eta.interval) != 2L ||
      any(!is.finite(eta.interval)) ||
      eta.interval[1] >= eta.interval[2]) {
    stop(
      "'eta.interval' must contain two finite increasing values.",
      call. = FALSE
    )
  }

  .check_positive_scalar(
    eta.max.abs,
    "eta.max.abs"
  )

  interval <- eta.interval

  repeat {

    z <- stats::optimize(
      f = function(eta) {
        .nll_LUB(
          eta = eta,
          x = x,
          size = size,
          tol = tol,
          maxit = maxit
        )
      },
      interval = interval,
      tol = sqrt(.Machine$double.eps)
    )

    width <- diff(interval)
    edge.tol <- 0.02 * width

    near.left <-
      (z$minimum - interval[1]) < edge.tol

    near.right <-
      (interval[2] - z$minimum) < edge.tol

    if (!near.left && !near.right) {
      return(z)
    }

    old <- interval

    if (near.left &&
        interval[1] > -eta.max.abs) {
      interval[1] <- max(
        -eta.max.abs,
        interval[1] - 5
      )
    }

    if (near.right &&
        interval[2] < eta.max.abs) {
      interval[2] <- min(
        eta.max.abs,
        interval[2] + 5
      )
    }

    if (identical(old, interval)) {
      return(z)
    }
  }
}


# -------------------------------------------------------------------------
# Profile-likelihood confidence interval
# -------------------------------------------------------------------------

.profile_ci_LUB <- function(
    beta.hat,
    logLik.hat,
    x,
    size,
    conf.level = 0.95,
    root.tol = 1e-10,
    tol = 1e-14,
    maxit = 10000L,
    max.expand = 60L) {

  if (!is.finite(beta.hat) ||
      beta.hat <= 0 ||
      !is.finite(logLik.hat)) {

    return(
      list(
        lower = NA_real_,
        upper = NA_real_,
        lower.boundary = NA,
        upper.boundary = NA,
        cutoff = NA_real_,
        target = NA_real_,
        ok = FALSE
      )
    )
  }

  qcrit <- stats::qchisq(
    conf.level,
    df = 1
  )

  cutoff <- qcrit / 2
  target <- logLik.hat - cutoff

  g.eta <- function(eta) {

    beta <- exp(eta)

    ll <- .loglik_LUB_beta(
      beta = beta,
      x = x,
      size = size,
      tol = tol,
      maxit = maxit
    )

    if (!is.finite(ll)) {
      return(-Inf)
    }

    ll - target
  }

  eta.hat <- log(beta.hat)


  # -----------------------------------------------------------------------
  # Lower endpoint
  # -----------------------------------------------------------------------

  eta.lower <- eta.hat
  g.lower <- cutoff
  lower.boundary <- FALSE

  for (i in seq_len(max.expand)) {

    eta.cand <- eta.hat - i

    if (eta.cand < -40) {
      eta.cand <- -40
      lower.boundary <- TRUE
    }

    g.cand <- g.eta(eta.cand)

    if (is.finite(g.cand) &&
        g.cand <= 0) {
      eta.lower <- eta.cand
      g.lower <- g.cand
      break
    }

    eta.lower <- eta.cand
    g.lower <- g.cand

    if (eta.cand <= -40) {
      lower.boundary <- TRUE
      break
    }
  }

  if (is.finite(g.lower) &&
      g.lower <= 0 &&
      eta.lower < eta.hat) {

    lower.try <- try(
      stats::uniroot(
        g.eta,
        interval = c(
          eta.lower,
          eta.hat
        ),
        tol = root.tol
      )$root,
      silent = TRUE
    )

    lower <- if (
      inherits(lower.try, "try-error")
    ) {
      NA_real_
    } else {
      exp(lower.try)
    }

  } else {

    lower <- exp(eta.lower)
    lower.boundary <- TRUE
  }


  # -----------------------------------------------------------------------
  # Upper endpoint
  # -----------------------------------------------------------------------

  eta.upper <- eta.hat
  g.upper <- cutoff
  upper.boundary <- FALSE

  for (i in seq_len(max.expand)) {

    eta.cand <- eta.hat + i

    if (eta.cand > 40) {
      eta.cand <- 40
      upper.boundary <- TRUE
    }

    g.cand <- g.eta(eta.cand)

    if (is.finite(g.cand) &&
        g.cand <= 0) {
      eta.upper <- eta.cand
      g.upper <- g.cand
      break
    }

    eta.upper <- eta.cand
    g.upper <- g.cand

    if (eta.cand >= 40) {
      upper.boundary <- TRUE
      break
    }
  }

  if (is.finite(g.upper) &&
      g.upper <= 0 &&
      eta.upper > eta.hat) {

    upper.try <- try(
      stats::uniroot(
        g.eta,
        interval = c(
          eta.hat,
          eta.upper
        ),
        tol = root.tol
      )$root,
      silent = TRUE
    )

    upper <- if (
      inherits(upper.try, "try-error")
    ) {
      NA_real_
    } else {
      exp(upper.try)
    }

  } else {

    upper <- exp(eta.upper)
    upper.boundary <- TRUE
  }

  ok <-
    is.finite(lower) &&
    is.finite(upper) &&
    lower < beta.hat &&
    beta.hat < upper

  list(
    lower = lower,
    upper = upper,
    lower.boundary = lower.boundary,
    upper.boundary = upper.boundary,
    cutoff = cutoff,
    target = target,
    ok = ok
  )
}


# -------------------------------------------------------------------------
# Likelihood profile
# -------------------------------------------------------------------------

#' Profile Likelihood for the LUB Parameter
#'
#' Evaluates the log-likelihood of the Lambert-uniform binomial model
#' over a user-specified grid of positive beta values.
#'
#' @param x Observed bounded counts.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param beta.grid Positive numeric grid of beta values.
#' @param tol Numerical tolerance passed to `dLUB()`.
#' @param maxit Maximum number of hypergeometric-series terms.
#'
#' @return An object of class `"profileLUB"`.
#'
#' @examples
#' x <- c(0, 0, 1, 0, 2, 1, 0, 0)
#' pr <- profileLUB(x, size = 4)
#' pr
#'
#' @export
profileLUB <- function(
    x,
    size,
    beta.grid = exp(
      seq(
        log(0.05),
        log(30),
        length.out = 300L
      )
    ),
    tol = 1e-14,
    maxit = 10000L) {

  x <- as.numeric(x)

  .check_LUB_sample(x, size)
  .check_LUB_control(tol, maxit)

  if (!is.numeric(beta.grid) ||
      length(beta.grid) < 2L ||
      any(!is.finite(beta.grid)) ||
      any(beta.grid <= 0)) {

    stop(
      paste0(
        "'beta.grid' must contain at least two ",
        "strictly positive finite values."
      ),
      call. = FALSE
    )
  }

  dat <- .LUB_frequency_data(
    x = x,
    size = size
  )

  logLik <- vapply(
    beta.grid,
    function(b) {

      lp <- dLUB(
        x = dat$support,
        size = size,
        beta = b,
        log = TRUE,
        tol = tol,
        maxit = maxit
      )

      if (any(!is.finite(lp))) {
        return(-Inf)
      }

      sum(dat$freq * lp)
    },
    numeric(1)
  )

  imax <- which.max(logLik)

  structure(
    list(
      beta = beta.grid,
      logLik = logLik,
      beta.max.grid = beta.grid[imax],
      logLik.max.grid = logLik[imax],
      size = size,
      n = dat$N,
      nobs = dat$N,
      frequencies = dat$full_freq,
      call = match.call()
    ),
    class = "profileLUB"
  )
}


# -------------------------------------------------------------------------
# Maximum-likelihood estimation
# -------------------------------------------------------------------------

#' Maximum Likelihood Estimation for the LUB Distribution
#'
#' Fits the Lambert-uniform binomial distribution by maximum likelihood.
#' Optimization is performed on the unconstrained parameter
#' \eqn{\eta=\log(\beta)}.
#'
#' A one-dimensional adaptive Brent search is used as the primary
#' optimization procedure. Optional multistart BFGS optimization can be
#' used as a numerical diagnostic.
#'
#' @param x Observed bounded counts.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param start Optional positive starting value for beta used by the
#'   diagnostic BFGS fits.
#' @param starts Positive deterministic starting values for diagnostic
#'   multistart BFGS optimization.
#' @param diagnostic_bfgs Logical; if `TRUE`, diagnostic multistart BFGS
#'   fits are also computed.
#' @param conf.level Confidence level for interval estimation.
#' @param profile.ci Logical; if `TRUE`, a profile-likelihood confidence
#'   interval is additionally computed.
#' @param profile.root.tol Root-finding tolerance for the profile-likelihood
#'   confidence interval.
#' @param maxit Maximum number of iterations for diagnostic BFGS.
#' @param reltol Relative tolerance for diagnostic BFGS.
#' @param tol Numerical tolerance passed to `dLUB()`.
#' @param hyp.maxit Maximum number of terms used by the internal
#'   hypergeometric-series evaluation.
#'
#' @details
#' The primary estimator is obtained by minimizing the negative
#' log-likelihood over \eqn{\eta=\log(\beta)} using an adaptive Brent
#' search.
#'
#' Standard errors are obtained from the observed information on the
#' eta scale and transformed to the beta scale using the delta method.
#' The default confidence interval is a transformed Wald interval.
#'
#' If every observation is zero, no finite interior maximum-likelihood
#' estimate exists. In that case, the likelihood is maximized in the
#' limit as \eqn{\beta\to\infty}, and the function returns an error.
#'
#' @return An object of class `"fitLUB_mle"`.
#'
#' @examples
#' x <- c(0, 0, 1, 0, 2, 1, 0, 0)
#' fit <- fitLUB_mle(
#'   x = x,
#'   size = 4,
#'   diagnostic_bfgs = FALSE
#' )
#' fit
#'
#' @export
fitLUB_mle <- function(
    x,
    size,
    start = NULL,
    starts = c(
      0.05,
      0.10,
      0.25,
      0.50,
      1,
      2,
      5,
      10,
      20
    ),
    diagnostic_bfgs = TRUE,
    conf.level = 0.95,
    profile.ci = FALSE,
    profile.root.tol = 1e-10,
    maxit = 5000L,
    reltol = 1e-10,
    tol = 1e-14,
    hyp.maxit = 10000L) {

  x <- as.numeric(x)

  .check_LUB_sample(x, size)
  .check_LUB_control(tol, hyp.maxit)
  .check_logical_scalar(
    diagnostic_bfgs,
    "diagnostic_bfgs"
  )
  .check_logical_scalar(
    profile.ci,
    "profile.ci"
  )

  # -----------------------------------------------------------------------
  # Boundary case: all observations are zero
  # -----------------------------------------------------------------------

  if (all(x == 0)) {
    stop(
      paste0(
        "No finite interior maximum-likelihood estimate exists when all ",
        "observations are zero: the likelihood is maximized in the limit ",
        "as beta tends to infinity."
      ),
      call. = FALSE
    )
  }

  # -----------------------------------------------------------------------
  # Validate control arguments
  # -----------------------------------------------------------------------

  if (!is.numeric(starts) ||
      length(starts) < 1L ||
      any(!is.finite(starts)) ||
      any(starts <= 0)) {

    stop(
      "'starts' must contain strictly positive finite values.",
      call. = FALSE
    )
  }

  if (!is.null(start)) {
    .check_positive_scalar(
      start,
      "start"
    )
  }

  if (length(conf.level) != 1L ||
      !is.numeric(conf.level) ||
      !is.finite(conf.level) ||
      conf.level <= 0 ||
      conf.level >= 1) {

    stop(
      "'conf.level' must be strictly between 0 and 1.",
      call. = FALSE
    )
  }

  .check_positive_scalar(
    profile.root.tol,
    "profile.root.tol"
  )

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

  .check_positive_scalar(
    reltol,
    "reltol"
  )

  beta.starts <- unique(
    c(
      start,
      starts,
      1
    )
  )

  beta.starts <- beta.starts[
    is.finite(beta.starts) &
      beta.starts > 0
  ]

  eta.starts <- log(beta.starts)

  all_fits <- list()
  kfit <- 0L


  # -----------------------------------------------------------------------
  # Optional multistart BFGS diagnostics
  # -----------------------------------------------------------------------

  if (isTRUE(diagnostic_bfgs)) {

    for (i in seq_along(beta.starts)) {

      kfit <- kfit + 1L

      z <- try(
        stats::optim(
          par = eta.starts[i],
          fn = .nll_LUB,
          x = x,
          size = size,
          method = "BFGS",
          control = list(
            maxit = maxit,
            reltol = reltol
          ),
          tol = tol,
          maxit = hyp.maxit
        ),
        silent = TRUE
      )

      if (inherits(z, "try-error")) {

        all_fits[[kfit]] <- list(
          method = "BFGS",
          start = beta.starts[i],
          convergence = NA_integer_,
          value = Inf,
          par = NA_real_,
          message = "optimization error"
        )

      } else {

        all_fits[[kfit]] <- list(
          method = "BFGS",
          start = beta.starts[i],
          convergence = z$convergence,
          value = z$value,
          par = unname(z$par),
          message = if (
            is.null(z$message)
          ) {
            ""
          } else {
            z$message
          }
        )
      }
    }
  }


  # -----------------------------------------------------------------------
  # Primary estimator: adaptive Brent search
  # -----------------------------------------------------------------------

  kfit <- kfit + 1L

  z.brent <- try(
    .brent_LUB(
      x = x,
      size = size,
      tol = tol,
      maxit = hyp.maxit
    ),
    silent = TRUE
  )

  if (inherits(z.brent, "try-error")) {

    all_fits[[kfit]] <- list(
      method = "Brent",
      start = NA_real_,
      convergence = NA_integer_,
      value = Inf,
      par = NA_real_,
      message = "optimization error"
    )

  } else {

    all_fits[[kfit]] <- list(
      method = "Brent",
      start = NA_real_,
      convergence = 0L,
      value = z.brent$objective,
      par = z.brent$minimum,
      message = ""
    )
  }


  # -----------------------------------------------------------------------
  # Validate numerical solutions
  # -----------------------------------------------------------------------

  valid <- vapply(
    all_fits,
    function(z) {

      if (!is.finite(z$value) ||
          is.na(z$convergence) ||
          z$convergence != 0L ||
          length(z$par) != 1L ||
          !is.finite(z$par)) {
        return(FALSE)
      }

      b <- exp(z$par)

      is.finite(b) &&
        b > 0
    },
    logical(1)
  )

  if (!any(valid)) {
    stop(
      "No valid converged solution was obtained.",
      call. = FALSE
    )
  }

  brent.id <- which(
    vapply(
      all_fits,
      function(z) {
        identical(
          z$method,
          "Brent"
        )
      },
      logical(1)
    )
  )

  brent.valid <- brent.id[
    valid[brent.id]
  ]

  if (length(brent.valid) >= 1L) {

    best.id <- brent.valid[1L]

  } else {

    valid.id <- which(valid)

    values <- vapply(
      all_fits[valid.id],
      function(z) z$value,
      numeric(1)
    )

    best.id <- valid.id[
      which.min(values)
    ]
  }

  best <- all_fits[[best.id]]

  eta.hat <- unname(best$par)
  beta.hat <- exp(eta.hat)

  par <- c(
    beta = beta.hat
  )


  # -----------------------------------------------------------------------
  # Final likelihood
  # -----------------------------------------------------------------------

  nll <- .nll_LUB(
    eta = eta.hat,
    x = x,
    size = size,
    tol = tol,
    maxit = hyp.maxit
  )

  if (!is.finite(nll)) {
    stop(
      "The final likelihood evaluation is not finite.",
      call. = FALSE
    )
  }

  logLik <- -nll


  # -----------------------------------------------------------------------
  # Observed Hessian on eta scale
  # -----------------------------------------------------------------------

  H <- try(
    stats::optimHess(
      par = eta.hat,
      fn = .nll_LUB,
      x = x,
      size = size,
      tol = tol,
      maxit = hyp.maxit
    ),
    silent = TRUE
  )

  hessian_finite <- FALSE
  hessian_ok <- FALSE
  eigenvalue <- NA_real_
  condition_number <- NA_real_

  if (!inherits(H, "try-error") &&
      length(H) == 1L &&
      all(is.finite(H))) {

    hessian_finite <- TRUE

    eigenvalue <- as.numeric(
      H[1, 1]
    )

    hessian_ok <-
      is.finite(eigenvalue) &&
      eigenvalue > 0

    if (hessian_ok) {
      condition_number <- 1
    }
  }


  # -----------------------------------------------------------------------
  # Covariance matrix and standard error
  # -----------------------------------------------------------------------

  Veta <- NULL
  Vpar <- NULL

  se.eta <- NA_real_

  se <- c(
    beta = NA_real_
  )

  vcov_ok <- FALSE
  se_ok <- FALSE

  if (hessian_ok) {

    Veta.try <- try(
      solve(H),
      silent = TRUE
    )

    if (!inherits(Veta.try, "try-error") &&
        all(is.finite(Veta.try)) &&
        Veta.try[1, 1] > 0) {

      Veta <- Veta.try

      dimnames(Veta) <- list(
        "eta",
        "eta"
      )

      se.eta <- sqrt(
        Veta[1, 1]
      )

      # beta = exp(eta)
      jac <- beta.hat

      Vpar <- matrix(
        jac^2 * Veta[1, 1],
        nrow = 1L,
        dimnames = list(
          "beta",
          "beta"
        )
      )

      se["beta"] <- sqrt(
        Vpar[1, 1]
      )

      vcov_ok <-
        is.finite(Vpar[1, 1]) &&
        Vpar[1, 1] > 0

      se_ok <-
        is.finite(se["beta"]) &&
        se["beta"] > 0
    }
  }


  # -----------------------------------------------------------------------
  # Transformed Wald confidence interval
  # -----------------------------------------------------------------------

  lower <- c(
    beta = NA_real_
  )

  upper <- c(
    beta = NA_real_
  )

  if (is.finite(se.eta) &&
      se.eta > 0) {

    zcrit <- stats::qnorm(
      1 - (1 - conf.level) / 2
    )

    lower["beta"] <- exp(
      eta.hat -
        zcrit * se.eta
    )

    upper["beta"] <- exp(
      eta.hat +
        zcrit * se.eta
    )
  }


  # -----------------------------------------------------------------------
  # Optional profile-likelihood confidence interval
  # -----------------------------------------------------------------------

  profile.info <- list(
    lower = NA_real_,
    upper = NA_real_,
    lower.boundary = NA,
    upper.boundary = NA,
    cutoff = NA_real_,
    target = NA_real_,
    ok = FALSE
  )

  if (isTRUE(profile.ci)) {

    profile.info <- .profile_ci_LUB(
      beta.hat = beta.hat,
      logLik.hat = logLik,
      x = x,
      size = size,
      conf.level = conf.level,
      root.tol = profile.root.tol,
      tol = tol,
      maxit = hyp.maxit
    )
  }

  lr.lower <- c(
    beta = profile.info$lower
  )

  lr.upper <- c(
    beta = profile.info$upper
  )

  coefficients <- cbind(
    Estimate = par,
    SE = se,
    Lower = lower,
    Upper = upper,
    LR.Lower = lr.lower,
    LR.Upper = lr.upper
  )


  # -----------------------------------------------------------------------
  # Information criteria
  # -----------------------------------------------------------------------

  N <- length(x)
  k <- 1L

  AIC <-
    -2 * logLik +
    2 * k

  AICc <- if (
    N > k + 1L
  ) {

    AIC +
      2 * k * (k + 1) /
      (N - k - 1)

  } else {

    NA_real_
  }

  BIC <-
    -2 * logLik +
    log(N) * k


  # -----------------------------------------------------------------------
  # Discrete-uniform reference: beta = 1
  # -----------------------------------------------------------------------

  logLik.uniform <-
    -N * log(size + 1)

  AIC.uniform <-
    -2 * logLik.uniform

  AICc.uniform <-
    AIC.uniform

  BIC.uniform <-
    -2 * logLik.uniform


  # -----------------------------------------------------------------------
  # Optimization diagnostics
  # -----------------------------------------------------------------------

  fits.summary <- do.call(
    rbind,
    lapply(
      all_fits,
      function(z) {

        bhat <- if (
          length(z$par) == 1L &&
          is.finite(z$par)
        ) {
          exp(z$par)
        } else {
          NA_real_
        }

        data.frame(
          method = z$method,
          start = z$start,
          convergence = z$convergence,
          nll = z$value,
          beta = bhat,
          message = z$message,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  finite.conv <-
    is.finite(fits.summary$nll) &
    !is.na(fits.summary$convergence) &
    fits.summary$convergence == 0L &
    is.finite(fits.summary$beta)

  best.nll.observed <- min(
    fits.summary$nll[
      finite.conv
    ]
  )

  nll.tolerance <-
    1e-7 *
    (
      1 +
        abs(best.nll.observed)
    )

  fits.summary$nll_gap <-
    fits.summary$nll -
    best.nll.observed

  fits.summary$near_optimal <-
    finite.conv &
    fits.summary$nll_gap <=
    nll.tolerance

  near.beta <-
    fits.summary$beta[
      fits.summary$near_optimal
    ]

  solution_spread <- if (
    length(near.beta) > 1L
  ) {
    max(near.beta) -
      min(near.beta)
  } else {
    0
  }


  # -----------------------------------------------------------------------
  # Output
  # -----------------------------------------------------------------------

  dat <- .LUB_frequency_data(
    x = x,
    size = size
  )

  out <- list(
    call = match.call(),
    par = par,
    se = se,
    coefficients = coefficients,

    eta = eta.hat,
    se.eta = se.eta,

    logLik = logLik,
    nll = nll,

    AIC = AIC,
    AICc = AICc,
    BIC = BIC,

    n = N,
    nobs = N,
    size = size,
    npar = k,

    method = best$method,
    convergence = best$convergence,
    message = best$message,
    success = isTRUE(
      best$convergence == 0L
    ),

    hessian = if (
      inherits(H, "try-error")
    ) {
      NULL
    } else {
      H
    },

    hessian_finite = hessian_finite,
    hessian_ok = hessian_ok,
    eigenvalue = eigenvalue,
    condition_number = condition_number,

    vcov.eta = Veta,
    vcov = Vpar,
    vcov_ok = vcov_ok,
    se_ok = se_ok,

    conf.level = conf.level,

    profile.ci = isTRUE(profile.ci),

    profile_ci = c(
      lower = profile.info$lower,
      upper = profile.info$upper
    ),

    profile_ci_ok =
      isTRUE(profile.info$ok),

    profile_ci_boundary = c(
      lower =
        isTRUE(profile.info$lower.boundary),
      upper =
        isTRUE(profile.info$upper.boundary)
    ),

    profile_lr_cutoff =
      profile.info$cutoff,

    profile_lr_target =
      profile.info$target,

    best.start = best$start,
    all_fits = fits.summary,
    solution_spread = solution_spread,
    nll_tolerance = nll.tolerance,
    diagnostic_bfgs =
      diagnostic_bfgs,

    frequencies =
      dat$full_freq,

    data = x,

    uniform = list(
      beta = 1,
      logLik = logLik.uniform,
      AIC = AIC.uniform,
      AICc = AICc.uniform,
      BIC = BIC.uniform
    )
  )

  class(out) <- "fitLUB_mle"

  out
}


# -------------------------------------------------------------------------
# Print method
# -------------------------------------------------------------------------

#' @export
print.fitLUB_mle <- function(
    x,
    digits = 4,
    ...) {

  cat(
    "Lambert-uniform binomial model fitted by maximum likelihood\n",
    "----------------------------------------------------------\n",
    sep = ""
  )

  cat(
    "Observations :",
    x$n,
    "\n"
  )

  cat(
    "Size         :",
    x$size,
    "\n\n"
  )

  print(
    round(
      x$coefficients,
      digits
    )
  )

  cat(
    "\nLog-likelihood :",
    format(
      x$logLik,
      digits = digits + 2
    ),
    "\n"
  )

  cat(
    "AIC            :",
    format(
      x$AIC,
      digits = digits + 2
    ),
    "\n"
  )

  cat(
    "AICc           :",
    format(
      x$AICc,
      digits = digits + 2
    ),
    "\n"
  )

  cat(
    "BIC            :",
    format(
      x$BIC,
      digits = digits + 2
    ),
    "\n"
  )

  cat(
    "Estimator      :",
    x$method,
    "\n"
  )

  cat(
    "Convergence    :",
    x$convergence,
    "\n"
  )

  cat(
    "Hessian OK     :",
    x$hessian_ok,
    "\n"
  )

  if (isTRUE(x$profile.ci)) {

    cat(
      "Profile CI     :",
      if (
        isTRUE(x$profile_ci_ok)
      ) {
        "available"
      } else {
        "boundary/diagnostic issue"
      },
      "\n"
    )
  }

  invisible(x)
}


# -------------------------------------------------------------------------
# Summary method
# -------------------------------------------------------------------------

#' @export
summary.fitLUB_mle <- function(
    object,
    ...) {

  out <- list(
    call = object$call,
    coefficients = object$coefficients,

    logLik = object$logLik,
    AIC = object$AIC,
    AICc = object$AICc,
    BIC = object$BIC,

    n = object$n,
    nobs = object$nobs,
    size = object$size,

    method = object$method,
    convergence = object$convergence,

    hessian_ok =
      object$hessian_ok,

    profile_ci =
      object$profile_ci,

    profile_ci_ok =
      object$profile_ci_ok,

    profile_ci_boundary =
      object$profile_ci_boundary
  )

  class(out) <-
    "summary.fitLUB_mle"

  out
}


# -------------------------------------------------------------------------
# Print summary
# -------------------------------------------------------------------------

#' @export
print.summary.fitLUB_mle <- function(
    x,
    digits = 4,
    ...) {

  cat(
    "Summary of Lambert-uniform binomial fit\n",
    "---------------------------------------\n",
    sep = ""
  )

  cat(
    "Observations :",
    x$n,
    "\n"
  )

  cat(
    "Size         :",
    x$size,
    "\n"
  )

  cat(
    "Estimator    :",
    x$method,
    "\n\n"
  )

  print(
    round(
      x$coefficients,
      digits
    )
  )

  cat(
    "\nLog-likelihood :",
    format(
      x$logLik,
      digits = digits + 2
    ),
    "\n"
  )

  cat(
    "AIC            :",
    format(
      x$AIC,
      digits = digits + 2
    ),
    "\n"
  )

  cat(
    "AICc           :",
    format(
      x$AICc,
      digits = digits + 2
    ),
    "\n"
  )

  cat(
    "BIC            :",
    format(
      x$BIC,
      digits = digits + 2
    ),
    "\n"
  )

  cat(
    "Convergence    :",
    x$convergence,
    "\n"
  )

  cat(
    "Hessian OK     :",
    x$hessian_ok,
    "\n"
  )

  invisible(x)
}


# -------------------------------------------------------------------------
# Print likelihood profile
# -------------------------------------------------------------------------

#' @export
print.profileLUB <- function(
    x,
    digits = 5,
    ...) {

  cat(
    "LUB log-likelihood profile\n",
    "--------------------------\n",
    sep = ""
  )

  cat(
    "Observations :",
    x$n,
    "\n"
  )

  cat(
    "Size         :",
    x$size,
    "\n"
  )

  cat(
    "Grid maximum : beta =",
    format(
      x$beta.max.grid,
      digits = digits
    ),
    "\n"
  )

  cat(
    "LogLik       :",
    format(
      x$logLik.max.grid,
      digits = digits
    ),
    "\n"
  )

  invisible(x)
}


# -------------------------------------------------------------------------
# Plot likelihood profile
# -------------------------------------------------------------------------

#' @export
plot.profileLUB <- function(
    x,
    ...,
    xlab = expression(beta),
    ylab = "Log-likelihood") {

  plot(
    x$beta,
    x$logLik,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    ...
  )

  invisible(x)
}
