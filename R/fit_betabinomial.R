# -------------------------------------------------------------------------
# Maximum-likelihood estimation for the beta-binomial distribution
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Sample validation
# -------------------------------------------------------------------------

.check_BB_sample <- function(
    x,
    size) {

  .check_BB_parameters(
    size = size,
    alpha = 1,
    beta = 1
  )

  if (!is.numeric(x) ||
      length(x) < 1L) {
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
# Negative log-likelihood on the log-parameter scale
# -------------------------------------------------------------------------

.nll_BB <- function(
    eta,
    x,
    size) {

  if (length(eta) != 2L ||
      !is.numeric(eta) ||
      any(!is.finite(eta))) {
    return(
      .Machine$double.xmax^0.25
    )
  }

  alpha <- exp(eta[1L])
  beta <- exp(eta[2L])

  if (!is.finite(alpha) ||
      !is.finite(beta) ||
      alpha <= 0 ||
      beta <= 0) {
    return(
      .Machine$double.xmax^0.25
    )
  }

  lp <- dBB(
    x = x,
    size = size,
    alpha = alpha,
    beta = beta,
    log = TRUE
  )

  if (any(!is.finite(lp))) {
    return(
      .Machine$double.xmax^0.25
    )
  }

  value <- -sum(lp)

  if (!is.finite(value)) {
    return(
      .Machine$double.xmax^0.25
    )
  }

  value
}


# -------------------------------------------------------------------------
# Maximum-likelihood estimation
# -------------------------------------------------------------------------

#' Maximum-Likelihood Estimation for the Beta-Binomial Distribution
#'
#' Fits the beta-binomial distribution to bounded count data using
#' multistart maximum-likelihood optimization.
#'
#' @param x Observed bounded counts.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param starts Optional two-column matrix containing positive starting
#'   values for `alpha` and `beta`.
#' @param conf.level Confidence level for transformed Wald intervals.
#' @param maxit Maximum number of BFGS iterations.
#' @param reltol Relative convergence tolerance for BFGS.
#'
#' @details
#' Optimization is performed over
#' \eqn{(\log\alpha,\log\beta)}, ensuring positivity of both parameters.
#' Multiple deterministic starting points are used and the converged
#' solution with the largest likelihood is retained.
#'
#' Standard errors are obtained from the observed Hessian on the
#' log-parameter scale and transformed to the original parameter scale by
#' the delta method. Confidence limits are obtained by transforming Wald
#' intervals constructed on the log scale.
#'
#' @return
#' A list containing parameter estimates, standard errors, confidence
#' limits, log-likelihood, information criteria, Hessian diagnostics,
#' optimization diagnostics, and model metadata.
#'
#' @examples
#' x <- c(0, 0, 1, 0, 2, 1, 0, 0)
#' fitBB_mle(
#'   x = x,
#'   size = 4
#' )
#'
#' @export
fitBB_mle <- function(
    x,
    size,
    starts = NULL,
    conf.level = 0.95,
    maxit = 10000L,
    reltol = 1e-12) {

  x <- as.numeric(x)

  .check_BB_sample(
    x = x,
    size = size
  )


  # -----------------------------------------------------------------------
  # Control arguments
  # -----------------------------------------------------------------------

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


  # -----------------------------------------------------------------------
  # Starting values
  # -----------------------------------------------------------------------

  if (is.null(starts)) {

    p0 <- min(
      max(
        mean(x) / size,
        1e-4
      ),
      1 - 1e-4
    )

    bin.var <-
      size *
      p0 *
      (1 - p0)

    rho0 <- if (
      length(x) > 1L &&
      is.finite(stats::var(x)) &&
      bin.var > 0 &&
      size > 1
    ) {

      (
        stats::var(x) /
          bin.var -
          1
      ) /
        (size - 1)

    } else {

      0.05
    }

    rho0 <- min(
      max(
        rho0,
        1e-3
      ),
      0.95
    )

    phi0 <-
      1 / rho0 -
      1

    starts <- rbind(
      c(
        p0 * phi0,
        (1 - p0) * phi0
      ),
      c(0.10, 0.10),
      c(0.25, 0.25),
      c(0.50, 0.50),
      c(1, 1),
      c(2, 2),
      c(5, 5),
      c(0.25, 1),
      c(1, 0.25),
      c(0.50, 2),
      c(2, 0.50),
      c(1, 5),
      c(5, 1)
    )
  }

  starts <- as.matrix(starts)

  if (!is.numeric(starts) ||
      nrow(starts) < 1L ||
      ncol(starts) != 2L ||
      any(!is.finite(starts)) ||
      any(starts <= 0)) {
    stop(
      paste0(
        "'starts' must be a numeric two-column matrix ",
        "of strictly positive finite values."
      ),
      call. = FALSE
    )
  }


  # -----------------------------------------------------------------------
  # Multistart BFGS
  # -----------------------------------------------------------------------

  fits <- vector(
    "list",
    nrow(starts)
  )

  for (i in seq_len(nrow(starts))) {

    st <- starts[i, ]

    z <- try(
      stats::optim(
        par = log(st),
        fn = .nll_BB,
        x = x,
        size = size,
        method = "BFGS",
        control = list(
          maxit = maxit,
          reltol = reltol
        )
      ),
      silent = TRUE
    )

    if (inherits(z, "try-error")) {

      fits[[i]] <- list(
        start.alpha = st[1L],
        start.beta = st[2L],
        convergence = NA_integer_,
        nll = Inf,
        eta = c(
          NA_real_,
          NA_real_
        ),
        alpha = NA_real_,
        beta = NA_real_,
        message = "optimization error"
      )

    } else {

      fits[[i]] <- list(
        start.alpha = st[1L],
        start.beta = st[2L],
        convergence = z$convergence,
        nll = z$value,
        eta = unname(z$par),
        alpha = exp(z$par[1L]),
        beta = exp(z$par[2L]),
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


  # -----------------------------------------------------------------------
  # Select best converged solution
  # -----------------------------------------------------------------------

  valid <- vapply(
    fits,
    function(z) {

      !is.na(z$convergence) &&
        z$convergence == 0L &&
        is.finite(z$nll) &&
        is.finite(z$alpha) &&
        z$alpha > 0 &&
        is.finite(z$beta) &&
        z$beta > 0
    },
    logical(1)
  )

  if (!any(valid)) {
    stop(
      "No converged beta-binomial fit was obtained.",
      call. = FALSE
    )
  }

  nlls <- vapply(
    fits,
    function(z) z$nll,
    numeric(1)
  )

  best.id <- which.min(
    ifelse(
      valid,
      nlls,
      Inf
    )
  )

  best <- fits[[best.id]]

  eta.hat <- best$eta

  par <- c(
    alpha = best$alpha,
    beta = best$beta
  )

  logLik <- -best$nll


  # -----------------------------------------------------------------------
  # Observed Hessian
  # -----------------------------------------------------------------------

  H <- try(
    stats::optimHess(
      par = eta.hat,
      fn = .nll_BB,
      x = x,
      size = size
    ),
    silent = TRUE
  )

  hessian_finite <- FALSE
  hessian_ok <- FALSE

  eigenvalues <- rep(
    NA_real_,
    2L
  )

  condition_number <- NA_real_

  vcov.eta <- matrix(
    NA_real_,
    2L,
    2L,
    dimnames = list(
      c("log_alpha", "log_beta"),
      c("log_alpha", "log_beta")
    )
  )

  vcov.par <- matrix(
    NA_real_,
    2L,
    2L,
    dimnames = list(
      names(par),
      names(par)
    )
  )

  se <- c(
    alpha = NA_real_,
    beta = NA_real_
  )

  if (!inherits(H, "try-error") &&
      all(dim(H) == c(2L, 2L)) &&
      all(is.finite(H))) {

    hessian_finite <- TRUE

    eig.try <- try(
      eigen(
        H,
        symmetric = TRUE,
        only.values = TRUE
      )$values,
      silent = TRUE
    )

    if (!inherits(
      eig.try,
      "try-error"
    ) &&
        all(is.finite(eig.try))) {

      eigenvalues <- eig.try

      hessian_ok <-
        all(eigenvalues > 0)

      if (hessian_ok) {

        condition_number <-
          max(eigenvalues) /
          min(eigenvalues)

        Veta.try <- try(
          solve(H),
          silent = TRUE
        )

        if (!inherits(
          Veta.try,
          "try-error"
        ) &&
            all(is.finite(Veta.try))) {

          vcov.eta <- Veta.try

          dimnames(vcov.eta) <- list(
            c(
              "log_alpha",
              "log_beta"
            ),
            c(
              "log_alpha",
              "log_beta"
            )
          )

          G <- diag(
            par,
            nrow = 2L
          )

          vcov.par <-
            G %*%
            vcov.eta %*%
            G

          dimnames(vcov.par) <- list(
            names(par),
            names(par)
          )

          var.diag <- diag(
            vcov.par
          )

          if (all(
            is.finite(var.diag)
          ) &&
              all(var.diag > 0)) {

            se <- sqrt(
              var.diag
            )

          } else {

            hessian_ok <- FALSE
          }
        } else {

          hessian_ok <- FALSE
        }
      }
    }
  }


  # -----------------------------------------------------------------------
  # Transformed Wald confidence intervals
  # -----------------------------------------------------------------------

  lower <- c(
    alpha = NA_real_,
    beta = NA_real_
  )

  upper <- c(
    alpha = NA_real_,
    beta = NA_real_
  )

  if (hessian_ok) {

    se.eta <- sqrt(
      diag(vcov.eta)
    )

    if (all(
      is.finite(se.eta)
    ) &&
        all(se.eta > 0)) {

      zcrit <- stats::qnorm(
        1 -
          (1 - conf.level) / 2
      )

      lower[] <- exp(
        eta.hat -
          zcrit * se.eta
      )

      upper[] <- exp(
        eta.hat +
          zcrit * se.eta
      )
    }
  }


  # -----------------------------------------------------------------------
  # Information criteria
  # -----------------------------------------------------------------------

  N <- length(x)
  k <- 2L

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
    k * log(N)


  # -----------------------------------------------------------------------
  # Optimization diagnostics
  # -----------------------------------------------------------------------

  all.fits <- do.call(
    rbind,
    lapply(
      fits,
      function(z) {

        data.frame(
          start.alpha =
            z$start.alpha,
          start.beta =
            z$start.beta,
          convergence =
            z$convergence,
          nll =
            z$nll,
          alpha =
            z$alpha,
          beta =
            z$beta,
          message =
            z$message,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  finite.conv <-
    !is.na(
      all.fits$convergence
    ) &
    all.fits$convergence == 0L &
    is.finite(
      all.fits$nll
    ) &
    is.finite(
      all.fits$alpha
    ) &
    is.finite(
      all.fits$beta
    )

  best.nll.observed <- min(
    all.fits$nll[
      finite.conv
    ]
  )

  nll.tolerance <-
    1e-7 *
    (
      1 +
        abs(
          best.nll.observed
        )
    )

  all.fits$nll_gap <-
    all.fits$nll -
    best.nll.observed

  all.fits$near_optimal <-
    finite.conv &
    all.fits$nll_gap <=
    nll.tolerance


  # -----------------------------------------------------------------------
  # Output
  # -----------------------------------------------------------------------

  out <- list(
    call = match.call(),

    par = par,
    se = se,

    coefficients = cbind(
      Estimate = par,
      SE = se,
      Lower = lower,
      Upper = upper
    ),

    logLik = logLik,
    nll = best$nll,

    AIC = AIC,
    AICc = AICc,
    BIC = BIC,

    convergence =
      best$convergence,

    success =
      isTRUE(
        best$convergence == 0L
      ),

    hessian = if (
      inherits(H, "try-error")
    ) {
      NULL
    } else {
      H
    },

    hessian_finite =
      hessian_finite,

    hessian_ok =
      hessian_ok,

    eigenvalues =
      eigenvalues,

    condition_number =
      condition_number,

    vcov.eta =
      vcov.eta,

    vcov =
      vcov.par,

    conf.level =
      conf.level,

    best.start = c(
      alpha =
        best$start.alpha,
      beta =
        best$start.beta
    ),

    all_fits =
      all.fits,

    nll_tolerance =
      nll.tolerance,

    size = size,
    n = N,
    nobs = N,
    npar = k,

    data = x
  )

  out
}
