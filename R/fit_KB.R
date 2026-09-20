# -------------------------------------------------------------------------
# Maximum-likelihood estimation for the
# Kumaraswamy-binomial distribution
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Sample validation
# -------------------------------------------------------------------------

.check_KB_sample <- function(
    x,
    size) {

  .check_KB(
    size = size,
    a = 1,
    b = 1
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
# Negative log-likelihood
# -------------------------------------------------------------------------

.KB_nll_eta <- function(
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

  a <- exp(
    eta[1L]
  )

  b <- exp(
    eta[2L]
  )

  if (!is.finite(a) ||
      !is.finite(b) ||
      a <= 0 ||
      b <= 0) {
    return(
      .Machine$double.xmax^0.25
    )
  }

  lp <- try(
    dKB(
      x = x,
      size = size,
      a = a,
      b = b,
      log = TRUE
    ),
    silent = TRUE
  )

  if (inherits(
    lp,
    "try-error"
  ) ||
      any(!is.finite(lp))) {
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

#' Maximum-Likelihood Estimation for the Kumaraswamy-Binomial Distribution
#'
#' Fits the Kumaraswamy-binomial distribution to bounded count data using
#' multistart maximum-likelihood optimization.
#'
#' @param x Observed bounded counts.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param starts Optional data frame or matrix containing positive starting
#'   values for `a` and `b`.
#' @param conf.level Confidence level for transformed Wald intervals.
#' @param maxit Maximum number of BFGS iterations.
#' @param reltol Relative convergence tolerance for BFGS.
#'
#' @details
#' Optimization is performed over
#' \eqn{(\log a,\log b)}, which automatically enforces positivity.
#' Multiple deterministic starting points are considered and the converged
#' solution with the smallest negative log-likelihood is retained.
#'
#' Standard errors are obtained from the numerical Hessian on the
#' log-parameter scale and transformed by the delta method.
#'
#' @return
#' A list containing parameter estimates, standard errors, confidence
#' intervals, log-likelihood, information criteria, Hessian diagnostics,
#' optimization diagnostics, and model metadata.
#'
#' @examples
#' x <- c(0, 0, 1, 0, 2, 1, 0, 0)
#' fitKB_mle(
#'   x = x,
#'   size = 4
#' )
#'
#' @export
fitKB_mle <- function(
    x,
    size,
    starts = NULL,
    conf.level = 0.95,
    maxit = 5000L,
    reltol = 1e-10) {

  x <- as.numeric(x)

  .check_KB_sample(
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

    starts <- expand.grid(
      a = c(
        0.25,
        0.50,
        1,
        2,
        5,
        10
      ),
      b = c(
        0.25,
        0.50,
        1,
        2,
        5,
        10
      ),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }

  starts <- as.data.frame(
    starts
  )

  if (!all(
    c("a", "b") %in%
      names(starts)
  )) {
    stop(
      "'starts' must contain columns named 'a' and 'b'.",
      call. = FALSE
    )
  }

  starts <- starts[
    is.finite(starts$a) &
      starts$a > 0 &
      is.finite(starts$b) &
      starts$b > 0,
    ,
    drop = FALSE
  ]

  if (nrow(starts) == 0L) {
    stop(
      "No valid starting values were supplied.",
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

  for (i in seq_len(
    nrow(starts)
  )) {

    eta0 <- log(
      c(
        starts$a[i],
        starts$b[i]
      )
    )

    z <- try(
      stats::optim(
        par = eta0,
        fn = .KB_nll_eta,
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

    if (inherits(
      z,
      "try-error"
    )) {

      fits[[i]] <- list(
        value = Inf,
        convergence = NA_integer_,
        par = c(
          NA_real_,
          NA_real_
        ),
        start.a = starts$a[i],
        start.b = starts$b[i],
        message = "optimization error"
      )

    } else {

      fits[[i]] <- list(
        value = z$value,
        convergence = z$convergence,
        par = unname(z$par),
        start.a = starts$a[i],
        start.b = starts$b[i],
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
  # Select best solution
  # -----------------------------------------------------------------------

  valid <- vapply(
    fits,
    function(z) {

      !is.na(z$convergence) &&
        z$convergence == 0L &&
        is.finite(z$value) &&
        length(z$par) == 2L &&
        all(is.finite(z$par))
    },
    logical(1)
  )

  if (!any(valid)) {
    stop(
      "No valid converged KB solution was obtained.",
      call. = FALSE
    )
  }

  ids <- which(valid)

  best.id <- ids[
    which.min(
      vapply(
        fits[ids],
        function(z) z$value,
        numeric(1)
      )
    )
  ]

  best <- fits[[best.id]]

  eta.hat <- best$par

  par <- c(
    a = exp(eta.hat[1L]),
    b = exp(eta.hat[2L])
  )

  nll <- .KB_nll_eta(
    eta = eta.hat,
    x = x,
    size = size
  )

  logLik <- -nll


  # -----------------------------------------------------------------------
  # Observed Hessian
  # -----------------------------------------------------------------------

  H <- try(
    stats::optimHess(
      par = eta.hat,
      fn = .KB_nll_eta,
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

  Veta <- matrix(
    NA_real_,
    2L,
    2L,
    dimnames = list(
      c("log_a", "log_b"),
      c("log_a", "log_b")
    )
  )

  Vpar <- matrix(
    NA_real_,
    2L,
    2L,
    dimnames = list(
      names(par),
      names(par)
    )
  )

  se <- c(
    a = NA_real_,
    b = NA_real_
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

          Veta <- Veta.try

          dimnames(Veta) <- list(
            c("log_a", "log_b"),
            c("log_a", "log_b")
          )

          J <- diag(
            par,
            nrow = 2L
          )

          Vpar <-
            J %*%
            Veta %*%
            J

          dimnames(Vpar) <- list(
            names(par),
            names(par)
          )

          vv <- diag(Vpar)

          if (all(
            is.finite(vv)
          ) &&
              all(vv > 0)) {

            se <- sqrt(vv)

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
    a = NA_real_,
    b = NA_real_
  )

  upper <- c(
    a = NA_real_,
    b = NA_real_
  )

  if (hessian_ok) {

    se.eta <- sqrt(
      diag(Veta)
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
    log(N) * k


  # -----------------------------------------------------------------------
  # Optimization diagnostics
  # -----------------------------------------------------------------------

  all_fits <- do.call(
    rbind,
    lapply(
      fits,
      function(z) {

        a.est <- if (
          length(z$par) == 2L &&
          all(is.finite(z$par))
        ) {
          exp(z$par[1L])
        } else {
          NA_real_
        }

        b.est <- if (
          length(z$par) == 2L &&
          all(is.finite(z$par))
        ) {
          exp(z$par[2L])
        } else {
          NA_real_
        }

        data.frame(
          start.a = z$start.a,
          start.b = z$start.b,
          convergence = z$convergence,
          nll = z$value,
          a = a.est,
          b = b.est,
          message = z$message,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  finite.conv <-
    !is.na(
      all_fits$convergence
    ) &
    all_fits$convergence == 0L &
    is.finite(
      all_fits$nll
    ) &
    is.finite(
      all_fits$a
    ) &
    is.finite(
      all_fits$b
    )

  best.nll.observed <- min(
    all_fits$nll[
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

  all_fits$nll_gap <-
    all_fits$nll -
    best.nll.observed

  all_fits$near_optimal <-
    finite.conv &
    all_fits$nll_gap <=
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

    eta = eta.hat,

    logLik = logLik,
    nll = nll,

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

    vcov.eta = Veta,
    vcov = Vpar,

    conf.level =
      conf.level,

    best.start = c(
      a = best$start.a,
      b = best$start.b
    ),

    all_fits =
      all_fits,

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
