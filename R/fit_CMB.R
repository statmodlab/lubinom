# -------------------------------------------------------------------------
# Maximum-likelihood estimation for the
# Conway-Maxwell-binomial distribution
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Sample validation
# -------------------------------------------------------------------------

.check_CMB_sample <- function(
    x,
    size) {

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
#
# eta[1] = logit(prob)
# eta[2] = log(nu)
# -------------------------------------------------------------------------

.CMB_nll_eta <- function(
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

  prob <- stats::plogis(
    eta[1L]
  )

  nu <- exp(
    eta[2L]
  )

  if (!is.finite(prob) ||
      prob <= 0 ||
      prob >= 1 ||
      !is.finite(nu) ||
      nu <= 0) {
    return(
      .Machine$double.xmax^0.25
    )
  }

  lp <- dCMB(
    x = x,
    size = size,
    prob = prob,
    nu = nu,
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

#' Maximum-Likelihood Estimation for the Conway-Maxwell-Binomial Distribution
#'
#' Fits the Conway-Maxwell-binomial distribution to bounded count data
#' using multistart maximum-likelihood optimization.
#'
#' @param x Observed bounded counts.
#' @param size Positive integer giving the upper endpoint of the support.
#' @param starts Optional data frame or matrix containing starting values
#'   for `prob` and `nu`.
#' @param conf.level Confidence level for transformed Wald intervals.
#' @param maxit Maximum number of BFGS iterations.
#' @param reltol Relative convergence tolerance for BFGS.
#'
#' @details
#' Optimization is performed using
#' \deqn{
#' \eta_1=\mathrm{logit}(prob),
#' \qquad
#' \eta_2=\log(\nu),
#' }
#' so that the parameter constraints \eqn{0<prob<1} and \eqn{\nu>0}
#' are automatically satisfied.
#'
#' Multiple deterministic starting values are considered and the
#' converged solution with the smallest negative log-likelihood is retained.
#'
#' Standard errors are obtained from the numerical Hessian on the
#' transformed scale and mapped to the original parameter scale by the
#' delta method. Confidence intervals are constructed on the transformed
#' scales and then mapped back to the original parameter space.
#'
#' @return
#' A list containing parameter estimates, standard errors, confidence
#' intervals, log-likelihood, information criteria, Hessian diagnostics,
#' optimization diagnostics, and model metadata.
#'
#' @examples
#' x <- c(0, 0, 1, 0, 2, 1, 0, 0)
#' fitCMB_mle(
#'   x = x,
#'   size = 4
#' )
#'
#' @export
fitCMB_mle <- function(
    x,
    size,
    starts = NULL,
    conf.level = 0.95,
    maxit = 5000L,
    reltol = 1e-10) {

  x <- as.numeric(x)

  .check_CMB_sample(
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
        0.02
      ),
      0.98
    )

    prob.starts <- unique(
      pmin(
        pmax(
          c(
            p0,
            0.5 * p0,
            1.5 * p0,
            0.10,
            0.25,
            0.50,
            0.75,
            0.90
          ),
          0.01
        ),
        0.99
      )
    )

    starts <- expand.grid(
      prob = prob.starts,
      nu = c(
        0.25,
        0.50,
        1,
        2,
        5
      ),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }

  starts <- as.data.frame(
    starts
  )

  if (!all(
    c("prob", "nu") %in%
      names(starts)
  )) {
    stop(
      "'starts' must contain columns named 'prob' and 'nu'.",
      call. = FALSE
    )
  }

  starts <- starts[
    is.finite(starts$prob) &
      starts$prob > 0 &
      starts$prob < 1 &
      is.finite(starts$nu) &
      starts$nu > 0,
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

    eta0 <- c(
      stats::qlogis(
        starts$prob[i]
      ),
      log(
        starts$nu[i]
      )
    )

    z <- try(
      stats::optim(
        par = eta0,
        fn = .CMB_nll_eta,
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
        start.prob =
          starts$prob[i],
        start.nu =
          starts$nu[i],
        message =
          "optimization error"
      )

    } else {

      fits[[i]] <- list(
        value = z$value,
        convergence =
          z$convergence,
        par = unname(z$par),
        start.prob =
          starts$prob[i],
        start.nu =
          starts$nu[i],
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
        is.finite(z$value) &&
        length(z$par) == 2L &&
        all(is.finite(z$par))
    },
    logical(1)
  )

  if (!any(valid)) {
    stop(
      "No valid converged CMB solution was obtained.",
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

  prob.hat <- stats::plogis(
    eta.hat[1L]
  )

  nu.hat <- exp(
    eta.hat[2L]
  )

  par <- c(
    prob = prob.hat,
    nu = nu.hat
  )

  nll <- .CMB_nll_eta(
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
      fn = .CMB_nll_eta,
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
      c(
        "logit_prob",
        "log_nu"
      ),
      c(
        "logit_prob",
        "log_nu"
      )
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
    prob = NA_real_,
    nu = NA_real_
  )

  if (!inherits(
    H,
    "try-error"
  ) &&
      all(
        dim(H) ==
          c(2L, 2L)
      ) &&
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
        all(
          is.finite(
            eig.try
          )
        )) {

      eigenvalues <- eig.try

      hessian_ok <-
        all(
          eigenvalues > 0
        )

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
            all(
              is.finite(
                Veta.try
              )
            )) {

          Veta <- Veta.try

          dimnames(Veta) <- list(
            c(
              "logit_prob",
              "log_nu"
            ),
            c(
              "logit_prob",
              "log_nu"
            )
          )

          J <- diag(
            c(
              prob.hat *
                (1 - prob.hat),
              nu.hat
            ),
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
    prob = NA_real_,
    nu = NA_real_
  )

  upper <- c(
    prob = NA_real_,
    nu = NA_real_
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

      lower["prob"] <-
        stats::plogis(
          eta.hat[1L] -
            zcrit *
            se.eta[1L]
        )

      upper["prob"] <-
        stats::plogis(
          eta.hat[1L] +
            zcrit *
            se.eta[1L]
        )

      lower["nu"] <-
        exp(
          eta.hat[2L] -
            zcrit *
            se.eta[2L]
        )

      upper["nu"] <-
        exp(
          eta.hat[2L] +
            zcrit *
            se.eta[2L]
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

        prob.est <- if (
          length(z$par) == 2L &&
          all(
            is.finite(
              z$par
            )
          )
        ) {
          stats::plogis(
            z$par[1L]
          )
        } else {
          NA_real_
        }

        nu.est <- if (
          length(z$par) == 2L &&
          all(
            is.finite(
              z$par
            )
          )
        ) {
          exp(
            z$par[2L]
          )
        } else {
          NA_real_
        }

        data.frame(
          start.prob =
            z$start.prob,
          start.nu =
            z$start.nu,
          convergence =
            z$convergence,
          nll =
            z$value,
          prob =
            prob.est,
          nu =
            nu.est,
          message =
            z$message,
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
      all_fits$prob
    ) &
    is.finite(
      all_fits$nu
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
      prob =
        best$start.prob,
      nu =
        best$start.nu
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
