# -------------------------------------------------------------------------
# Comparison of bounded-count models
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Internal helper: safely fit a model
# -------------------------------------------------------------------------

.compare_fit_try <- function(
    expr,
    suppress_warnings = TRUE) {

  if (isTRUE(suppress_warnings)) {

    suppressWarnings(
      try(
        expr,
        silent = TRUE
      )
    )

  } else {

    try(
      expr,
      silent = TRUE
    )
  }
}


# -------------------------------------------------------------------------
# Internal helper: extract parameter estimates
# -------------------------------------------------------------------------

.compare_extract_estimates <- function(
    fit,
    model) {

  par <- fit$par

  if (is.null(par) ||
      length(par) == 0L) {

    return(
      data.frame(
        Model = character(0),
        Parameter = character(0),
        Estimate = numeric(0),
        SE = numeric(0),
        stringsAsFactors = FALSE
      )
    )
  }

  par_names <- names(par)

  if (is.null(par_names)) {

    par_names <- paste0(
      "par",
      seq_along(par)
    )
  }


  # -----------------------------------------------------------------------
  # Standard errors
  # -----------------------------------------------------------------------

  se <- rep(
    NA_real_,
    length(par)
  )

  names(se) <- par_names

  if (!is.null(fit$se)) {

    fit_se <- fit$se

    if (!is.null(names(fit_se))) {

      common <- intersect(
        par_names,
        names(fit_se)
      )

      se[common] <-
        as.numeric(
          fit_se[common]
        )

    } else if (
      length(fit_se) ==
      length(par)
    ) {

      se <-
        as.numeric(
          fit_se
        )

      names(se) <-
        par_names
    }
  }


  data.frame(
    Model = model,
    Parameter = par_names,
    Estimate = as.numeric(par),
    SE = as.numeric(se),
    stringsAsFactors = FALSE
  )
}


# -------------------------------------------------------------------------
# Internal helper: identify best model
# -------------------------------------------------------------------------

.compare_best_model <- function(
    comparison,
    criterion) {

  values <- comparison[[criterion]]

  valid <-
    is.finite(values)

  if (!any(valid)) {

    return(
      NA_character_
    )
  }

  ids <-
    which(valid)

  comparison$Model[
    ids[
      which.min(
        values[ids]
      )
    ]
  ]
}


# -------------------------------------------------------------------------
# Main comparison function
# -------------------------------------------------------------------------

#' Comparison of bounded-count models
#'
#' Fits and compares the binomial, beta-binomial,
#' Kumaraswamy-binomial, Conway-Maxwell-binomial, and
#' Lambert-uniform binomial distributions for bounded count data.
#'
#' @param x Integer vector of observations in \code{0, ..., size}.
#' @param size Positive integer giving the fixed upper bound of the support.
#' @param lub.control Optional named list of arguments passed to
#'   \code{\link{fitLUB_mle}}.
#' @param suppress_warnings Logical; if \code{TRUE}, warnings generated
#'   during model fitting are suppressed.
#'
#' @details
#' The B, BB, KB, CMB, and LUB models are fitted by maximum likelihood.
#' Model comparison is based on the maximized log-likelihood, Akaike
#' information criterion (AIC), corrected Akaike information criterion
#' (AICc), and Bayesian information criterion (BIC).
#'
#' Models are ordered according to their AIC values.
#'
#' @return An object of class \code{"compareLUBmodels"} containing the
#'   model-comparison table, parameter estimates and standard errors,
#'   fitted model objects, and the best model according to AIC, AICc,
#'   and BIC.
#'
#' @examples
#' x <- rep(
#'   0:4,
#'   times = c(69, 30, 9, 2, 0)
#' )
#'
#' fits <- compareLUBmodels(
#'   x = x,
#'   size = 4
#' )
#'
#' summary(fits)
#'
#' @seealso \code{\link{fitLUB_mle}}
#'
#' @export
compareLUBmodels <- function(
    x,
    size,
    lub.control = list(),
    suppress_warnings = TRUE) {


  # -----------------------------------------------------------------------
  # Input validation
  # -----------------------------------------------------------------------

  .check_LUB_sample(
    x = x,
    size = size
  )

  x <-
    as.integer(x)


  if (!is.list(
    lub.control
  )) {

    stop(
      "'lub.control' must be a list.",
      call. = FALSE
    )
  }


  if (length(
    suppress_warnings
  ) != 1L ||
      !is.logical(
        suppress_warnings
      ) ||
      is.na(
        suppress_warnings
      )) {

    stop(
      "'suppress_warnings' must be TRUE or FALSE.",
      call. = FALSE
    )
  }


  # -----------------------------------------------------------------------
  # Validate names supplied through lub.control
  # -----------------------------------------------------------------------

  if (length(
    lub.control
  ) > 0L) {

    if (is.null(
      names(
        lub.control
      )
    ) ||
        any(
          names(
            lub.control
          ) == ""
        )) {

      stop(
        "'lub.control' must be a named list.",
        call. = FALSE
      )
    }


    reserved <-
      intersect(
        names(
          lub.control
        ),
        c(
          "x",
          "size"
        )
      )

    if (length(
      reserved
    ) > 0L) {

      stop(
        paste0(
          "'lub.control' must not contain: ",
          paste(
            reserved,
            collapse = ", "
          ),
          "."
        ),
        call. = FALSE
      )
    }
  }


  n <-
    length(x)


  # -----------------------------------------------------------------------
  # Fit binomial model
  # -----------------------------------------------------------------------

  fit.bin <-
    .compare_fit_try(
      fitBIN_mle(
        x = x,
        size = size
      ),
      suppress_warnings =
        suppress_warnings
    )


  # -----------------------------------------------------------------------
  # Fit beta-binomial model
  # -----------------------------------------------------------------------

  fit.bb <-
    .compare_fit_try(
      fitBB_mle(
        x = x,
        size = size
      ),
      suppress_warnings =
        suppress_warnings
    )


  # -----------------------------------------------------------------------
  # Fit Kumaraswamy-binomial model
  # -----------------------------------------------------------------------

  fit.kb <-
    .compare_fit_try(
      fitKB_mle(
        x = x,
        size = size
      ),
      suppress_warnings =
        suppress_warnings
    )


  # -----------------------------------------------------------------------
  # Fit Conway-Maxwell-binomial model
  # -----------------------------------------------------------------------

  fit.cmb <-
    .compare_fit_try(
      fitCMB_mle(
        x = x,
        size = size
      ),
      suppress_warnings =
        suppress_warnings
    )


  # -----------------------------------------------------------------------
  # Fit LUB model
  # -----------------------------------------------------------------------

  fit.lub <-
    .compare_fit_try(
      do.call(
        fitLUB_mle,
        c(
          list(
            x = x,
            size = size
          ),
          lub.control
        )
      ),
      suppress_warnings =
        suppress_warnings
    )


  # -----------------------------------------------------------------------
  # Collect fitted models
  # -----------------------------------------------------------------------

  fits <- list(
    B = fit.bin,
    BB = fit.bb,
    KB = fit.kb,
    CMB = fit.cmb,
    LUB = fit.lub
  )


  # -----------------------------------------------------------------------
  # Check fitting failures
  # -----------------------------------------------------------------------

  failed <-
    vapply(
      fits,
      inherits,
      logical(1),
      what = "try-error"
    )

  if (any(
    failed
  )) {

    stop(
      paste0(
        "Model fitting failed for: ",
        paste(
          names(
            fits
          )[failed],
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }


  # -----------------------------------------------------------------------
  # Model names
  # -----------------------------------------------------------------------

  model_names <-
    names(
      fits
    )


  # -----------------------------------------------------------------------
  # Number of estimated parameters
  # -----------------------------------------------------------------------

  npar_values <-
    vapply(
      fits,
      function(fit) {

        length(
          fit$par
        )
      },
      integer(1)
    )


  # -----------------------------------------------------------------------
  # Maximized log-likelihood
  # -----------------------------------------------------------------------

  logLik_values <-
    vapply(
      fits,
      function(fit) {

        as.numeric(
          fit$logLik
        )[1L]
      },
      numeric(1)
    )


  # -----------------------------------------------------------------------
  # Information criteria
  # -----------------------------------------------------------------------

  AIC_values <-
    vapply(
      fits,
      function(fit) {

        as.numeric(
          fit$AIC
        )[1L]
      },
      numeric(1)
    )


  AICc_values <-
    vapply(
      fits,
      function(fit) {

        as.numeric(
          fit$AICc
        )[1L]
      },
      numeric(1)
    )


  BIC_values <-
    vapply(
      fits,
      function(fit) {

        as.numeric(
          fit$BIC
        )[1L]
      },
      numeric(1)
    )


  # -----------------------------------------------------------------------
  # Convergence information
  # -----------------------------------------------------------------------

  convergence_values <-
    vapply(
      fits,
      function(fit) {

        if (is.null(
          fit$convergence
        )) {

          return(
            NA_integer_
          )
        }

        as.integer(
          fit$convergence
        )[1L]
      },
      integer(1)
    )


  # -----------------------------------------------------------------------
  # Model-comparison table
  #
  # No DeltaAIC, DeltaAICc, or DeltaBIC values are reported.
  # -----------------------------------------------------------------------

  comparison <- data.frame(
    Model = model_names,
    npar = npar_values,
    logLik = logLik_values,
    AIC = AIC_values,
    AICc = AICc_values,
    BIC = BIC_values,
    convergence = convergence_values,
    stringsAsFactors = FALSE
  )


  # -----------------------------------------------------------------------
  # Order models according to AIC
  # -----------------------------------------------------------------------

  comparison <-
    comparison[
      order(
        comparison$AIC,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

  rownames(
    comparison
  ) <- NULL


  # -----------------------------------------------------------------------
  # Parameter estimates
  # -----------------------------------------------------------------------

  estimates <-
    do.call(
      rbind,
      Map(
        .compare_extract_estimates,
        fit = fits,
        model = model_names
      )
    )

  rownames(
    estimates
  ) <- NULL


  # -----------------------------------------------------------------------
  # Best models according to information criteria
  # -----------------------------------------------------------------------

  best.AIC <-
    .compare_best_model(
      comparison = comparison,
      criterion = "AIC"
    )

  best.AICc <-
    .compare_best_model(
      comparison = comparison,
      criterion = "AICc"
    )

  best.BIC <-
    .compare_best_model(
      comparison = comparison,
      criterion = "BIC"
    )


  # -----------------------------------------------------------------------
  # Output
  # -----------------------------------------------------------------------

  out <- list(
    comparison = comparison,
    estimates = estimates,
    fits = fits,
    n = n,
    size = size,
    best_AIC = best.AIC,
    best_AICc = best.AICc,
    best_BIC = best.BIC,
    suppress_warnings =
      suppress_warnings,
    lub.control =
      lub.control,
    call =
      match.call()
  )


  class(
    out
  ) <- c(
    "compareLUBmodels",
    "list"
  )


  out
}


# -------------------------------------------------------------------------
# Summary method
# -------------------------------------------------------------------------

#' Summary method for LUB model comparisons
#'
#' @param object An object of class \code{"compareLUBmodels"}.
#' @param digits Number of digits used for printing.
#' @param ... Additional arguments.
#'
#' @return The input object, returned invisibly.
#'
#' @method summary compareLUBmodels
#' @export
summary.compareLUBmodels <- function(
    object,
    digits = 4,
    ...) {


  if (!inherits(
    object,
    "compareLUBmodels"
  )) {

    stop(
      "'object' must be of class 'compareLUBmodels'.",
      call. = FALSE
    )
  }


  if (length(
    digits
  ) != 1L ||
      !is.numeric(
        digits
      ) ||
      !is.finite(
        digits
      ) ||
      digits < 0) {

    stop(
      "'digits' must be a non-negative finite number.",
      call. = FALSE
    )
  }

  digits <-
    as.integer(
      digits
    )


  cat("\n")
  cat(
    "Comparison of Bounded-Count Distributions\n"
  )
  cat(
    "-----------------------------------------\n\n"
  )


  # -----------------------------------------------------------------------
  # Call and data information
  # -----------------------------------------------------------------------

  cat(
    "Call:\n"
  )

  print(
    object$call
  )

  cat("\n")

  cat(
    "Sample size:",
    object$n,
    "\n"
  )

  cat(
    "Support    : 0, ... ,",
    object$size,
    "\n"
  )


  # -----------------------------------------------------------------------
  # Model comparison
  # -----------------------------------------------------------------------

  cat("\n")

  cat(
    "Model comparison:\n\n"
  )


  comparison <-
    object$comparison[
      ,
      c(
        "Model",
        "npar",
        "logLik",
        "AIC",
        "AICc",
        "BIC",
        "convergence"
      ),
      drop = FALSE
    ]


  numeric_cols <-
    vapply(
      comparison,
      is.numeric,
      logical(1)
    )

  comparison[
    numeric_cols
  ] <-
    lapply(
      comparison[
        numeric_cols
      ],
      function(z) {

        round(
          z,
          digits
        )
      }
    )


  print(
    comparison,
    row.names = FALSE
  )


  # -----------------------------------------------------------------------
  # Parameter estimates
  # -----------------------------------------------------------------------

  cat("\n")

  cat(
    "Parameter estimates:\n\n"
  )


  estimates <-
    object$estimates


  estimates$Estimate <-
    round(
      estimates$Estimate,
      digits
    )

  estimates$SE <-
    round(
      estimates$SE,
      digits
    )


  print(
    estimates,
    row.names = FALSE
  )


  # -----------------------------------------------------------------------
  # Best models
  # -----------------------------------------------------------------------

  cat("\n")

  cat(
    "Best model according to AIC :",
    object$best_AIC,
    "\n"
  )

  cat(
    "Best model according to AICc:",
    object$best_AICc,
    "\n"
  )

  cat(
    "Best model according to BIC :",
    object$best_BIC,
    "\n"
  )


  invisible(
    object
  )
}


# -------------------------------------------------------------------------
# Print method
# -------------------------------------------------------------------------

#' Print method for LUB model comparisons
#'
#' @param x An object of class \code{"compareLUBmodels"}.
#' @param digits Number of digits used for printing.
#' @param ... Additional arguments.
#'
#' @return The input object, returned invisibly.
#'
#' @method print compareLUBmodels
#' @export
print.compareLUBmodels <- function(
    x,
    digits = 4,
    ...) {


  summary.compareLUBmodels(
    object = x,
    digits = digits,
    ...
  )


  invisible(
    x
  )
}
