# -------------------------------------------------------------------------
# Monte Carlo simulation study for the LUB model
#
# Finite-sample performance of the maximum likelihood estimator
# -------------------------------------------------------------------------

library(lubinom)


# -------------------------------------------------------------------------
# Expected event proportion
#
# E(X) / size = E(P)
# -------------------------------------------------------------------------

mean_prop_LUB <- function(beta) {

  stopifnot(
    length(beta) == 1L,
    is.numeric(beta),
    is.finite(beta),
    beta > 0
  )

  a <- 1 - beta

  # Stable evaluation around beta = 1.
  if (abs(a) < 1e-7) {

    return(
      0.5 +
        a / 6 +
        a^2 / 24 +
        a^3 / 120
    )
  }

  (expm1(a) - a) / a^2
}


# -------------------------------------------------------------------------
# Safe summary helpers
# -------------------------------------------------------------------------

safe_mean <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(x)
}


safe_sd <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) < 2L) {
    return(NA_real_)
  }

  stats::sd(x)
}


safe_rate <- function(x) {

  x <- x[
    !is.na(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(x)
}


# -------------------------------------------------------------------------
# Simulation settings
# -------------------------------------------------------------------------

beta_values <- c(
  0.5,
  1,
  3,
  8,
  13
)

size_values <- c(
  4,
  8,
  13
)

N_values <- c(
  25,
  50,
  100,
  200
)

R <- 2000L

conf.level <- 0.95

seed <- 123L

set.seed(seed)


# -------------------------------------------------------------------------
# Parameter regimes
# -------------------------------------------------------------------------

paper_table_regimes <- data.frame(

  beta = beta_values,

  MeanProp = vapply(
    beta_values,
    mean_prop_LUB,
    numeric(1)
  ),

  Regime = c(
    "Relatively high event proportion",
    "Discrete-uniform reference case",
    "Moderate event proportion",
    "Low event proportion",
    "Strong concentration at low counts"
  ),

  stringsAsFactors = FALSE
)


# -------------------------------------------------------------------------
# Simulation design
# -------------------------------------------------------------------------

design <- expand.grid(
  beta = beta_values,
  size = size_values,
  N = N_values,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

design <- design[
  order(
    design$size,
    design$beta,
    design$N
  ),
  ,
  drop = FALSE
]

row.names(design) <- NULL


# -------------------------------------------------------------------------
# Design checks
# -------------------------------------------------------------------------

stopifnot(
  nrow(design) ==
    length(beta_values) *
    length(size_values) *
    length(N_values),

  all(design$beta > 0),

  all(design$size > 0),

  all(design$N > 0)
)


# -------------------------------------------------------------------------
# Storage
# -------------------------------------------------------------------------

scenario_results <- vector(
  "list",
  nrow(design)
)

raw_results <- vector(
  "list",
  nrow(design)
)


# -------------------------------------------------------------------------
# Monte Carlo experiment
# -------------------------------------------------------------------------

cat(
  "\nRunning Monte Carlo simulation...\n\n"
)

for (s in seq_len(
  nrow(design)
)) {

  beta.true <- design$beta[s]
  size <- design$size[s]
  N <- design$N[s]


  # -----------------------------------------------------------------------
  # Progress
  # -----------------------------------------------------------------------

  cat(
    sprintf(
      paste0(
        "\rRunning scenario %2d/%2d | ",
        "beta = %4.1f | size = %2d | N = %3d"
      ),
      s,
      nrow(design),
      beta.true,
      size,
      N
    )
  )

  flush.console()


  # -----------------------------------------------------------------------
  # Replication-level storage
  # -----------------------------------------------------------------------

  sim <- data.frame(

    beta_hat = rep(
      NA_real_,
      R
    ),

    se_beta = rep(
      NA_real_,
      R
    ),

    wald_lower = rep(
      NA_real_,
      R
    ),

    wald_upper = rep(
      NA_real_,
      R
    ),

    cover_wald = rep(
      NA,
      R
    ),

    all_zero = rep(
      FALSE,
      R
    ),

    fit_attempted = rep(
      FALSE,
      R
    ),

    fit_ok = rep(
      FALSE,
      R
    ),

    hessian_ok = rep(
      FALSE,
      R
    ),

    se_ok = rep(
      FALSE,
      R
    ),

    stringsAsFactors = FALSE
  )


  # -----------------------------------------------------------------------
  # Replications
  # -----------------------------------------------------------------------

  for (r in seq_len(R)) {

    x <- rLUB(
      n = N,
      size = size,
      beta = beta.true
    )


    # ---------------------------------------------------------------------
    # Upper-boundary case
    #
    # For an all-zero sample, the likelihood is maximized only in the
    # limit beta -> infinity. Hence, no finite interior MLE exists.
    # ---------------------------------------------------------------------

    sim$all_zero[r] <-
      all(x == 0L)

    if (sim$all_zero[r]) {
      next
    }


    # ---------------------------------------------------------------------
    # Maximum likelihood estimation
    #
    # Brent is used as the official one-dimensional optimizer.
    # BFGS and profile-likelihood diagnostics are disabled because the
    # simulation study concerns the finite-sample behavior of the
    # primary estimator and its transformed Wald interval.
    # ---------------------------------------------------------------------

    sim$fit_attempted[r] <- TRUE

    fit <- try(
      fitLUB_mle(
        x = x,
        size = size,
        diagnostic_bfgs = FALSE,
        profile.ci = FALSE,
        conf.level = conf.level
      ),
      silent = TRUE
    )

    if (inherits(
      fit,
      "try-error"
    )) {
      next
    }


    # ---------------------------------------------------------------------
    # Numerical diagnostics
    # ---------------------------------------------------------------------

    sim$fit_ok[r] <-
      isTRUE(
        fit$success
      )

    sim$hessian_ok[r] <-
      isTRUE(
        fit$hessian_ok
      )

    sim$se_ok[r] <-
      isTRUE(
        fit$se_ok
      )

    if (!sim$fit_ok[r]) {
      next
    }


    # ---------------------------------------------------------------------
    # Point estimate
    # ---------------------------------------------------------------------

    beta.hat <- unname(
      fit$par["beta"]
    )

    if (!is.finite(beta.hat) ||
        beta.hat <= 0) {
      next
    }

    sim$beta_hat[r] <-
      beta.hat


    # ---------------------------------------------------------------------
    # Standard error
    # ---------------------------------------------------------------------

    se.beta <- unname(
      fit$se["beta"]
    )

    if (is.finite(se.beta) &&
        se.beta > 0) {

      sim$se_beta[r] <-
        se.beta
    }


    # ---------------------------------------------------------------------
    # Transformed Wald confidence interval
    # ---------------------------------------------------------------------

    if (
      isTRUE(
        fit$se_ok
      ) &&
      !is.null(
        fit$coefficients
      ) &&
      "beta" %in%
        rownames(
          fit$coefficients
        ) &&
      all(
        c(
          "Lower",
          "Upper"
        ) %in%
          colnames(
            fit$coefficients
          )
      )
    ) {

      ci <- fit$coefficients[
        "beta",
        c(
          "Lower",
          "Upper"
        )
      ]

      if (all(
        is.finite(ci)
      )) {

        sim$wald_lower[r] <-
          ci["Lower"]

        sim$wald_upper[r] <-
          ci["Upper"]

        sim$cover_wald[r] <-
          ci["Lower"] <= beta.true &&
          beta.true <= ci["Upper"]
      }
    }
  }


  # -----------------------------------------------------------------------
  # Valid finite interior estimates
  # -----------------------------------------------------------------------

  valid.beta <-
    !sim$all_zero &
    sim$fit_ok &
    is.finite(
      sim$beta_hat
    )

  valid.se <-
    valid.beta &
    sim$se_ok &
    is.finite(
      sim$se_beta
    ) &
    sim$se_beta > 0

  valid.ci <-
    valid.se &
    !is.na(
      sim$cover_wald
    )


  # -----------------------------------------------------------------------
  # Finite-sample summaries
  # -----------------------------------------------------------------------

  beta.estimates <-
    sim$beta_hat[
      valid.beta
    ]

  beta.errors <-
    beta.estimates -
    beta.true

  bias <-
    safe_mean(
      beta.errors
    )

  rmse <- if (
    length(
      beta.errors[
        is.finite(
          beta.errors
        )
      ]
    ) > 0L
  ) {

    sqrt(
      safe_mean(
        beta.errors^2
      )
    )

  } else {

    NA_real_
  }


  # -----------------------------------------------------------------------
  # Diagnostic denominators
  # -----------------------------------------------------------------------

  attempted <-
    sim$fit_attempted

  successful <-
    attempted &
    sim$fit_ok


  # -----------------------------------------------------------------------
  # Scenario summary
  # -----------------------------------------------------------------------

  scenario_results[[s]] <- data.frame(

    beta = beta.true,

    mean_proportion =
      mean_prop_LUB(
        beta.true
      ),

    size = size,

    N = N,

    R = R,


    # ---------------------------------------------------------------------
    # Finite-sample performance
    # ---------------------------------------------------------------------

    mean_beta_hat =
      safe_mean(
        beta.estimates
      ),

    bias = bias,

    rmse = rmse,

    ESE =
      safe_sd(
        beta.estimates
      ),

    ASE =
      safe_mean(
        sim$se_beta[
          valid.se
        ]
      ),

    CP =
      safe_rate(
        sim$cover_wald[
          valid.ci
        ]
      ),


    # ---------------------------------------------------------------------
    # Internal diagnostics
    #
    # Retained for reproducibility but not included in the manuscript
    # performance table.
    # ---------------------------------------------------------------------

    n_boundary =
      sum(
        sim$all_zero
      ),

    upper_boundary_rate =
      mean(
        sim$all_zero
      ),

    n_fit_attempted =
      sum(
        attempted
      ),

    n_fit_success =
      sum(
        successful
      ),

    fit_rate = if (
      any(attempted)
    ) {

      mean(
        sim$fit_ok[
          attempted
        ]
      )

    } else {

      NA_real_
    },

    hessian_rate = if (
      any(successful)
    ) {

      mean(
        sim$hessian_ok[
          successful
        ]
      )

    } else {

      NA_real_
    },

    se_rate = if (
      any(successful)
    ) {

      mean(
        sim$se_ok[
          successful
        ]
      )

    } else {

      NA_real_
    },

    ci_rate = if (
      any(valid.se)
    ) {

      mean(
        !is.na(
          sim$cover_wald[
            valid.se
          ]
        )
      )

    } else {

      NA_real_
    },

    interior_rate =
      mean(
        valid.beta
      ),

    stringsAsFactors = FALSE
  )


  # -----------------------------------------------------------------------
  # Store replication-level results
  # -----------------------------------------------------------------------

  raw_results[[s]] <-
    sim
}


cat(
  "\n\nSimulation completed.\n"
)


# -------------------------------------------------------------------------
# Combined results
# -------------------------------------------------------------------------

simulation_results <- do.call(
  rbind,
  scenario_results
)

row.names(
  simulation_results
) <- NULL


# -------------------------------------------------------------------------
# Global numerical checks
# -------------------------------------------------------------------------

stopifnot(
  nrow(simulation_results) ==
    nrow(design),

  all(
    simulation_results$
      upper_boundary_rate >= 0 &
    simulation_results$
      upper_boundary_rate <= 1
  ),

  all(
    simulation_results$
      interior_rate >= 0 &
    simulation_results$
      interior_rate <= 1
  )
)

finite_cp <- is.finite(
  simulation_results$CP
)

stopifnot(
  all(
    simulation_results$CP[
      finite_cp
    ] >= 0 &
    simulation_results$CP[
      finite_cp
    ] <= 1
  )
)


# -------------------------------------------------------------------------
# Manuscript table 1: parameter regimes
# -------------------------------------------------------------------------

paper_table_regimes$MeanProp <-
  round(
    paper_table_regimes$MeanProp,
    3
  )


# -------------------------------------------------------------------------
# Manuscript table 2: finite-sample performance
#
# Bias = empirical bias of beta_hat
# RMSE = root mean squared error
# ESE  = empirical standard error
# ASE  = average estimated standard error
# CP   = empirical coverage probability of the transformed
#        95% Wald confidence interval
#
# Boundary samples with all observations equal to zero are excluded
# because no finite interior MLE exists in those samples.
# -------------------------------------------------------------------------

paper_table_mle <- simulation_results[
  ,
  c(
    "beta",
    "mean_proportion",
    "size",
    "N",
    "bias",
    "rmse",
    "ESE",
    "ASE",
    "CP"
  ),
  drop = FALSE
]

names(
  paper_table_mle
) <- c(
  "beta",
  "MeanProp",
  "n",
  "N",
  "Bias",
  "RMSE",
  "ESE",
  "ASE",
  "CP"
)


# -------------------------------------------------------------------------
# Manuscript formatting
# -------------------------------------------------------------------------

paper_table_mle$MeanProp <-
  round(
    paper_table_mle$MeanProp,
    3
  )

paper_table_mle$Bias <-
  round(
    paper_table_mle$Bias,
    3
  )

paper_table_mle$RMSE <-
  round(
    paper_table_mle$RMSE,
    3
  )

paper_table_mle$ESE <-
  round(
    paper_table_mle$ESE,
    3
  )

paper_table_mle$ASE <-
  round(
    paper_table_mle$ASE,
    3
  )

paper_table_mle$CP <-
  round(
    paper_table_mle$CP,
    3
  )


# -------------------------------------------------------------------------
# Display manuscript results
# -------------------------------------------------------------------------

cat(
  "\n\nParameter regimes used in the simulation study\n\n"
)

print(
  paper_table_regimes,
  row.names = FALSE
)


cat(
  paste0(
    "\n\nFinite-sample performance of the ",
    "LUB maximum likelihood estimator\n"
  )
)


# -------------------------------------------------------------------------
# Split manuscript table by group size
# -------------------------------------------------------------------------

paper_tables_by_size <- split(
  paper_table_mle,
  paper_table_mle$n
)

for (nn in names(
  paper_tables_by_size
)) {

  tab <-
    paper_tables_by_size[[nn]]

  tab$n <- NULL

  cat(
    "\n\nn =",
    nn,
    "\n\n"
  )

  print(
    tab,
    row.names = FALSE
  )
}


# -------------------------------------------------------------------------
# Display boundary diagnostics
# -------------------------------------------------------------------------

boundary_table <- simulation_results[
  simulation_results$
    upper_boundary_rate > 0,
  c(
    "beta",
    "size",
    "N",
    "n_boundary",
    "upper_boundary_rate"
  ),
  drop = FALSE
]

if (nrow(
  boundary_table
) > 0L) {

  cat(
    "\n\nBoundary samples observed during the simulation\n\n"
  )

  print(
    boundary_table,
    row.names = FALSE
  )

} else {

  cat(
    "\n\nNo boundary samples were observed.\n"
  )
}
