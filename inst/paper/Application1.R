# -------------------------------------------------------------------------
# Application 1: family mortality data
# -------------------------------------------------------------------------
#
# Response:
#   Number of child deaths in families with exactly four children.
#
# Models:
#   B   = Binomial
#   BB  = Beta-binomial
#   KB  = Kumaraswamy-binomial
#   CMB = Conway-Maxwell-binomial
#   LUB = Lambert-uniform binomial
#
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------------

library(ggplot2)
library(tidyr)
library(dplyr)
library(grid)
library(lubinom)


# -------------------------------------------------------------------------
# Data
# -------------------------------------------------------------------------

size <- 4

freq <- c(
  69,
  30,
  9,
  2,
  0
)

x <- rep(
  0:size,
  times = freq
)

N <- length(x)

grid <- 0:size


# -------------------------------------------------------------------------
# Data checks
# -------------------------------------------------------------------------

stopifnot(
  length(freq) == size + 1L,
  all(freq >= 0),
  all(freq == floor(freq)),
  sum(freq) == N,
  length(x) > 0L,
  all(x >= 0),
  all(x <= size),
  all(x == floor(x))
)


# -------------------------------------------------------------------------
# Observed frequencies
# -------------------------------------------------------------------------

observed_df <- data.frame(
  x = grid,
  Frequency = freq,
  Proportion = freq / N
)

print(
  observed_df,
  row.names = FALSE
)


# -------------------------------------------------------------------------
# Model fitting and numerical comparison
# -------------------------------------------------------------------------

comp <- compareLUBmodels(
  x = x,
  size = size,
  lub.control = list(
    diagnostic_bfgs = TRUE,
    profile.ci = TRUE
  )
)

summary(comp)


# -------------------------------------------------------------------------
# Extract fitted models
# -------------------------------------------------------------------------

fit.bin <- comp$fits$B
fit.bb  <- comp$fits$BB
fit.kb  <- comp$fits$KB
fit.cmb <- comp$fits$CMB
fit.lub <- comp$fits$LUB


# -------------------------------------------------------------------------
# Fit checks
# -------------------------------------------------------------------------

fits <- list(
  B = fit.bin,
  BB = fit.bb,
  KB = fit.kb,
  CMB = fit.cmb,
  LUB = fit.lub
)

stopifnot(
  all(
    vapply(
      fits,
      function(z) {
        is.list(z)
      },
      logical(1)
    )
  ),

  all(
    vapply(
      fits,
      function(z) {
        length(z$par) >= 1L
      },
      logical(1)
    )
  ),

  all(
    vapply(
      fits,
      function(z) {
        is.finite(z$logLik)
      },
      logical(1)
    )
  )
)


# -------------------------------------------------------------------------
# LUB parametric-bootstrap goodness-of-fit
# -------------------------------------------------------------------------

gof.lub <- bootstrapGOF_LUB(
  x = x,
  size = size,
  B = 2000,
  seed = 123,
  progress = TRUE
)

print(gof.lub)


# -------------------------------------------------------------------------
# Common graphical specification
# -------------------------------------------------------------------------

model_levels <- c(
  "B",
  "BB",
  "KB",
  "CMB",
  "LUB"
)

line_types <- c(
  B   = "dashed",
  BB  = "dotdash",
  KB  = "solid",
  CMB = "longdash",
  LUB = "solid"
)

point_shapes <- c(
  B   = 1,
  BB  = 2,
  KB  = 0,
  CMB = 5,
  LUB = 16
)

line_colors <- c(
  B   = "#7A7A7A",
  BB  = "#4C78A8",
  KB  = "#72B7B2",
  CMB = "#B279A2",
  LUB = "#E45756"
)


# -------------------------------------------------------------------------
# Common theme
# -------------------------------------------------------------------------

theme_app <- theme_bw(
  base_size = 14
) +
  theme(

    panel.grid.minor =
      element_blank(),

    panel.grid.major.x =
      element_line(
        colour = "grey88",
        linewidth = 0.35
      ),

    panel.grid.major.y =
      element_line(
        colour = "grey88",
        linewidth = 0.35
      ),

    panel.border =
      element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.8
      ),

    axis.title =
      element_text(
        size = 14
      ),

    axis.text =
      element_text(
        size = 12,
        colour = "black"
      ),

    legend.position =
      "bottom",

    legend.title =
      element_blank(),

    legend.text =
      element_text(
        size = 11
      ),

    legend.key.width =
      unit(
        1.2,
        "cm"
      )
  )


# -------------------------------------------------------------------------
# Figure 1: empirical and fitted probability mass functions
# -------------------------------------------------------------------------

plot_prob_df <- data.frame(
  x = grid,
  Observed = freq / N
)

plot_prob_long <- data.frame(

  x = rep(
    grid,
    times = length(model_levels)
  ),

  Model = factor(
    rep(
      model_levels,
      each = length(grid)
    ),
    levels = model_levels
  ),

  Expected = c(

    dBIN(
      grid,
      size = size,
      prob = fit.bin$par["prob"]
    ),

    dBB(
      grid,
      size = size,
      alpha = fit.bb$par["alpha"],
      beta = fit.bb$par["beta"]
    ),

    dKB(
      grid,
      size = size,
      a = fit.kb$par["a"],
      b = fit.kb$par["b"]
    ),

    dCMB(
      grid,
      size = size,
      prob = fit.cmb$par["prob"],
      nu = fit.cmb$par["nu"]
    ),

    dLUB(
      grid,
      size = size,
      beta = fit.lub$par["beta"]
    )
  )
)


# -------------------------------------------------------------------------
# PMF checks
# -------------------------------------------------------------------------

stopifnot(
  all(
    is.finite(
      plot_prob_long$Expected
    )
  ),

  all(
    plot_prob_long$Expected >= 0
  ),

  all(
    plot_prob_long$Expected <= 1
  )
)

pmf_sums <- aggregate(
  Expected ~ Model,
  data = plot_prob_long,
  FUN = sum
)

stopifnot(
  max(
    abs(
      pmf_sums$Expected - 1
    )
  ) < 1e-8
)


# -------------------------------------------------------------------------
# Plot fitted PMFs
# -------------------------------------------------------------------------

p_prob <- ggplot() +

  geom_segment(
    data = plot_prob_df,
    aes(
      x = x,
      xend = x,
      y = 0,
      yend = Observed
    ),
    linewidth = 1.1,
    colour = "grey35"
  ) +

  geom_line(
    data = plot_prob_long,
    aes(
      x = x,
      y = Expected,
      colour = Model,
      linetype = Model,
      group = Model
    ),
    linewidth = 0.9
  ) +

  geom_point(
    data = plot_prob_long,
    aes(
      x = x,
      y = Expected,
      colour = Model,
      shape = Model
    ),
    size = 2.7,
    stroke = 0.8
  ) +

  scale_colour_manual(
    values = line_colors
  ) +

  scale_linetype_manual(
    values = line_types
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  scale_x_continuous(
    breaks = grid,
    labels = grid,
    limits = c(
      min(grid) - 0.2,
      max(grid) + 0.2
    )
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.03
      )
    )
  ) +

  coord_cartesian(
    ylim = c(
      -0.025,
      0.68
    )
  ) +

  labs(
    x = "Number of child deaths",
    y = paste0(
      "Probability of a family experiencing\n",
      "x child deaths"
    ),
    colour = NULL,
    linetype = NULL,
    shape = NULL
  ) +

  theme_app

print(p_prob)


# -------------------------------------------------------------------------
# Figure 2: threshold-conditioned probability
#
# h(x) = P(X = x | X >= x)
# -------------------------------------------------------------------------

grid_h <- 0:size

hazard_df <- data.frame(

  x = rep(
    grid_h,
    times = length(model_levels)
  ),

  Model = factor(
    rep(
      model_levels,
      each = length(grid_h)
    ),
    levels = model_levels
  ),

  Hazard = c(

    hBIN(
      grid_h,
      size = size,
      prob = fit.bin$par["prob"]
    ),

    hBB(
      grid_h,
      size = size,
      alpha = fit.bb$par["alpha"],
      beta = fit.bb$par["beta"]
    ),

    hKB(
      grid_h,
      size = size,
      a = fit.kb$par["a"],
      b = fit.kb$par["b"]
    ),

    hCMB(
      grid_h,
      size = size,
      prob = fit.cmb$par["prob"],
      nu = fit.cmb$par["nu"]
    ),

    hLUB(
      grid_h,
      size = size,
      beta = fit.lub$par["beta"]
    )
  )
)


# -------------------------------------------------------------------------
# Hazard checks
# -------------------------------------------------------------------------

stopifnot(
  all(
    is.finite(
      hazard_df$Hazard
    )
  ),

  all(
    hazard_df$Hazard >= 0
  ),

  all(
    hazard_df$Hazard <= 1
  )
)

hazard_endpoint <- hazard_df[
  hazard_df$x == size,
  ,
  drop = FALSE
]

stopifnot(
  max(
    abs(
      hazard_endpoint$Hazard - 1
    )
  ) < 1e-10
)


# -------------------------------------------------------------------------
# Plot threshold-conditioned probability
# -------------------------------------------------------------------------

p_hazard <- ggplot(
  hazard_df,
  aes(
    x = x,
    y = Hazard,
    colour = Model,
    linetype = Model,
    shape = Model,
    group = Model
  )
) +

  geom_line(
    linewidth = 0.9
  ) +

  geom_point(
    size = 2.7,
    stroke = 0.8
  ) +

  scale_colour_manual(
    values = line_colors
  ) +

  scale_linetype_manual(
    values = line_types
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  scale_x_continuous(
    breaks = grid_h,
    labels = grid_h,
    limits = c(
      min(grid_h) - 0.15,
      max(grid_h) + 0.15
    )
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.03
      )
    )
  ) +

  coord_cartesian(
    ylim = c(
      0.55,
      1.00
    )
  ) +

  labs(
    x = "Number of child deaths",
    y = paste0(
      "Probability of exactly x child deaths\n",
      "given at least x child deaths"
    ),
    colour = NULL,
    linetype = NULL,
    shape = NULL
  ) +

  theme_app

print(p_hazard)


# -------------------------------------------------------------------------
# Figure 3: mean residual count
#
# m(x) = E(X - x | X >= x)
# -------------------------------------------------------------------------

grid_mrc <- 0:size

mrc_df <- data.frame(

  x = rep(
    grid_mrc,
    times = length(model_levels)
  ),

  Model = factor(
    rep(
      model_levels,
      each = length(grid_mrc)
    ),
    levels = model_levels
  ),

  MRC = c(

    mrlBIN(
      grid_mrc,
      size = size,
      prob = fit.bin$par["prob"]
    ),

    mrlBB(
      grid_mrc,
      size = size,
      alpha = fit.bb$par["alpha"],
      beta = fit.bb$par["beta"]
    ),

    mrlKB(
      grid_mrc,
      size = size,
      a = fit.kb$par["a"],
      b = fit.kb$par["b"]
    ),

    mrlCMB(
      grid_mrc,
      size = size,
      prob = fit.cmb$par["prob"],
      nu = fit.cmb$par["nu"]
    ),

    mrlLUB(
      grid_mrc,
      size = size,
      beta = fit.lub$par["beta"]
    )
  )
)


# -------------------------------------------------------------------------
# Mean-residual-count checks
# -------------------------------------------------------------------------

stopifnot(
  all(
    is.finite(
      mrc_df$MRC
    )
  ),

  all(
    mrc_df$MRC >= -1e-10
  )
)

mrc_endpoint <- mrc_df[
  mrc_df$x == size,
  ,
  drop = FALSE
]

stopifnot(
  max(
    abs(
      mrc_endpoint$MRC
    )
  ) < 1e-10
)


# -------------------------------------------------------------------------
# Plot mean residual count
# -------------------------------------------------------------------------

p_mrc <- ggplot(
  mrc_df,
  aes(
    x = x,
    y = MRC,
    colour = Model,
    linetype = Model,
    shape = Model,
    group = Model
  )
) +

  geom_line(
    linewidth = 0.9
  ) +

  geom_point(
    size = 2.7,
    stroke = 0.8
  ) +

  scale_colour_manual(
    values = line_colors
  ) +

  scale_linetype_manual(
    values = line_types
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  scale_x_continuous(
    breaks = grid_mrc,
    labels = grid_mrc,
    limits = c(
      min(grid_mrc) - 0.15,
      max(grid_mrc) + 0.15
    )
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.03
      )
    )
  ) +

  coord_cartesian(
    ylim = c(
      -0.03,
      0.52
    )
  ) +

  labs(
    x = "Number of child deaths",
    y = paste0(
      "Expected additional child deaths\n",
      "given at least x child deaths"
    ),
    colour = NULL,
    linetype = NULL,
    shape = NULL
  ) +

  theme_app

print(p_mrc)

