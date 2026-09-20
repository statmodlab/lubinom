# -------------------------------------------------------------------------
# Skewness and excess kurtosis of the LUB distribution
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(grid)
library(lubinom)


# -------------------------------------------------------------------------
# Parameter grid
# -------------------------------------------------------------------------

beta_grid <- seq(
  0.1,
  15,
  length.out = 500
)

size_values <- c(
  4,
  8,
  13
)


# -------------------------------------------------------------------------
# First four raw moments computed from the LUB PMF
# -------------------------------------------------------------------------

raw_moments_LUB <- function(
    size,
    beta) {

  x <- 0:size

  probs <- dLUB(
    x = x,
    size = size,
    beta = beta
  )

  if (any(!is.finite(probs)) ||
      any(probs < 0)) {
    stop(
      "Invalid LUB probabilities encountered.",
      call. = FALSE
    )
  }

  total <- sum(probs)

  if (!is.finite(total) ||
      total <= 0) {
    stop(
      "Invalid LUB probability normalization.",
      call. = FALSE
    )
  }

  probs <- probs / total

  c(
    m1 = sum(
      x * probs
    ),

    m2 = sum(
      x^2 * probs
    ),

    m3 = sum(
      x^3 * probs
    ),

    m4 = sum(
      x^4 * probs
    )
  )
}


# -------------------------------------------------------------------------
# Shape measures
# -------------------------------------------------------------------------

shape_LUB <- function(
    size,
    beta) {

  moments <- raw_moments_LUB(
    size = size,
    beta = beta
  )

  m1 <- moments["m1"]
  m2 <- moments["m2"]
  m3 <- moments["m3"]
  m4 <- moments["m4"]

  mu2 <-
    m2 -
    m1^2

  mu3 <-
    m3 -
    3 * m1 * m2 +
    2 * m1^3

  mu4 <-
    m4 -
    4 * m1 * m3 +
    6 * m1^2 * m2 -
    3 * m1^4

  if (!is.finite(mu2) ||
      mu2 <= 0) {
    stop(
      "Non-positive variance encountered.",
      call. = FALSE
    )
  }

  skewness <-
    mu3 /
    mu2^(3 / 2)

  excess_kurtosis <-
    mu4 /
    mu2^2 -
    3

  data.frame(
    size = size,
    beta = beta,
    skewness = skewness,
    excess_kurtosis = excess_kurtosis,
    stringsAsFactors = FALSE
  )
}


# -------------------------------------------------------------------------
# Compute shape grid
# -------------------------------------------------------------------------

shape_df <- bind_rows(
  lapply(
    size_values,
    function(size_value) {

      bind_rows(
        lapply(
          beta_grid,
          function(beta_value) {

            shape_LUB(
              size = size_value,
              beta = beta_value
            )
          }
        )
      )
    }
  )
)


# -------------------------------------------------------------------------
# Numerical checks
# -------------------------------------------------------------------------

stopifnot(
  nrow(shape_df) ==
    length(beta_grid) *
    length(size_values),

  all(
    is.finite(
      shape_df$skewness
    )
  ),

  all(
    is.finite(
      shape_df$excess_kurtosis
    )
  )
)


# -------------------------------------------------------------------------
# Format plotting data
# -------------------------------------------------------------------------

shape_df <- shape_df %>%
  mutate(
    size = factor(
      size,
      levels = size_values
    )
  )


# -------------------------------------------------------------------------
# Discrete-uniform reference case: beta = 1
#
# For a discrete uniform distribution on {0, ..., size},
#
# skewness = 0
#
# excess kurtosis
#   = -6 * (size^2 + 2 * size + 2)
#       / (5 * size * (size + 2))
# -------------------------------------------------------------------------

uniform_reference <- data.frame(
  size = size_values,

  skewness = 0,

  excess_kurtosis =
    -6 *
    (
      size_values^2 +
      2 * size_values +
      2
    ) /
    (
      5 *
      size_values *
      (size_values + 2)
    )
)

uniform_computed <- bind_rows(
  lapply(
    size_values,
    function(size_value) {

      shape_LUB(
        size = size_value,
        beta = 1
      )
    }
  )
)

stopifnot(
  max(
    abs(
      uniform_computed$skewness -
      uniform_reference$skewness
    )
  ) < 1e-12,

  max(
    abs(
      uniform_computed$excess_kurtosis -
      uniform_reference$excess_kurtosis
    )
  ) < 1e-12
)


# -------------------------------------------------------------------------
# Display reference values
# -------------------------------------------------------------------------

cat(
  "\nValues at beta = 1\n\n"
)

print(
  uniform_computed,
  row.names = FALSE,
  digits = 6
)


# -------------------------------------------------------------------------
# Graphical specification
# -------------------------------------------------------------------------

size_colors <- c(
  "4"  = "#4C78A8",
  "8"  = "#72B7B2",
  "13" = "#E45756"
)

size_linetypes <- c(
  "4"  = "solid",
  "8"  = "dashed",
  "13" = "dotdash"
)


# -------------------------------------------------------------------------
# Common theme
# -------------------------------------------------------------------------

base_theme <- theme_bw(
  base_size = 14
) +
  theme(

    panel.grid.minor =
      element_blank(),

    panel.grid.major =
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
        size = 13
      ),

    axis.text =
      element_text(
        size = 11,
        colour = "black"
      ),

    legend.position =
      "bottom",

    legend.direction =
      "horizontal",

    legend.title =
      element_text(
        size = 11
      ),

    legend.text =
      element_text(
        size = 10
      ),

    legend.box =
      "horizontal",

    legend.key.width =
      unit(
        1.1,
        "cm"
      )
  )


# -------------------------------------------------------------------------
# Skewness
# -------------------------------------------------------------------------

p_skew <- ggplot(
  shape_df,
  aes(
    x = beta,
    y = skewness,
    colour = size,
    linetype = size,
    group = size
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5,
    colour = "grey45"
  ) +

  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey45"
  ) +

  geom_line(
    linewidth = 0.95
  ) +

  scale_colour_manual(
    name = "n",
    values = size_colors
  ) +

  scale_linetype_manual(
    name = "n",
    values = size_linetypes
  ) +

  scale_x_continuous(
    breaks = c(
      0,
      1,
      3,
      5,
      8,
      10,
      13,
      15
    )
  ) +

  labs(
    x = expression(beta),
    y = "Skewness"
  ) +

  base_theme


# -------------------------------------------------------------------------
# Excess kurtosis
# -------------------------------------------------------------------------

p_kurt <- ggplot(
  shape_df,
  aes(
    x = beta,
    y = excess_kurtosis,
    colour = size,
    linetype = size,
    group = size
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5,
    colour = "grey45"
  ) +

  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey45"
  ) +

  geom_line(
    linewidth = 0.95
  ) +

  scale_colour_manual(
    name = "n",
    values = size_colors
  ) +

  scale_linetype_manual(
    name = "n",
    values = size_linetypes
  ) +

  scale_x_continuous(
    breaks = c(
      0,
      1,
      3,
      5,
      8,
      10,
      13,
      15
    )
  ) +

  labs(
    x = expression(beta),
    y = "Excess kurtosis"
  ) +

  base_theme


# -------------------------------------------------------------------------
# Display
# -------------------------------------------------------------------------

print(p_skew)
print(p_kurt)

