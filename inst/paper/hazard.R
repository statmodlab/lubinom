# -------------------------------------------------------------------------
# LUB: hazard function and mean residual count
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------------

library(ggplot2)
library(lubinom)


# -------------------------------------------------------------------------
# Parameter values
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
  13
)

beta_levels <- as.character(
  beta_values
)


# -------------------------------------------------------------------------
# Graphical specification
# -------------------------------------------------------------------------

beta_colors <- c(
  "0.5" = "#7A7A7A",
  "1"   = "#4C78A8",
  "3"   = "#72B7B2",
  "8"   = "#B279A2",
  "13"  = "#E45756"
)

# beta = 1 is highlighted as the discrete-uniform reference case.
beta_linetypes <- c(
  "0.5" = "solid",
  "1"   = "dashed",
  "3"   = "solid",
  "8"   = "solid",
  "13"  = "solid"
)

beta_shapes <- c(
  "0.5" = 1,
  "1"   = 0,
  "3"   = 2,
  "8"   = 5,
  "13"  = 16
)


# -------------------------------------------------------------------------
# Build plotting data
# -------------------------------------------------------------------------

dat <- do.call(
  rbind,
  lapply(
    size_values,
    function(size_value) {

      do.call(
        rbind,
        lapply(
          beta_values,
          function(beta_value) {

            x_grid <- 0:size_value

            data.frame(
              size = size_value,

              beta = factor(
                as.character(beta_value),
                levels = beta_levels
              ),

              x = x_grid,

              hazard = hLUB(
                x = x_grid,
                size = size_value,
                beta = beta_value
              ),

              mrc = mrlLUB(
                x = x_grid,
                size = size_value,
                beta = beta_value
              ),

              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

row.names(dat) <- NULL


# -------------------------------------------------------------------------
# Numerical checks
# -------------------------------------------------------------------------

stopifnot(
  nrow(dat) > 0L,

  all(
    is.finite(
      dat$hazard
    )
  ),

  all(
    is.finite(
      dat$mrc
    )
  ),

  all(
    dat$hazard >= 0
  ),

  all(
    dat$hazard <= 1
  ),

  all(
    dat$mrc >= -1e-12
  )
)


# -------------------------------------------------------------------------
# Endpoint identities
#
# h(size) = 1
# m(size) = 0
# -------------------------------------------------------------------------

endpoint_data <- dat[
  dat$x == dat$size,
  ,
  drop = FALSE
]

stopifnot(
  max(
    abs(
      endpoint_data$hazard - 1
    )
  ) < 1e-12,

  max(
    abs(
      endpoint_data$mrc
    )
  ) < 1e-12
)


# -------------------------------------------------------------------------
# Discrete-uniform reference case: beta = 1
#
# S(x) = (size - x + 1) / (size + 1)
#
# h(x) = 1 / (size - x + 1)
#
# m(x) = (size - x) / 2
# -------------------------------------------------------------------------

uniform_data <- dat[
  dat$beta == "1",
  ,
  drop = FALSE
]

hazard_reference <-
  1 /
  (
    uniform_data$size -
      uniform_data$x +
      1
  )

mrc_reference <-
  (
    uniform_data$size -
      uniform_data$x
  ) /
  2

stopifnot(
  max(
    abs(
      uniform_data$hazard -
        hazard_reference
    )
  ) < 1e-12,

  max(
    abs(
      uniform_data$mrc -
        mrc_reference
    )
  ) < 1e-12
)


# -------------------------------------------------------------------------
# Common theme
# -------------------------------------------------------------------------

theme_lub <- theme_bw(
  base_size = 14
) +
  theme(

    panel.grid.major =
      element_line(
        colour = "grey88",
        linewidth = 0.35
      ),

    panel.grid.minor =
      element_blank(),

    panel.border =
      element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.8
      ),

    axis.text =
      element_text(
        size = 12,
        colour = "black"
      ),

    axis.title =
      element_text(
        size = 14
      ),

    legend.background =
      element_rect(
        fill = scales::alpha(
          "white",
          0.85
        ),
        colour = NA
      ),

    legend.key =
      element_blank(),

    legend.title =
      element_text(
        hjust = 0.5
      ),

    plot.margin =
      margin(
        8,
        8,
        8,
        8
      )
  )


# -------------------------------------------------------------------------
# Hazard-function panel
# -------------------------------------------------------------------------

plot_hazard <- function(
    size_value,
    legend_position,
    legend_justification) {

  df_plot <- dat[
    dat$size == size_value,
    ,
    drop = FALSE
  ]

  if (nrow(df_plot) == 0L) {
    stop(
      "No observations found for size = ",
      size_value,
      ".",
      call. = FALSE
    )
  }

  ggplot(
    df_plot,
    aes(
      x = x,
      y = hazard,
      colour = beta,
      linetype = beta,
      shape = beta,
      group = beta
    )
  ) +

    geom_line(
      linewidth = 0.9
    ) +

    geom_point(
      size = 2.6
    ) +

    scale_colour_manual(
      values = beta_colors,
      name = expression(beta)
    ) +

    scale_linetype_manual(
      values = beta_linetypes,
      name = expression(beta)
    ) +

    scale_shape_manual(
      values = beta_shapes,
      name = expression(beta)
    ) +

    scale_x_continuous(
      breaks = 0:size_value,
      limits = c(
        -0.15,
        size_value + 0.15
      )
    ) +

    scale_y_continuous(
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      expand = expansion(
        mult = c(
          0,
          0.02
        )
      )
    ) +

    coord_cartesian(
      ylim = c(
        0,
        1.02
      )
    ) +

    labs(
      x = "x",
      y = "Hazard function"
    ) +

    theme_lub +

    theme(
      legend.position =
        legend_position,

      legend.justification =
        legend_justification
    )
}


# -------------------------------------------------------------------------
# Mean-residual-count panel
#
# m(x) = E(X - x | X >= x)
# -------------------------------------------------------------------------

plot_mrc <- function(
    size_value,
    ymin,
    ymax,
    legend_position,
    legend_justification) {

  df_plot <- dat[
    dat$size == size_value,
    ,
    drop = FALSE
  ]

  if (nrow(df_plot) == 0L) {
    stop(
      "No observations found for size = ",
      size_value,
      ".",
      call. = FALSE
    )
  }

  ggplot(
    df_plot,
    aes(
      x = x,
      y = mrc,
      colour = beta,
      linetype = beta,
      shape = beta,
      group = beta
    )
  ) +

    geom_line(
      linewidth = 0.9
    ) +

    geom_point(
      size = 2.6
    ) +

    scale_colour_manual(
      values = beta_colors,
      name = expression(beta)
    ) +

    scale_linetype_manual(
      values = beta_linetypes,
      name = expression(beta)
    ) +

    scale_shape_manual(
      values = beta_shapes,
      name = expression(beta)
    ) +

    scale_x_continuous(
      breaks = 0:size_value,
      limits = c(
        -0.15,
        size_value + 0.15
      )
    ) +

    scale_y_continuous(
      breaks = pretty(
        c(
          0,
          ymax
        )
      ),
      expand = expansion(
        mult = c(
          0,
          0.02
        )
      )
    ) +

    coord_cartesian(
      ylim = c(
        ymin,
        ymax
      )
    ) +

    labs(
      x = "x",
      y = "Mean residual count"
    ) +

    theme_lub +

    theme(
      legend.position =
        legend_position,

      legend.justification =
        legend_justification
    )
}


# -------------------------------------------------------------------------
# Hazard-function figures
# -------------------------------------------------------------------------

hazard_n4 <- plot_hazard(
  size_value = 4,
  legend_position = c(
    0.96,
    0.04
  ),
  legend_justification = c(
    1,
    0
  )
)

hazard_n13 <- plot_hazard(
  size_value = 13,
  legend_position = c(
    0.04,
    0.96
  ),
  legend_justification = c(
    0,
    1
  )
)


# -------------------------------------------------------------------------
# Mean-residual-count figures
# -------------------------------------------------------------------------

mrc_n4 <- plot_mrc(
  size_value = 4,
  ymin = -0.05,
  ymax = 2,
  legend_position = c(
    0.96,
    0.96
  ),
  legend_justification = c(
    1,
    1
  )
)

mrc_n13 <- plot_mrc(
  size_value = 13,
  ymin = -0.20,
  ymax = 7,
  legend_position = c(
    0.96,
    0.96
  ),
  legend_justification = c(
    1,
    1
  )
)


# -------------------------------------------------------------------------
# Display
# -------------------------------------------------------------------------

print(hazard_n4)
print(hazard_n13)

print(mrc_n4)
print(mrc_n13)
