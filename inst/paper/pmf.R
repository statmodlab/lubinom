# -------------------------------------------------------------------------
# PMF figures for the LUB distribution
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(grid)
library(scales)
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
# Legend labels
# -------------------------------------------------------------------------

beta_labels <- c(
  "0.5" = "0.5",
  "1"   = "1",
  "3"   = "3",
  "8"   = "8",
  "13"  = "13"
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

plot_df <- bind_rows(
  lapply(
    size_values,
    function(size_value) {

      x_grid <- 0:size_value

      bind_rows(
        lapply(
          beta_values,
          function(beta_value) {

            data.frame(
              size = size_value,
              x = x_grid,
              beta = factor(
                as.character(beta_value),
                levels = beta_levels
              ),
              pmf = dLUB(
                x = x_grid,
                size = size_value,
                beta = beta_value
              )
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
  nrow(plot_df) > 0L,
  all(is.finite(plot_df$pmf)),
  all(plot_df$pmf >= 0),
  all(plot_df$pmf <= 1)
)

pmf_sums <- aggregate(
  pmf ~ size + beta,
  data = plot_df,
  FUN = sum
)

stopifnot(
  all(
    is.finite(
      pmf_sums$pmf
    )
  ),
  max(
    abs(
      pmf_sums$pmf - 1
    )
  ) < 1e-10
)


# -------------------------------------------------------------------------
# Check discrete-uniform reference case
# -------------------------------------------------------------------------

uniform_check <- plot_df[
  plot_df$beta == "1",
  ,
  drop = FALSE
]

uniform_diff <- vapply(
  split(
    uniform_check,
    uniform_check$size
  ),
  function(dat) {

    max(
      abs(
        dat$pmf -
          1 / (
            dat$size[1L] + 1
          )
      )
    )
  },
  numeric(1)
)

stopifnot(
  max(uniform_diff) <
    1e-12
)


# -------------------------------------------------------------------------
# Common theme
# -------------------------------------------------------------------------

theme_lub_pmf <- theme_bw(
  base_size = 14
) +
  theme(

    panel.grid.minor =
      element_blank(),

    panel.grid.major.x =
      element_line(
        colour = "grey90",
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

    plot.title =
      element_text(
        hjust = 0.5,
        size = 14
      ),

    legend.position =
      c(
        0.96,
        0.96
      ),

    legend.justification =
      c(
        1,
        1
      ),

    legend.background =
      element_rect(
        fill = alpha(
          "white",
          0.82
        ),
        colour = NA
      ),

    legend.key =
      element_rect(
        fill = "transparent",
        colour = NA
      ),

    legend.title =
      element_text(
        size = 11,
        hjust = 0.5
      ),

    legend.text =
      element_text(
        size = 10
      ),

    legend.key.width =
      unit(
        0.9,
        "cm"
      ),

    legend.key.height =
      unit(
        0.45,
        "cm"
      ),

    legend.spacing.y =
      unit(
        0.05,
        "cm"
      ),

    plot.margin =
      margin(
        t = 8,
        r = 8,
        b = 8,
        l = 8
      )
  )


# -------------------------------------------------------------------------
# Function to create each PMF panel
# -------------------------------------------------------------------------

plot_lub_pmf <- function(
    size_value) {

  df_plot <- plot_df[
    plot_df$size == size_value,
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

  ymax <- max(
    df_plot$pmf,
    na.rm = TRUE
  ) * 1.08

  ggplot(
    df_plot,
    aes(
      x = x,
      y = pmf,
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
      size = 2.7,
      stroke = 0.8
    ) +

    scale_colour_manual(
      name = expression(beta),
      values = beta_colors,
      breaks = beta_levels,
      labels = beta_labels
    ) +

    scale_linetype_manual(
      name = expression(beta),
      values = beta_linetypes,
      breaks = beta_levels,
      labels = beta_labels
    ) +

    scale_shape_manual(
      name = expression(beta),
      values = beta_shapes,
      breaks = beta_levels,
      labels = beta_labels
    ) +

    scale_x_continuous(
      breaks = 0:size_value,
      limits = c(
        -0.15,
        size_value + 0.15
      )
    ) +

    scale_y_continuous(
      expand = expansion(
        mult = c(
          0,
          0.02
        )
      )
    ) +

    coord_cartesian(
      ylim = c(
        -0.025,
        ymax
      )
    ) +

    labs(
      x = "Number of events",
      y = "Probability mass"
    ) +

    guides(
      colour = guide_legend(
        title = expression(beta),
        title.position = "top",
        title.hjust = 0.5,
        ncol = 1,
        byrow = TRUE
      ),

      linetype = guide_legend(
        title = expression(beta),
        title.position = "top",
        title.hjust = 0.5,
        ncol = 1,
        byrow = TRUE
      ),

      shape = guide_legend(
        title = expression(beta),
        title.position = "top",
        title.hjust = 0.5,
        ncol = 1,
        byrow = TRUE
      )
    ) +

    theme_lub_pmf
}


# -------------------------------------------------------------------------
# Figures
# -------------------------------------------------------------------------

p_n4 <- plot_lub_pmf(
  size_value = 4
)

p_n13 <- plot_lub_pmf(
  size_value = 13
)


# -------------------------------------------------------------------------
# Display
# -------------------------------------------------------------------------

print(p_n4)
print(p_n13)
