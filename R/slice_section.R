################################################################################
#
# Cross-section of one posterior slice (R/slices.R), showing Laplace's
# "height x width". For a Gaussian slice, the rectangle of peak height and
# width sqrt(2 pi) sigma has exactly the slice's area, so the two tails
# (under the bell, outside the rectangle) fill the two corners (inside the
# rectangle, above the bell). Both are drawn in the same colour.
#
################################################################################

## ----- Configuration ---------------------------------------------------------
library(tidyverse)

col_inside <- "#F18F00" # bell inside the rectangle, orange as in R/slices.R
col_swap <- "#00A6AA" # tails and corners, equal areas, turquoise as tails
col_curve <- "#b10f2e" # slice and its peak, red as the spine in R/slices.R
col_arrow <- "black" # height and width
col_box <- "gray45"
col_move <- "#00777A" # tails-into-corners arrows (dark turquoise)

zmax <- 3.6
half_width <- sqrt(2 * pi) / 2 # rectangle area = peak x width = 1
peak <- dnorm(0)

## ----- Regions ---------------------------------------------------------------
bell <- tibble(x = seq(-zmax, zmax, length.out = 801), y = dnorm(x))
inside <- bell |> filter(abs(x) <= half_width)
tails <- bell |>
  filter(abs(x) >= half_width) |>
  mutate(side = if_else(x < 0, "left", "right"))

## ----- Plot ------------------------------------------------------------------
y_width <- -0.045 # baseline offset of the width arrow
arrow_both <- arrow(ends = "both", length = unit(0.18, "cm"), type = "closed")
arrow_move <- arrow(length = unit(0.16, "cm"), type = "closed")

ggplot() +
  geom_area(data = inside, aes(x, y), fill = col_inside, alpha = 0.5) +
  geom_ribbon(
    data = inside,
    aes(x, ymin = y, ymax = peak),
    fill = col_swap,
    alpha = 0.45
  ) +
  geom_area(
    data = tails,
    aes(x, y, group = side),
    fill = col_swap,
    alpha = 0.45
  ) +
  annotate(
    "rect",
    xmin = -half_width,
    xmax = half_width,
    ymin = 0,
    ymax = peak,
    fill = NA,
    colour = col_box,
    linewidth = 0.7
  ) +
  geom_line(data = bell, aes(x, y), colour = col_curve, linewidth = 1.2) +
  annotate(
    "segment",
    x = -zmax,
    xend = zmax,
    y = 0,
    yend = 0,
    colour = "gray40"
  ) +
  annotate(
    "segment",
    x = 0,
    xend = 0,
    y = 0,
    yend = peak - 0.012,
    colour = col_arrow,
    linewidth = 0.7,
    arrow = arrow_both
  ) +
  annotate("point", x = 0, y = peak, colour = col_curve, size = 2.6) +
  annotate(
    "text",
    x = -0.08,
    y = 0.3 * peak,
    label = "height",
    colour = col_arrow,
    size = 6.5,
    hjust = 1
  ) +
  annotate(
    "segment",
    x = -half_width,
    xend = half_width,
    y = y_width,
    yend = y_width,
    colour = col_arrow,
    linewidth = 0.7,
    arrow = arrow_both
  ) +
  annotate(
    "text",
    x = 0,
    y = y_width - 0.035,
    label = "width",
    colour = col_arrow,
    size = 6.5
  ) +
  annotate(
    "text",
    x = zmax,
    y = y_width,
    label = "vartheta[-j]",
    parse = TRUE,
    colour = "gray40",
    size = 6,
    hjust = 1
  ) +
  annotate(
    "curve",
    x = -1.85,
    xend = -0.98,
    y = 0.035,
    yend = 0.33,
    curvature = -0.35,
    colour = col_move,
    linewidth = 0.6,
    linetype = "dashed",
    arrow = arrow_move
  ) +
  annotate(
    "curve",
    x = 1.85,
    xend = 0.98,
    y = 0.035,
    yend = 0.33,
    curvature = 0.35,
    colour = col_move,
    linewidth = 0.6,
    linetype = "dashed",
    arrow = arrow_move
  ) +
  coord_cartesian(
    xlim = c(-zmax, zmax),
    ylim = c(y_width - 0.06, peak + 0.01),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void() +
  theme(plot.margin = margin(4, 12, 4, 4))
