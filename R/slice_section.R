################################################################################
#
# Cross-section of the central posterior slice in R/slices.R, showing Laplace's
# "height x width". For a Gaussian slice, the rectangle of peak height and
# width sqrt(2 pi) sigma has exactly the slice's area, so the two tails
# (under the bell, outside the rectangle) fill the two corners (inside the
# rectangle, above the bell). Tails and corners share the tail colour.
#
################################################################################

## ----- Configuration ---------------------------------------------------------
library(tidyverse)

col_inside <- "#f1a700" # bell inside the rectangle, peak colour of R/slices.R
col_swap <- "#3ca9a1" # tails and corners, tail colour of R/slices.R
col_peak <- "#b10f2e" # red as the spine in R/slices.R
col_ink <- "black" # slice, height and width
col_curve <- "gray30" # outline of the slice, as in R/slice_section_anim.R
col_box <- "gray45"
col_move <- "gray25" # tails-into-corners arrows
the_alpha <- 0.79

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
  geom_area(data = inside, aes(x, y), fill = col_inside, alpha = the_alpha) +
  geom_ribbon(data = inside, aes(x, ymin = y, ymax = peak), fill = col_swap, alpha = the_alpha) +
  geom_area(data = tails, aes(x, y, group = side), fill = col_swap, alpha = the_alpha) +
  geom_line(data = bell, aes(x, y), colour = col_curve, linewidth = 0.7) +
  annotate(
    "rect",
    xmin = -half_width,
    xmax = half_width,
    ymin = 0,
    ymax = peak,
    fill = NA,
    colour = "black",
    linewidth = 0.7
  ) +
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
    colour = col_ink,
    linewidth = 0.7,
    arrow = arrow_both
  ) +
  annotate("point", x = 0, y = peak, colour = col_peak, size = 2.6) +
  annotate(
    "text",
    x = -0.08,
    y = 0.3 * peak,
    label = "height",
    colour = col_ink,
    size = 6.5,
    hjust = 1
  ) +
  annotate(
    "segment",
    x = -half_width,
    xend = half_width,
    y = y_width,
    yend = y_width,
    colour = col_ink,
    linewidth = 0.7,
    arrow = arrow_both
  ) +
  annotate(
    "text",
    x = 0,
    y = y_width - 0.035,
    label = "width",
    colour = col_ink,
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
  # annotate(
  #   "curve",
  #   x = -1.85,
  #   xend = -0.98,
  #   y = 0.035,
  #   yend = 0.33,
  #   curvature = -0.35,
  #   colour = col_move,
  #   linewidth = 0.6,
  #   linetype = "dashed",
  #   arrow = arrow_move
  # ) +
  # annotate(
  #   "curve",
  #   x = 1.85,
  #   xend = 0.98,
  #   y = 0.035,
  #   yend = 0.33,
  #   curvature = 0.35,
  #   colour = col_move,
  #   linewidth = 0.6,
  #   linetype = "dashed",
  #   arrow = arrow_move
  # ) +
  coord_cartesian(
    xlim = c(-zmax, zmax),
    ylim = c(y_width - 0.06, peak + 0.01),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void() +
  theme(plot.margin = margin(4, 12, 4, 4))
