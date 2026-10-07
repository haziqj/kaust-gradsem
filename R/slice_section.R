################################################################################
#
# Cross-section of the central posterior slice in R/slices.R, showing Laplace's
# "height x width". For a Gaussian slice, the rectangle of peak height and
# width sqrt(2 pi) sigma has exactly the slice's area, so the two tails
# (under the bell, outside the rectangle) fill the two corners (inside the
# rectangle, above the bell). Each corner is drawn as its tail flipped about
# the rectangle's side, in the tail's colours.
#
################################################################################

## ----- Configuration ---------------------------------------------------------
library(tidyverse)

# Same palette and interpolation as scale_fill_gradientn() in R/slices.R. There
# the fill is the density relative to the joint mode, which is also the peak of
# this (central) slice.
kaust_cols <- c("#00A6AA", "#CDCE00", "#F0B500", "#F18F00")
ombre <- scales::gradient_n_pal(kaust_cols)

col_peak <- "#b10f2e" # red as the spine in R/slices.R
col_ink <- "black" # slice, height and width
col_box <- "gray45"
col_move <- "gray25" # tails-into-corners arrows
corner_white <- 0.4 # share of white mixed into the corner colours

zmax <- 3.6
half_width <- sqrt(2 * pi) / 2 # rectangle area = peak x width = 1
peak <- dnorm(0)

## ----- Regions ---------------------------------------------------------------
# Thin vertical strips carry the gradient, as in ggridges. Each strip overlaps
# the next by half a step, which hides antialiasing seams.
step <- 0.01
mix_white <- function(col, w) {
  rgb(t(col2rgb(col) * (1 - w) + 255 * w), maxColorValue = 255)
}

bell <- tibble(x = seq(-zmax, zmax, by = step), y = dnorm(x))
under <- bell |> mutate(ymin = 0, ymax = y, fill = ombre(y / peak))
corners <- bell |>
  filter(abs(x) <= half_width) |>
  mutate(
    ymin = y,
    ymax = peak,
    fill = mix_white(ombre(dnorm(2 * half_width - abs(x)) / peak), corner_white)
  )
strips <- bind_rows(under, corners)

## ----- Plot ------------------------------------------------------------------
y_width <- -0.045 # baseline offset of the width arrow
arrow_both <- arrow(ends = "both", length = unit(0.18, "cm"), type = "closed")
arrow_move <- arrow(length = unit(0.16, "cm"), type = "closed")

ggplot() +
  geom_rect(
    data = strips,
    aes(xmin = x - step / 2, xmax = x + step, ymin = ymin, ymax = ymax),
    fill = strips$fill
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
  geom_line(data = bell, aes(x, y), colour = col_ink, linewidth = 1) +
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
