################################################################################
#
# Animated version of R/slice_section.R. An orange slice is pulled out of a
# point, as if lifted from the posterior in R/slices.R. Its peak height and
# width are then marked out, sweeping out the height x width rectangle. The two
# tails turn turquoise and pour into the empty corners of the rectangle, which
# ends up as one orange block of the same area. Run from the project root.
# Writes figures/slice_section.gif.
#
################################################################################

## ----- Configuration ---------------------------------------------------------
library(tidyverse)

col_fill <- "#f1a700" # the slice, peak colour of R/slices.R
col_swap <- "#3ca9a1" # tails and corners, tail colour of R/slices.R
col_peak <- "#b10f2e" # red as the spine in R/slices.R
col_ink <- "black" # height, width and the rectangle
col_curve <- "gray30" # outline of the slice, still visible once the tails drain
col_axis <- "gray40"
the_alpha <- 0.79 # fills are pre-blended onto white, so overlaps stay flat

zmax <- 3.6
half_width <- sqrt(2 * pi) / 2 # rectangle area = peak x width = 1
peak <- dnorm(0)
y_width <- -0.045 # baseline offset of the width arrow
y_lo <- y_width - 0.06 # bottom of the panel
grow_from <- c(x = zmax, y = y_lo) # bottom right, nearest the 3D plot
curve_fade <- 0.3 # final opacity of the outline, so the rectangle reads as one

fps <- 25
gif_file <- "figures/slice_section.gif"

## ----- Timeline --------------------------------------------------------------
# fmt: skip
timeline <- tribble(
  ~step,      ~secs,
  "blank",    0.2,
  "grow",     1.4,  # slice grows out of grow_from
  "label",    0.3,  # vartheta_{-j} fades in
  "pause",    0.3,
  "height",   0.7,  # height arrow rises from the baseline
  "dot",      0.25, # peak dot pops in
  "pause",    0.4,
  "width",    0.9,  # width arrow sweeps out the rectangle
  "pause",    0.6,
  "recolour", 0.6,  # tails turn turquoise
  "pause",    0.3,
  "pour",     1.8,  # tails drain into the corners
  "pause",    0.5,
  "merge",    0.9,  # corners turn orange, outline fades
  "hold",     3     # rest on the rectangle before looping
) |>
  mutate(start = cumsum(secs) - secs)

# Cubic ease-in-out, as in R/variational_anim.R
ease <- function(t) ifelse(t < 0.5, 4 * t^3, 1 - (-2 * t + 2)^3 / 2)

progress <- function(t, name) {
  row <- filter(timeline, step == name)
  ease(pmin(1, pmax(0, (t - row$start) / row$secs)))
}

frames <- tibble(t = seq(0, sum(timeline$secs), by = 1 / fps))
for (name in setdiff(timeline$step, c("blank", "pause", "hold"))) {
  frames[[name]] <- progress(frames$t, name)
}

## ----- Pour ------------------------------------------------------------------
# Each tail drains from its far end towards the rectangle, while the matching
# corner fills from the rectangle's edge towards the centre at the same rate
# of area. The turquoise is seen to pass through the edge into the corner.
# Right-hand side only, mirrored for the left.
tail_area <- pnorm(zmax) - pnorm(half_width)
corner_grid <- tibble(
  x = seq(0, half_width, length.out = 2001),
  area = peak * (half_width - x) - (pnorm(half_width) - pnorm(x))
)
corner_area <- max(corner_grid$area)

frames <- frames |>
  mutate(
    x_drain = qnorm(pnorm(half_width) + (1 - pour) * tail_area),
    x_fill = approx(
      corner_grid$area,
      corner_grid$x,
      xout = pour * corner_area
    )$y
  )

## ----- Helpers ---------------------------------------------------------------
mix_col <- function(from, to, p) {
  rgb(t((1 - p) * col2rgb(from) + p * col2rgb(to)), maxColorValue = 255)
}
orange <- mix_col("white", col_fill, the_alpha)
turquoise <- mix_col("white", col_swap, the_alpha)

# Blend in polar Luv (as hcl()), which passes through yellow-green like the
# gradient in R/slices.R. A straight RGB blend goes through a muddy olive.
to_hcl <- function(col) {
  luv <- convertColor(t(col2rgb(col)) / 255, from = "sRGB", to = "Luv")
  c(
    l = luv[1],
    c = sqrt(luv[2]^2 + luv[3]^2),
    h = atan2(luv[3], luv[2]) * 180 / pi
  )
}
blend_hcl <- function(from, to, p) {
  if (p == 0) {
    return(from)
  }
  if (p == 1) {
    return(to)
  }
  a <- to_hcl(from)
  b <- to_hcl(to)
  dh <- (b[["h"]] - a[["h"]] + 180) %% 360 - 180 # shorter way round
  hcl(
    h = a[["h"]] + p * dh,
    c = (1 - p) * a[["c"]] + p * b[["c"]],
    l = (1 - p) * a[["l"]] + p * b[["l"]],
    fixup = TRUE
  )
}

# Region under the slice, or between the slice and its peak, from a to b
under <- function(a, b, id) {
  xs <- seq(a, b, length.out = 300)
  tibble(id = id, x = c(xs, rev(xs)), y = c(dnorm(xs), rep(0, length(xs))))
}
above <- function(a, b, id) {
  xs <- seq(a, b, length.out = 300)
  tibble(id = id, x = c(xs, rev(xs)), y = c(dnorm(xs), rep(peak, length(xs))))
}

# Scale towards grow_from, so the slice emerges from a single point
shrink <- function(df, g) {
  mutate(
    df,
    x = grow_from[["x"]] + g * (x - grow_from[["x"]]),
    y = grow_from[["y"]] + g * (y - grow_from[["y"]])
  )
}

arrow_both <- arrow(ends = "both", length = unit(0.18, "cm"), type = "closed")
bell <- tibble(x = seq(-zmax, zmax, length.out = 801), y = dnorm(x))
baseline <- tibble(x = c(-zmax, zmax), y = 0)

## ----- Frame -----------------------------------------------------------------
draw_frame <- function(f) {
  p <- ggplot() +
    coord_cartesian(
      xlim = c(-zmax, zmax),
      ylim = c(y_lo, peak + 0.01),
      expand = FALSE,
      clip = "off"
    ) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 12, 4, 4)
    )
  if (f$grow == 0) {
    return(p)
  }

  ## ----- Fills ---------------------------------------------------------------
  # Orange under the rectangle once it is full, so no white shows through the
  # anti-aliased edges as the corners turn orange.
  if (f$pour == 1) {
    p <- p +
      annotate(
        "rect",
        xmin = -half_width,
        xmax = half_width,
        ymin = 0,
        ymax = peak,
        fill = orange
      )
  }
  tail_col <- blend_hcl(orange, turquoise, f$recolour)
  corner_col <- blend_hcl(turquoise, orange, f$merge)
  body <- under(-f$x_drain, f$x_drain, "body")
  tails <- bind_rows(
    under(-f$x_drain, -half_width, "left"),
    under(half_width, f$x_drain, "right")
  )
  p <- p +
    geom_polygon(
      data = shrink(body, f$grow),
      aes(x, y, group = id),
      fill = orange
    )
  if (f$recolour > 0 && f$pour < 1) {
    p <- p + geom_polygon(data = tails, aes(x, y, group = id), fill = tail_col)
  }
  if (f$pour > 0) {
    corners <- bind_rows(
      above(-half_width, -f$x_fill, "left"),
      above(f$x_fill, half_width, "right")
    )
    p <- p +
      geom_polygon(data = corners, aes(x, y, group = id), fill = corner_col)
  }

  ## ----- Lines and labels ----------------------------------------------------
  p <- p +
    geom_line(
      data = shrink(baseline, f$grow),
      aes(x, y),
      colour = col_axis,
      alpha = pmin(1, 3 * f$grow)
    ) +
    geom_line(
      data = shrink(bell, f$grow),
      aes(x, y),
      colour = col_curve,
      linewidth = 0.7,
      alpha = pmin(1, 3 * f$grow) * (1 - (1 - curve_fade) * f$merge)
    ) +
    annotate(
      "text",
      x = zmax,
      y = y_width,
      label = "vartheta[-j]",
      parse = TRUE,
      colour = col_axis,
      size = 6,
      hjust = 1,
      alpha = f$label
    )
  if (f$width > 0) {
    p <- p +
      annotate(
        "rect",
        xmin = -half_width * f$width,
        xmax = half_width * f$width,
        ymin = 0,
        ymax = peak,
        fill = NA,
        colour = col_ink,
        linewidth = 0.7
      )
  }
  # Arrows wait until the shaft is longer than its two heads
  if (f$width > 0.12) {
    p <- p +
      annotate(
        "segment",
        x = -half_width * f$width,
        xend = half_width * f$width,
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
        size = 6.5,
        alpha = f$width
      )
  }
  if (f$height > 0.08) {
    p <- p +
      annotate(
        "segment",
        x = 0,
        xend = 0,
        y = 0,
        yend = f$height * (peak - 0.012),
        colour = col_ink,
        linewidth = 0.7,
        arrow = arrow_both
      ) +
      annotate(
        "text",
        x = -0.08,
        y = 0.3 * peak,
        label = "height",
        colour = col_ink,
        size = 6.5,
        hjust = 1,
        alpha = f$height
      )
  }
  if (f$dot > 0) {
    p <- p +
      annotate("point", x = 0, y = peak, colour = col_peak, size = 2.6 * f$dot)
  }
  p
}

## ----- Render ----------------------------------------------------------------
# Same size as the static figure (fig-width 5.4, fig-height 3.5). Loops, and
# the slide restarts it from the beginning when its fragment shows.
png_dir <- tempfile("slice_section_")
dir.create(png_dir)
png_files <- file.path(png_dir, sprintf("frame%04d.png", seq_len(nrow(frames))))
for (i in seq_len(nrow(frames))) {
  ragg::agg_png(
    png_files[i],
    width = 5.4,
    height = 3.5,
    units = "in",
    res = 160
  )
  print(draw_frame(frames[i, ]))
  dev.off()
}
gifski::gifski(
  png_files,
  gif_file,
  width = 5.4 * 160,
  height = 3.5 * 160,
  delay = 1 / fps,
  loop = TRUE
)
unlink(png_dir, recursive = TRUE)
