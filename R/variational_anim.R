library(tidyverse)
library(sn)
library(gganimate)

# Animated version of R/variational.R: the VB-corrected Gaussian (turquoise)
# peels away from the Laplace Gaussian (orange) at the mode and slides to the
# mean, while the overlap with the truth and the mean shift delta update.
# Run from the project root; writes figures/variational.gif.

# Same setup as R/variational.R
xi <- 0; omega <- 1; alpha <- 8
dskew <- function(x) dsn(x, xi, omega, alpha)
mode_skew <- optimize(function(x) -1 * dskew(x), c(-5, 5))$minimum
delta <- alpha / sqrt(1 + alpha^2)
mean_skew <- xi + omega * sqrt(2 / pi) * delta
var_skew <- omega^2 * (1 - 2 * delta^2 / pi)
sd_laplace <- sqrt(var_skew) / 1.1

col_lap <- "#e38b0b"
col_vbc <- "#00A6AA"

## ---- Overlap ------------------------------------------------------------------
# Overlap = area shared by the two densities = 1 - total variation distance,
# i.e. int min(p, q) = 1 - 1/2 int |p - q|  (in [0, 1]; higher = closer)
overlap <- function(centre) {
  integrate(function(x) pmin(dskew(x), dnorm(x, centre, sd_laplace)), -Inf, Inf)$value
}

## ---- Frames -------------------------------------------------------------------
# Ease-in-out so the Gaussian accelerates away from the mode and settles on the mean
n_move <- 40  # frames of motion
n_hold <- 6   # extra frames at the mean, during which the Mean line finishes fading in
ease <- function(t) ifelse(t < 0.5, 4 * t^3, 1 - (-2 * t + 2)^3 / 2)
frames <-
  tibble(
    frame = seq_len(n_move + n_hold),
    progress = c(ease(seq(0, 1, length.out = n_move)), rep(1, n_hold)),
    centre = mode_skew + (mean_skew - mode_skew) * progress
  ) |>
  mutate(ovl = map_dbl(centre, overlap))

# Mean line fades in: starts once the shift is 95% complete, fully opaque by the last frame
fade_start <- min(which(frames$progress >= 0.95))
frames <- mutate(frames, mean_alpha = pmin(1, pmax(0, (frame - fade_start) / (max(frame) - fade_start))))

x <- seq(-1.5, 3, length.out = 600)
truth_df <- tibble(x = x, Density = dskew(x))
lap_df <- tibble(x = x, Density = dnorm(x, mode_skew, sd_laplace))
vbc_df <-
  frames |>
  select(frame, centre) |>
  crossing(x = x) |>
  mutate(Density = dnorm(x, centre, sd_laplace))
# Arrow and delta appear only once there is room for them
arrow_df <- filter(frames, centre - mode_skew > 0.06)
delta_df <- filter(frames, centre - mode_skew > 0.15)
mean_df <- filter(frames, mean_alpha > 0)

## ---- Plot ---------------------------------------------------------------------
p <-
  ggplot() +
  geom_area(data = truth_df, aes(x, Density), fill = "gray") +
  # VB curve drawn first so the Laplace curve is on top: it starts orange and peels away
  geom_line(data = vbc_df, aes(x, Density, group = frame), colour = col_vbc, linewidth = 1.1) +
  geom_line(data = lap_df, aes(x, Density), colour = col_lap, linewidth = 1.1) +
  geom_vline(xintercept = mode_skew, linetype = "dashed", colour = col_lap, linewidth = 1.1) +
  geom_vline(data = mean_df, aes(xintercept = mean_skew, alpha = mean_alpha),
             linetype = "dashed", colour = col_vbc, linewidth = 1.1) +
  annotate("label", x = -0.75, y = 0.45, colour = col_lap, size = 6.5,
           label = scales::percent(frames$ovl[1], accuracy = 0.1)) +
  geom_label(data = frames, aes(x = 1.80, y = 0.45, label = scales::percent(ovl, accuracy = 0.1)),
             colour = col_vbc, size = 6.5) +
  annotate("text", mode_skew, 0.85, label = "Mode ", hjust = 1, size = 6, colour = col_lap) +
  geom_text(data = mean_df, aes(x = mean_skew, y = 0.85, alpha = mean_alpha),
            label = " Mean", hjust = 0, size = 6, colour = col_vbc) +
  geom_segment(data = arrow_df, aes(x = mode_skew, xend = centre, y = 0.85, yend = 0.85),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), colour = "black") +
  geom_text(data = delta_df, aes(x = (mode_skew + centre) / 2, y = 0.85), label = "delta",
            parse = TRUE, vjust = -0.3, size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
  scale_alpha_identity() +
  theme_void() +
  transition_manual(frame)

# Same aspect as the static figure (fig-width 7, fig-height 2.8)
animate(
  p,
  nframes = n_move + n_hold + 15 + 50, fps = 20,  # nframes includes the pauses
  start_pause = 15, end_pause = 50,
  width = 7, height = 2.8, units = "in", res = 160,
  renderer = gifski_renderer("figures/variational.gif", loop = TRUE)
)
