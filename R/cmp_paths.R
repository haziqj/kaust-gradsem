library(tidyverse)

# How well does the conditional mean path (CMP) track the "height" term of the
# marginal Laplace approximation? For the marginal of theta_1 we need, at each
# theta_1, the conditional mode theta_2*(theta_1) = argmax_{theta_2} pi(theta_1, theta_2).
# The CMP replaces this optimisation by the straight line implied by the Gaussian
# (Laplace) approximation at the joint mode:
#   theta_2 = theta_2* + (Omega_21 / Omega_11) (theta_1 - theta_1*).
# It is exact for a Gaussian, but can be far off for non-Gaussian shapes.
#
# Everything is computed in the natural coordinates; only the picture is rotated
# (theme_void, no axes) so that the CMP sits at 45 degrees.
# Run from the project root; writes figures/cmp_funnel.png and figures/cmp_banana.png.
# (The Gaussian figure is now two-panel, in R/cmp_trumpet.R.)

col_cmp  <- "#00A6AA"  # CMP (turquoise)
col_true <- "#b10f2e"  # true conditional-mode path (KAUST red)
levels_drop <- c(0.5, 2, 4.5, 8, 12.5)  # log-density contours at f* - c ("1 to 5 SDs" for a Gaussian)

rotate <- function(df, angle) {
  # rotate columns t1, t2 counter-clockwise by `angle` (radians) into x, y
  mutate(df, x = cos(angle) * t1 - sin(angle) * t2, y = sin(angle) * t1 + cos(angle) * t2)
}

cmp_paths <- function(logf, start, t1_lim, t2_lim, n_grid = 400) {
  # Joint mode and Laplace covariance
  m <- optim(start, function(t) -logf(t), method = "BFGS", control = list(reltol = 1e-12))$par
  Omega <- solve(-numDeriv::hessian(logf, m))
  fmax <- logf(m)

  # Log-density contours (as paths, so they can be rotated)
  g1 <- seq(t1_lim[1], t1_lim[2], length.out = n_grid)
  g2 <- seq(t2_lim[1], t2_lim[2], length.out = n_grid)
  z <- outer(g1, g2, Vectorize(function(a, b) logf(c(a, b))))
  cl <- contourLines(g1, g2, z, levels = fmax - rev(levels_drop))
  contours <- imap_dfr(cl, ~ tibble(t1 = .x$x, t2 = .x$y, id = .y))

  # Paths for the marginal of theta_1: CMP (straight line) and true conditional mode
  t1 <- seq(t1_lim[1], t1_lim[2], length.out = 300)
  cmp <- tibble(t1 = t1, t2 = m[2] + Omega[2, 1] / Omega[1, 1] * (t1 - m[1]))
  cond <- tibble(t1 = t1, t2 = map_dbl(t1, function(a) {
    optimize(function(b) logf(c(a, b)), t2_lim + c(-20, 20), maximum = TRUE)$maximum
  }))

  list(mode = tibble(t1 = m[1], t2 = m[2]), Omega = Omega,
       contours = contours, cmp = cmp, cond = cond,
       cmp_angle = atan(Omega[2, 1] / Omega[1, 1]))
}

# angle = on-screen angle of the CMP, in degrees
draw_cmp <- function(res, xlim, ylim, angle = 45, rot = angle * pi / 180 - res$cmp_angle) {
  ggplot() +
    geom_path(data = rotate(res$contours, rot), aes(x, y, group = id),
              colour = "gray55", linewidth = 0.5) +
    # geom_path (not geom_line) keeps the data order after rotation
    geom_path(data = rotate(res$cmp, rot), aes(x, y), colour = col_cmp, linewidth = 1.4) +
    geom_path(data = rotate(res$cond, rot), aes(x, y), colour = col_true,
              linewidth = 1.1, linetype = "dashed") +
    geom_point(data = rotate(res$mode, rot), aes(x, y), size = 2.5) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void()
}

## ---- 2. Neal's funnel: CMP is way off --------------------------------------------
# v ~ N(0, 3^2), x | v ~ N(0, e^v). theta_1 = x (the marginal we want), theta_2 = v.
# The Laplace covariance at the mode (the neck) is diagonal, so the CMP is flat,
# while the conditional mode climbs up into the bowl of the funnel.
logf_funnel <- function(t) -t[2]^2 / 18 - t[2] / 2 - t[1]^2 / (2 * exp(t[2]))

res_funnel <- cmp_paths(logf_funnel, start = c(0.1, -4), t1_lim = c(-12, 12), t2_lim = c(-11, 8))
p_funnel <-
  draw_cmp(res_funnel, xlim = c(-8, 6), ylim = c(-8, 8)) +
  annotate("text", 3.4, 0.2, label = "CMP", colour = col_cmp, size = 5, hjust = 0) +
  annotate("text", -7.5, 0.3, label = "conditional\nmode", colour = col_true, size = 5,
           hjust = 0, lineheight = 0.9)

## ---- 3. Banana: CMP is the tangent to a curved ridge, then drifts off ------------
# theta_1 ~ N(0, 1.5^2), theta_2 | theta_1 ~ N(b theta_1^2, 0.5^2). The conditional
# mode is the parabola theta_2 = b theta_1^2; the CMP is its tangent at the mode.
b <- 0.4
logf_banana <- function(t) -t[1]^2 / (2 * 1.5^2) - (t[2] - b * t[1]^2)^2 / (2 * 0.5^2)

res_banana <- cmp_paths(logf_banana, start = c(0.1, 0.1), t1_lim = c(-12, 12), t2_lim = c(-8, 16))
p_banana <-
  draw_cmp(res_banana, xlim = c(-8, 7), ylim = c(-2.2, 3.2), angle = 30) +
  annotate("text", 5.6, 2.3, label = "CMP", colour = col_cmp, size = 5, hjust = 0) +
  annotate("text", -7.8, 0.6, label = "conditional\nmode", colour = col_true, size = 5,
           hjust = 0, lineheight = 0.9)

# Landscape (banana): shown under the lemma box on the Step 3 slide
ggsave("figures/cmp_funnel.png", p_funnel, width = 4.5, height = 5.5, dpi = 200, bg = "white")
ggsave("figures/cmp_banana.png", p_banana, width = 8, height = 2.9, dpi = 200, bg = "white")
