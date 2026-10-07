library(tidyverse)
library(patchwork)

# Is the height term enough? Two-panel figures (joint | marginal of theta_j):
#   - Gaussian: CMP = conditional mode, and the width is constant along it, so the
#     normalised height alone already equals the true marginal.
#   - Trumpet (a mild Neal's funnel): CMP = conditional mode again, but the width
#     changes along it, so the height alone is shifted; the Lemma 2 volume
#     correction (height x width) recovers the true marginal.
# gamma'(0) is the slope of -1/2 log|H_{-j}| along the CMP at the mode (Lemma 2),
# estimated numerically exactly as the method would.
# Run from the project root; writes figures/cmp_gaussian.png and figures/cmp_trumpet.png.

col_cmp   <- "#00A6AA"  # CMP (KAUST turquoise)
col_true  <- "#E07B00"  # true conditional-mode path (KAUST orange, darkened)
col_h     <- "#b10f2e"  # height only (KAUST red, as the spine in R/slices.R)
col_hw    <- "#5284C4"  # height x width (KAUST blue, as the marginal in R/slices.R)
col_width <- "gray30"   # conditional width intervals

theme_panel <- function() {
  theme_classic(base_size = 13) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "gray40", linewidth = 0.4,
                               arrow = arrow(length = unit(2, "mm"), type = "closed")),
      axis.title = element_text(colour = "gray25"),
      plot.title = element_text(face = "bold", size = 13, colour = "gray20"),
      plot.title.position = "plot"
    )
}

# Vertical legend, shown between the two panels. Items sit left of centre (cx) so the
# extra white space on the right visually attaches the legend to the joint panel.
cx <- 0.36; half <- 0.26
p_legend <-
  ggplot() +
  annotate("segment", x = cx - half, xend = cx + half, y = 0.90, yend = 0.90,
           colour = col_cmp, linewidth = 1.3) +
  annotate("text", cx, 0.82, label = "CMP", colour = col_cmp, size = 4.2) +
  annotate("segment", x = cx - half, xend = cx + half, y = 0.60, yend = 0.60, colour = col_true,
           linewidth = 0.9, linetype = "dashed") +
  annotate("text", cx, 0.49, label = "conditional\nmode", colour = col_true, size = 4.2,
           lineheight = 0.9) +
  annotate("segment", x = cx, xend = cx, y = 0.16, yend = 0.28, colour = col_width, linewidth = 0.5,
           arrow = arrow(ends = "both", length = unit(0.08, "cm"), angle = 90)) +
  annotate("text", cx, 0.08, label = "\u00b11 SD", colour = col_width, size = 4.2) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void()

two_panel <- function(logf, start, xl, yl, width_at, true_marg, marg_lim, ymax,
                      levels_drop = c(0.5, 2, 4.5, 8, 12.5), labels) {
  # Joint mode, Laplace covariance, CMP for the marginal of theta_j (= theta_1)
  m <- optim(start, function(t) -logf(t), method = "BFGS", control = list(reltol = 1e-12))$par
  Omega <- solve(-numDeriv::hessian(logf, m))
  cmp <- function(t1) m[2] + Omega[2, 1] / Omega[1, 1] * (t1 - m[1])
  H22 <- function(t1) -numDeriv::hessian(logf, c(t1, cmp(t1)))[2, 2]

  # Lemma 2: slope of -1/2 log|H_{-j}| along the CMP at the mode
  eps <- 1e-3
  gamma1 <- (-0.5 * log(H22(m[1] + eps)) + 0.5 * log(H22(m[1] - eps))) / (2 * eps)

  # Left: contours, conditional widths, CMP and true conditional mode
  g1 <- seq(xl[1], xl[2], length.out = 400)
  g2 <- seq(yl[1], yl[2], length.out = 400)
  z <- outer(g1, g2, Vectorize(function(u, v) logf(c(u, v))))
  cl <- contourLines(g1, g2, z, levels = logf(m) - rev(levels_drop))
  contours <- imap_dfr(cl, ~ tibble(t1 = .x$x, t2 = .x$y, id = .y))
  # +/- 1 conditional SD of theta_{-j} along the CMP
  widths <- tibble(t1 = width_at, hw = 1 / sqrt(map_dbl(width_at, H22)), mid = cmp(width_at))
  cond <- tibble(t1 = g1, t2 = map_dbl(g1, function(a) {
    optimize(function(b) logf(c(a, b)), yl + c(-20, 20), maximum = TRUE)$maximum
  }))

  p_left <-
    ggplot() +
    geom_path(data = contours, aes(t1, t2, group = id), colour = "gray78", linewidth = 0.45) +
    geom_segment(data = widths, aes(x = t1, xend = t1, y = mid - hw, yend = mid + hw),
                 colour = col_width, linewidth = 0.5,
                 arrow = arrow(ends = "both", length = unit(0.08, "cm"), angle = 90)) +
    geom_line(data = tibble(t1 = g1, t2 = cmp(g1)), aes(t1, t2), colour = col_cmp, linewidth = 1.3) +
    geom_line(data = cond, aes(t1, t2), colour = col_true, linewidth = 0.9, linetype = "dashed") +
    annotate("point", m[1], m[2], size = 2.2) +
    coord_cartesian(xlim = xl, ylim = yl, expand = FALSE, clip = "off") +
    labs(title = "Joint", x = expression(vartheta[j]), y = expression(vartheta[-j])) +
    theme_panel()

  # Right: marginal of theta_j -- truth, normalised height, height x width
  x <- seq(marg_lim[1], marg_lim[2], length.out = 600)
  normalise <- function(y) y / sum(y * diff(x)[1])
  height <- exp(map_dbl(x, ~ logf(c(.x, cmp(.x)))) - logf(m))
  marg <- tibble(x = x, truth = true_marg(x), height_only = normalise(height),
                 corrected = normalise(height * exp(gamma1 * (x - m[1]))))

  p_right <-
    ggplot(marg, aes(x)) +
    geom_area(aes(y = truth), fill = "gray82") +
    geom_line(aes(y = corrected), colour = col_hw, linewidth = 1.2) +
    geom_line(aes(y = height_only), colour = col_h, linewidth = 1.1, linetype = "dashed") +
    labels$right +
    coord_cartesian(xlim = marg_lim, ylim = c(0, ymax), expand = FALSE, clip = "off") +
    labs(title = "Marginal", x = expression(vartheta[j]), y = "normalised height") +
    theme_panel()

  message(sprintf("mode = (%.2f, %.2f); gamma'(0) = %.3f", m[1], m[2], gamma1))
  # middle panel doubles as a vertical legend
  p_left + p_legend + p_right + plot_layout(widths = c(1, 0.45, 1))
}

## ---- 1. Bivariate normal (rho = 0.5): constant width, height alone suffices -----
rho <- 0.5
xl_g <- c(-3.2, 3.2); yl_g <- c(-3.2, 3.2)  # joint-panel ranges
Sinv <- solve(matrix(c(1, rho, rho, 1), 2))
logf_gauss <- function(t) -0.5 * drop(t %*% Sinv %*% t)

p_gauss <- two_panel(
  logf_gauss, start = c(0.1, 0.1), xl = xl_g, yl = yl_g,
  width_at = c(-2, -1, 0, 1, 2), true_marg = dnorm, marg_lim = c(-3.5, 3.5), ymax = 0.47,
  levels_drop = c(0.5, 2, 4.5, 8),
  labels = list(
    right = list(
      annotate("text", -1.35, 0.40, label = "height only", colour = col_h, size = 4.2, hjust = 1),
      annotate("text", 1.35, 0.40, label = "height × width", colour = col_hw, size = 4.2, hjust = 0),
      annotate("text", 3.4, 0.03, label = "truth", colour = "gray45", size = 4, hjust = 1)
    )
  )
)

## ---- 2. Trumpet (mild Neal's funnel, sheared): width changes, height alone is shifted
a <- 1.2
# Sheared so its ridge (CMP = conditional mode) sits at the same on-screen angle as the
# Gaussian's: theta_{-j} | theta_j ~ N(c theta_j, exp(a theta_j)). The shear only slides each
# vertical slice, so the widths, gamma'(0) and the marginal of theta_j are unchanged.
# c accounts for the different panel ranges, so the two CMPs look parallel.
xl_t <- c(-3, 3); yl_t <- c(-6, 6)  # joint-panel ranges
c_shear <- rho * (diff(yl_t) / diff(xl_t)) / (diff(yl_g) / diff(xl_g))
logf_trumpet <- function(t) -t[1]^2 / 2 - a * t[1] / 2 - (t[2] - c_shear * t[1])^2 * exp(-a * t[1]) / 2

p_trumpet <- two_panel(
  logf_trumpet, start = c(0.1, 0.1), xl = xl_t, yl = yl_t,
  width_at = c(-2, -1, 0, 1, 2), true_marg = dnorm, marg_lim = c(-3.5, 3.5), ymax = 0.47,
  labels = list(
    right = list(
      annotate("text", -1.55, 0.40, label = "height only", colour = col_h, size = 4.2, hjust = 1),
      annotate("text", 1.15, 0.40, label = "height × width", colour = col_hw, size = 4.2, hjust = 0),
      annotate("text", 3.4, 0.03, label = "truth", colour = "gray45", size = 4, hjust = 1)
    )
  )
)

ggsave("figures/cmp_gaussian.png", p_gauss, width = 9.8, height = 3.2, dpi = 200, bg = "white")
ggsave("figures/cmp_trumpet.png", p_trumpet, width = 9.8, height = 3.2, dpi = 200, bg = "white")
