library(tidyverse)
library(sn)
library(patchwork)

# Laplace as a quadratic fit on the log scale: the osculating circle at the mode
# has radius 1/H = Omega (the Laplace variance) when the axes are equally scaled.
# Run from the project root. Used on the [Step 1] Laplace (cont.) slide.

# Skew-normal "posterior" (unnormalised; same as R/variational_circle.R)
xi <- 0; omega <- 1; alpha <- 8
dskew <- function(x) 2.5 * dsn(x, xi, omega, alpha)
f <- function(x) log(dskew(x))  # log-posterior

# Laplace: mode and curvature of the LOG-density
mode_skew <- optimize(function(x) -dskew(x), c(-5, 5))$minimum
fmax <- f(mode_skew)
H <- as.numeric(-numDeriv::hessian(f, mode_skew))  # curvature = -f''
Omega <- 1 / H                                      # Laplace variance
quad <- function(x) fmax - H / 2 * (x - mode_skew)^2
lap  <- function(x) exp(quad(x))                    # = dskew(mode) * exp(-H (x - mode)^2 / 2)

col_lap <- "#e38b0b"
col_rad <- "#b10f2e"  # circle / radius (KAUST red)
lab_size <- 4          # small labels: Mode, quadratic, exp(quadratic)

## ---- Left: log scale (equal axes, so the circle is a circle) -----------------
# Panel sizing: the left panel has equal axes, so its y-range is set to fill the
# full figure height at its share of the width (fig size = slide chunk: 8 x 3.5 in).
fig_w <- 8; fig_h <- 3.5
widths <- c(1.1, 1.45)
xl <- c(-0.6, 1.1)
left_w <- fig_w * widths[1] / sum(widths)
yl <- c(fmax + 0.35 - diff(xl) * fig_h / left_w, fmax + 0.35)
log_df <- tibble(x = seq(xl[1], xl[2], length.out = 600), f = f(x), q = quad(x))

L0 <-
  ggplot(log_df, aes(x, f)) +
  geom_line(colour = "gray45", linewidth = 1.2) +
  geom_vline(xintercept = mode_skew, linetype = "dashed", colour = col_lap, linewidth = 1) +
  annotate("text", mode_skew, yl[2] - 0.03, label = " Mode", hjust = 0, vjust = 1,
           size = lab_size, colour = col_lap) +
  annotate("text", xl[1], yl[2], label = "'log scale:'~f(vartheta)", parse = TRUE,
           hjust = 0, vjust = 1, size = 5, colour = "gray30") +
  coord_equal(xlim = xl, ylim = yl, expand = FALSE, clip = "on") +
  theme_void()

circle_df <- tibble(
  t = seq(0, 2 * pi, length.out = 300),
  x = mode_skew + Omega * cos(t),
  y = (fmax - Omega) + Omega * sin(t)
)
L1 <-
  L0 +
  geom_path(data = circle_df, aes(x, y), colour = col_rad, linewidth = 1, inherit.aes = FALSE) +
  annotate("segment", x = mode_skew, xend = mode_skew + Omega,
           y = fmax - Omega, yend = fmax - Omega, colour = col_rad, linewidth = 0.6) +
  # small "r" on the radius; its definition sits in a free corner (same colour, no leader line)
  annotate("text", x = mode_skew + Omega / 2, y = fmax - Omega, label = "r",
           vjust = -0.3, size = 3.5, colour = col_rad) +
  annotate("text", x = mode_skew + Omega + 0.05, y = fmax + 0.07,
           label = "paste(r, ' = 1/', H, ' = ', Omega)", parse = TRUE,
           hjust = 0, size = lab_size, colour = col_rad)

L2 <-
  L1 +
  geom_line(aes(y = q), colour = col_lap, linewidth = 1, linetype = "solid") +
  annotate("text", x = mode_skew - sqrt(2 * 0.75 / H) - 0.04, y = fmax - 0.75,
           label = "quadratic", hjust = 1, size = lab_size, colour = col_lap)

## ---- Right: density scale -------------------------------------------------------
dens_df <- tibble(x = seq(-1.5, 3, length.out = 1000), Truth = dskew(x), Laplace = lap(x))

R0 <-
  ggplot(dens_df, aes(x, Truth)) +
  geom_area(fill = "gray") +
  geom_vline(xintercept = mode_skew, linetype = "dashed", colour = col_lap, linewidth = 1) +
  annotate("text", mode_skew, 2.1, label = " Mode", hjust = 0, size = lab_size, colour = col_lap) +
  annotate("text", 3, 2.1, label = "'density:'~pi(vartheta~'|'~bold(y))", parse = TRUE,
           hjust = 1, vjust = 1, size = 5, colour = "gray30") +
  coord_cartesian(xlim = c(-1.5, 3), ylim = c(0, 2.15), expand = FALSE) +
  theme_void()

exp_lab_y <- 1.5
exp_lab_x <- uniroot(function(x) dskew(x) - (exp_lab_y - 0.15), c(mode_skew, 3))$root + 0.04

R3 <-
  R0 +
  geom_line(aes(y = Laplace), colour = col_lap, linewidth = 1) +
  # just clear of the grey true density: start where the truth drops below the
  # label (with a little padding), at a height close to the Gaussian
  annotate("text", x = exp_lab_x, y = exp_lab_y,
           label = "exp(quadratic)", hjust = 0, size = lab_size, colour = col_lap)

side <- function(l, r) l + r + plot_layout(widths = widths)
p0 <- side(L0, R0)
p1 <- side(L1, R0)
p2 <- side(L2, R0)
p3 <- side(L2, R3)

save(p0, p1, p2, p3, file = "R/p_laplace_logscale.RData")
