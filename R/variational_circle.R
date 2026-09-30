library(tidyverse)
library(sn)

# 1. Parameters for the skew normal
xi <- 0 # location
omega <- 1 # scale
alpha <- 8 # skewness (make it really skewed)

# 2. Functions for skew normal
dskew <- function(x) 2.5*dsn(x, xi, omega, alpha) 
pskew <- function(x) psn(x, xi, omega, alpha)

# 3. Compute mode and mean of the skew normal
# Mode: numerical optimization
mode_skew <- optimize(function(x) -1 * dskew(x), c(-5, 5))$minimum
delta <- alpha / sqrt(1 + alpha^2)
mean_skew <- xi + omega * sqrt(2 / pi) * delta

# 4. Laplace: normal at mode, variance = theoretical variance of skew normal
# Analytical variance of skew normal
var_skew <- omega^2 * (1 - 2 * delta^2 / pi)
sd_laplace <- sqrt(var_skew)
# var_skew <- 1 / numDeriv::hessian(function(x) -1 * dskew(x), mode_skew)
# sd_laplace <- sqrt(as.numeric(var_skew))

## ---- Jensen-Shannon Error ---------------------------------------------------
# 1. Define a helper function for KL Divergence
# KL(P || Q) = Integral[ P(x) * log(P(x) / Q(x)) ]
calc_kl_base2 <- function(p_func, q_func, lower = -Inf, upper = Inf) {
  integrand <- function(x) {
    px <- p_func(x)
    qx <- q_func(x)
    term <- ifelse(px < 1e-10, 0, px * log(px / (qx + 1e-20), base = 2))
    return(term)
  }
  integrate(integrand, lower, upper)$value
}

# 2. Define a helper function for Jensen-Shannon Divergence
# JSD(P || Q) = 0.5 * KL(P || M) + 0.5 * KL(Q || M)
# where M = 0.5 * (P + Q)
calc_js_base2 <- function(p_func, q_func) {
  m_func <- function(x) 0.5 * (p_func(x) + q_func(x))
  0.5 * calc_kl_base2(p_func, m_func) + 0.5 * calc_kl_base2(q_func, m_func)
}

# 3. Define the density functions for your specific approximations
f_tru <- function(x) dskew(x)
f_lap <- function(x) dnorm(x, mean = mode_skew, sd = sd_laplace)
f_vbc <- function(x) dnorm(x, mean = mean_skew, sd = sd_laplace)

# 4. Compute the JS Errors
js_error_lap <- calc_js_base2(f_tru, f_lap)
js_error_vbc <- calc_js_base2(f_tru, f_vbc)

# 5. Display Results
js_df <- tibble(
  Type = c("Laplace", "VB Correction"),
  JS_div = c(js_error_lap, js_error_vbc),
  JS_dis = sqrt(JS_div),
  JS_sim = 1 - JS_dis
)

## ----- Plot ------------------------------------------------------------------
cor_scl <- dskew(mode_skew) / dnorm(mode_skew, mode_skew, sd_laplace)

plot_df <-
  tibble(
    x = seq(-1.5, 3, length.out = 1000),
    Truth = dskew(x),
    Laplace = cor_scl*dnorm(x, mean = mode_skew, sd = sd_laplace),
    `VB Correction` = dnorm(x, mean = mean_skew, sd = sd_laplace)
  ) |>
  pivot_longer(-x, names_to = "Type", values_to = "Density") |>
  mutate(Type = factor(Type, levels = c("Truth", "Laplace", "VB Correction")))

ggplot(plot_df, aes(x, Density, fill = Type)) +
  geom_area(position = "identity") +
  geom_line(
    data = filter(plot_df, Type != "Truth"),
    aes(x, Density, color = Type),
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = mode_skew,
    linetype = "dashed",
    color = "#e38b0b",
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = mean_skew,
    linetype = "dashed",
    color = "#00A6AA",
    linewidth = 1.1
  ) +
  annotate(
    "label",
    x = -0.75,
    y = 0.45,
    label = scales::percent(js_df$JS_sim[1], accuracy = 0.1),
    col = "#e38b0b",
    size = 6.5
  ) +
  annotate(
    "label",
    x = 1.80,
    y = 0.45,
    label = scales::percent(js_df$JS_sim[2], accuracy = 0.1),
    col = "#00A6AA",
    size = 6.5
  ) +
  annotate(
    "text",
    mode_skew,
    0.85,
    label = "Mode ",
    hjust = 1,
    size = 6,
    col = "#e38b0b"
  ) +
  annotate(
    "text",
    mean_skew,
    0.85,
    label = " Mean",
    hjust = 0,
    size = 6,
    col = "#00A6AA"
  ) +
  annotate(
    "segment",
    x = mode_skew,
    xend = mean_skew,
    y = 0.85,
    yend = 0.85,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    color = "black" # Using black to distinguish structural elements from data
  ) +
  annotate(
    "text",
    x = (mode_skew + mean_skew) / 2,
    y = 0.85,
    label = expression(delta),
    vjust = -0.3, # Pushes the text slightly below the arrow
    size = 6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
  scale_colour_manual(values = c("#e38b0b", "#00A6AA")) +
  scale_fill_manual(values = c("gray", NA, NA)) +
  theme_void() +
  theme(legend.position = "none")

## ----- Step by step ----------------------------------------------------------

# Just the true density
p0 <-
  filter(plot_df, Type == "Truth") |> 
  ggplot(aes(x, Density)) +
  geom_area(fill = "gray") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
  theme_void() +
  coord_equal(xlim = c(-1.5, 3), ylim = c(0, 2.15))

#  True density + vertical line at mode
p1 <-
  p0 +
  geom_vline(xintercept = mode_skew, linetype = "dashed", color = "#e38b0b", linewidth = 1) +
  annotate("text", mode_skew, 2.1, label = "Mode ", hjust = 1, size = 6, col = "#e38b0b") 

# p1 + a circle tangent to the peak from below. radius = sd_laplace (inverse curvature).
# Center is shifted down by sd_laplace so the top of the circle kisses the mode.
radius <- sd_laplace^2 / dskew(mode_skew)
circle_df <- tibble(
  t = seq(0, 2 * pi, length.out = 300),
  x = mode_skew + radius * cos(t),
  y = (dskew(mode_skew) - radius) + radius * sin(t)
) 
p2 <-
  p1 +
  geom_path(data = circle_df, aes(x, y), color = "#e38b0b", linewidth = 1,
            inherit.aes = FALSE) +
  annotate("segment",
           x = mode_skew, xend = mode_skew + radius,
           y = dskew(mode_skew) - radius, yend = dskew(mode_skew) - radius,
           color = "#e38b0b", linewidth = 0.6) +
  annotate("text",
           x = 0.3+mode_skew + radius / 2,
           y = 0.1+dskew(mode_skew) - radius,
           label = "r == 1/kappa",
           parse = TRUE,
           vjust = -0.4,
           color = "#e38b0b",
           size = 5); p2

# p1 + gaussian density at the mode
p3 <-
  p2 +
  geom_line(data = filter(plot_df, Type == "Laplace"), aes(x, Density), color = "#e38b0b", linewidth = 1) +
  annotate("text", -0.7, 0.85, label = "Laplace", hjust = 1, size = 6, col = "#e38b0b"); p3

save(p0, p1, p2, p3, file = "R/p_laplace_circle.RData")
