library(tidyverse)
library(INLAvaan)

# Per-person latent variables for the Kidney example, fitted as in the rest of the
# talk (meanstructure = TRUE, std.lv = TRUE).
#
# predict(fit, type = "lv") returns posterior draws of eta_i for every individual.
# Each person is summarised by their posterior mean (dot) and a 95% posterior
# ellipse (halo; Gaussian, chi^2_2 radius). The red line is the structural slope
# beta from the SEM itself (posterior draws), not a regression on the dots.
# Three slides build up the picture:
#   p_lv1: posterior means only
#   p_lv2: + 95% halos, three people outlined
#   p_lv3: halos + structural regression line with 95% credible band (no outlines)
# Run from the project root; writes R/p_latent.RData (p_lv1, p_lv2, p_lv3, and
# lv_stats with the numbers quoted on the slide).

source("R/kidney.R")  # Kidney data

mod <- "
  # Measurement model
  GlyCon =~ y1 + y2 + y3
  KdnHlt =~ y4 + y5 + y6

  # Structural model
  KdnHlt ~ GlyCon
"
fit <- asem(mod, Kidney, meanstructure = TRUE, std.lv = TRUE)

# Posterior draws of the latent variables -------------------------------------
set.seed(221)
eta_samp <- predict(fit, type = "lv", nsamp = 1000)
arr <- simplify2array(lapply(eta_samp, as.matrix))  # n x 2 x nsamp
n <- dim(arr)[1]

mu_df <- tibble(
  id = seq_len(n),
  x = apply(arr[, 1, ], 1, mean),
  y = apply(arr[, 2, ], 1, mean)
)

r95 <- sqrt(qchisq(0.95, df = 2))
tt <- seq(0, 2 * pi, length.out = 72)
circ <- rbind(cos(tt), sin(tt))
ell <- map_dfr(seq_len(n), function(i) {
  L <- t(chol(cov(t(arr[i, , ]))))
  pts <- t(c(mu_df$x[i], mu_df$y[i]) + r95 * L %*% circ)
  tibble(id = i, x = pts[, 1], y = pts[, 2])
})

# Three people to show individual halos: two at either end of eta1 (not
# outliers), and one off the trend in the upper left, clear of the other two
hl <- map_int(
  quantile(mu_df$x, c(0.06, 0.94)),
  ~ which.min(abs(mu_df$x - .x))
)
hl <- c(hl, which.min((mu_df$x + 1)^2 + (mu_df$y - 1.3)^2))

# Structural slope ------------------------------------------------------------
beta <- sampling(fit, nsamp = 1000)[, "KdnHlt~GlyCon"]
band <- tibble(x = seq(-3.2, 3.2, length.out = 100)) |>
  mutate(
    fit = mean(beta) * x,
    lo = map_dbl(x, ~ quantile(beta * .x, 0.025)),
    hi = map_dbl(x, ~ quantile(beta * .x, 0.975))
  )

# Plots -----------------------------------------------------------------------
col_pt <- "#131516"
col_halo <- "grey60"   # neutral gray
col_hl <- "#e38b0b"    # orange
col_line <- "#b10f2e"  # KAUST red

base <- ggplot() +
  coord_equal(xlim = c(-3.2, 3.2), ylim = c(-3.4, 3.6), expand = FALSE) +
  labs(
    x = expression("Glycemic control" ~ (eta[1])),
    y = expression("Kidney dysfunction" ~ (eta[2]))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA)
  )

lyr_halo <- list(
  geom_polygon(
    data = ell, aes(x, y, group = id),
    fill = col_halo, alpha = 0.05, colour = NA
  )
)
lyr_pts <- list(
  geom_point(data = mu_df, aes(x, y), colour = col_pt, size = 1, alpha = 0.75)
)
lyr_hl <- list(
  geom_polygon(
    data = filter(ell, id %in% hl), aes(x, y, group = id),
    fill = NA, colour = col_hl, linewidth = 0.8
  ),
  geom_point(data = filter(mu_df, id %in% hl), aes(x, y), colour = col_hl, size = 2.6)
)
lyr_line <- list(
  geom_ribbon(data = band, aes(x, ymin = lo, ymax = hi), fill = col_line, alpha = 0.2),
  geom_line(data = band, aes(x, fit), colour = col_line, linewidth = 1)
)

# Labels go on top of everything, on an opaque box so the halo doesn't show through
lbl_hl <- annotate(
  "label", x = mu_df$x[hl[2]], y = min(ell$y[ell$id == hl[2]]) - 0.08,
  label = "95% posterior region", colour = col_hl, size = 4, vjust = 1,
  fill = "white", border.colour = NA
)
lbl_line <- annotate(
  "label", x = 3.15, y = mean(beta) * 3.15 + 0.5, hjust = 1,
  label = sprintf("beta == %.2f", mean(beta)), parse = TRUE,
  colour = col_line, size = 4.2, fill = "white", border.colour = NA
)

p_lv1 <- base + lyr_pts
p_lv2 <- base + lyr_halo + lyr_pts + lyr_hl + lbl_hl
p_lv3 <- base + lyr_halo + lyr_line + lyr_pts + lbl_line

lv_stats <- list(
  beta = mean(beta),
  beta_lo = unname(quantile(beta, 0.025)),
  beta_hi = unname(quantile(beta, 0.975)),
  post_sd = apply(apply(arr, c(1, 2), sd), 2, median)  # typical posterior SD, per factor
)

save(p_lv1, p_lv2, p_lv3, lv_stats, file = "R/p_latent.RData")
