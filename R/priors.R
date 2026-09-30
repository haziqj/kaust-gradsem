library(tidyverse)
library(INLAvaan)

# Effect of strong priors on the Kidney example, fitted as in the rest of the talk
# (meanstructure = TRUE, std.lv = TRUE, so the latent variances are fixed at 1).
#
# The Kidney indicators live on very different scales, so one prior for all loadings
# cannot make sense. Strong priors are applied only to indicators on a comparable scale,
# y1, y2 and y6 (all other parameters keep the default priors):
#   - loadings (true 0.87, 1.1, 1.8):             N(4, 0.1^2)
#   - residual variances (true 0.33, 0.70, 3.1):  Gamma(100, 50)  (mean 2, SD 0.2)
# Each panel has its own strong-prior fit, so the two sets of priors do not interact.
# Run from the project root; writes R/p_priors.RData (p1: loadings, p2: residual variances).

source("R/kidney.R")  # Kidney data

mod <- "
  # Measurement model
  GlyCon =~ y1 + y2 + y3
  KdnHlt =~ y4 + y5 + y6

  # Structural model
  KdnHlt ~ GlyCon
"

mod_lam <- "
  GlyCon =~ prior('normal(2,0.1)')*y1 + prior('normal(2,0.1)')*y2 + y3
  KdnHlt =~ y4 + y5 + prior('normal(2,0.1)')*y6
  KdnHlt ~ GlyCon
"

mod_theta <- "
  GlyCon =~ y1 + y2 + y3
  KdnHlt =~ y4 + y5 + y6
  KdnHlt ~ GlyCon
  y1 ~~ prior('gamma(100,50)')*y1
  y2 ~~ prior('gamma(100,50)')*y2
  y6 ~~ prior('gamma(100,50)')*y6
"

fit_kidney <- function(model) asem(model, Kidney, meanstructure = TRUE, std.lv = TRUE)
fit0      <- fit_kidney(mod)        # default priors
fit_lam   <- fit_kidney(mod_lam)    # strong priors on three loadings
fit_theta <- fit_kidney(mod_theta)  # strong priors on three residual variances

# Select parameters by name: the fits order their parameters differently
pdf_of <- function(fit) INLAvaan:::get_inlavaan_internal(fit)$pdf_data
nm_lam   <- c("GlyCon=~y1", "GlyCon=~y2", "KdnHlt=~y6")
nm_theta <- c("y1~~y1", "y2~~y2", "y6~~y6")

# Posterior densities of the selected parameters under both priors
pdf_df <- function(nms, fit_strong) {
  bind_rows(
    list(
      weak = bind_rows(pdf_of(fit0)[nms], .id = "param"),
      strong = bind_rows(pdf_of(fit_strong)[nms], .id = "param")
    ),
    .id = "prior"
  ) |>
    mutate(param = factor(param, levels = nms))
}

p1 <-
  pdf_df(nm_lam, fit_lam) |>
  ggplot(aes(x, y, col = prior)) +
  geom_line() +
  facet_wrap(~param, nrow = 1, scales = "free") +
  theme_minimal() +
  scale_colour_manual(values = c("#00A6AA", "#e38b0b")) +
  labs(x = NULL, y = "Posterior density")

p2 <-
  pdf_df(nm_theta, fit_theta) |>
  ggplot(aes(x, y, col = prior)) +
  geom_line() +
  facet_wrap(~param, nrow = 1, scales = "free") +
  theme_minimal() +
  scale_colour_manual(values = c("#00A6AA", "#e38b0b")) +
  labs(x = NULL, y = "Posterior density")

save(p1, p2, file = "R/p_priors.RData")
