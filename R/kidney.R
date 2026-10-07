library(lavaan)
library(INLAvaan)

# Song & Lee (2004) inspired example: glycemic control and kidney health
#
# Observed variables (y1-y6), in the units shown on the slides:
#   y1: HbA1c (%)                            -- mean ~6.5, SD ~1.21
#   y2: Fasting plasma glucose (mmol/L)      -- mean ~6.1, SD ~1.23
#   y3: Fasting insulin level (µU/mL)        -- mean ~12,  SD ~5.6
#   y4: Plasma creatinine (µmol/L)           -- mean ~88,  SD ~29
#   y5: Albumin-creatinine ratio (mg/g)      -- mean ~60,  SD ~26
#   y6: Blood urea nitrogen (mmol/L)         -- mean ~7.9, SD ~2.35
#
# Latent factors (eta1 var = 1; eta2 residual var = 1, as in a std.lv = TRUE fit):
#   eta1: Glycemic Control   (higher = worse glycaemia)
#   eta2: Kidney Dysfunction (higher = more dysfunction)
#
# var(yi) = lambda_i^2 * var(eta) + theta_i, with theta_i set so that every
# indicator has reliability (R^2) 0.8.
# Structural path 0.60: var(eta2) = 1 + 0.6^2 = 1.36, so eta1 explains ~26%
# of the variance in eta2 (standardised coefficient ~0.51).

n <- 250

mod_tru <- "
  # -- Measurement model --
  # eta1: Glycemic Control
  eta1 =~ 1.08*y1 + 1.1*y2 + 5*y3

  # eta2: Kidney Dysfunction
  eta2 =~ 22*y4 + 20*y5 + 1.8*y6

  # -- Structural model --
  eta2 ~ 0.60*eta1

  # -- Latent (residual) variances --
  eta1 ~~ 1*eta1
  eta2 ~~ 1*eta2

  # -- Residual variances (theta) --
  # theta = 0.25 * lambda^2 * var(eta), so reliability = 1 / 1.25 = 0.8
  # (var(eta1) = 1, var(eta2) = 1.36)
  y1 ~~ 0.29*y1   # 0.25 * 1.08^2
  y2 ~~ 0.30*y2   # 0.25 * 1.1^2
  y3 ~~ 6.25*y3   # 0.25 * 5^2
  y4 ~~  165*y4   # 0.25 * 22^2  * 1.36
  y5 ~~  136*y5   # 0.25 * 20^2  * 1.36
  y6 ~~  1.1*y6   # 0.25 * 1.8^2 * 1.36

  # -- Intercepts --
  y1 ~ 6.5*1
  y2 ~ 6.1*1
  y3 ~  12*1
  y4 ~  88*1
  y5 ~  60*1
  y6 ~ 7.9*1
"

Kidney <- simulateData(mod_tru, sample.nobs = n, seed = 300926)

# # Fitted model (first indicator scaling; intercepts free)
# mod <- "
#   # Measurement model
#   GlyCon =~ y1 + y2 + y3
#   KdnHlt =~ y4 + y5 + y6
# 
#   # Structural model
#   KdnHlt ~ GlyCon
# "
# 
# fit <- asem(mod, Kidney)

