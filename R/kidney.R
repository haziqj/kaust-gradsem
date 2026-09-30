library(lavaan)
library(INLAvaan)

# Song & Lee (2004) inspired example: glycemic control and kidney health
#
# Observed variables (y1-y6), in the units shown on the slides:
#   y1: HbA1c (%)                            -- mean ~6.5, SD ~1.04
#   y2: Fasting plasma glucose (mmol/L)      -- mean ~6.1, SD ~1.38
#   y3: Fasting insulin level (µU/mL)        -- mean ~12,  SD ~6
#   y4: Plasma creatinine (µmol/L)           -- mean ~88,  SD ~34
#   y5: Albumin-creatinine ratio (mg/g)      -- mean ~60,  SD ~28
#   y6: Blood urea nitrogen (mmol/L)         -- mean ~7.9, SD ~2.74
#
# Latent factors (eta1 var = 1; eta2 residual var = 1, as in a std.lv = TRUE fit):
#   eta1: Glycemic Control   (higher = worse glycaemia)
#   eta2: Kidney Dysfunction (higher = more dysfunction)
#
# Loadings chosen so that var(yi) = lambda_i^2 * var(eta) + theta_i
# giving R^2 ≈ 0.58-0.71 for each indicator.
# Structural path 0.60: var(eta2) = 1 + 0.6^2 = 1.36, so eta1 explains ~26%
# of the variance in eta2 (standardised coefficient ~0.51).

n <- 250

mod_tru <- "
  # -- Measurement model --
  # eta1: Glycemic Control
  eta1 =~ 0.87*y1 + 1.1*y2 + 5*y3

  # eta2: Kidney Dysfunction
  eta2 =~ 22*y4 + 20*y5 + 1.8*y6

  # -- Structural model --
  eta2 ~ 0.60*eta1

  # -- Latent (residual) variances --
  eta1 ~~ 1*eta1
  eta2 ~~ 1*eta2

  # -- Residual variances (theta) --
  # var(yi) = lambda^2 + theta  =>  theta = var(yi) - lambda^2
  # (y4-y6 set assuming var(KdnHlt) = 1; it is now 1.36, so their SDs are
  # ~34, ~28, ~2.74 and R^2 ~0.58, 0.71, 0.59)
  y1 ~~ 0.33*y1   # 1.04^2 - 0.87^2
  y2 ~~ 0.70*y2   # 1.38^2 - 1.1^2
  y3 ~~   11*y3   #    6^2 -   5^2
  y4 ~~  477*y4   #   31^2 -  22^2
  y5 ~~  225*y5   #   25^2 -  20^2
  y6 ~~  3.1*y6   # 2.52^2 - 1.8^2

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

