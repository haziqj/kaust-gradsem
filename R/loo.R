library(INLAvaan)

# LOO model comparison for the Kidney example (refit-free, via INLAvaan::loo()).
# Kidney is simulated from the two-factor model below (see R/kidney.R), so the base
# model is the true one; the comparison checks that LOO confirms it against two
# plausible rivals.
# Run from the project root.

source("R/kidney.R")  # Kidney data

# Base: the model fitted throughout the talk
mod_base <- "
  GlyCon =~ y1 + y2 + y3
  KdnHlt =~ y4 + y5 + y6
  KdnHlt ~ GlyCon
"

# Rival 1: a single 'metabolic' factor behind all six markers (simpler)
mod_1f <- "
  Metab =~ y1 + y2 + y3 + y4 + y5 + y6
"

# Rival 2: base + fasting insulin also reflecting kidney function (richer).
# Biologically motivated -- the kidneys clear insulin -- but not in the data.
mod_xload <- "
  GlyCon =~ y1 + y2 + y3
  KdnHlt =~ y4 + y5 + y6 + y3
  KdnHlt ~ GlyCon
"

fit_kidney <- function(model) asem(model, Kidney, meanstructure = TRUE, std.lv = TRUE)
fit_base  <- fit_kidney(mod_base)
fit_1f    <- fit_kidney(mod_1f)
fit_xload <- fit_kidney(mod_xload)

# Per-model LOO (ELPD, p_loo, LOOIC)
loo(fit_base)

# Paired comparison: first argument is the baseline
compare(fit_base, fit_1f, fit_xload, loo = TRUE)
