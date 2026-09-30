library(tidyverse)
library(INLAvaan)
library(patchwork)
source("R/kidney.R")

mod <- "
  # Measurement model
  GlyCon =~ y1 + y2 + y3
  KdnHlt =~ y4 + y5 + y6

  # Structural model
  KdnHlt ~ GlyCon
"
tmp <- capture.output({
  fit <<- asem(mod, Kidney[1:24, ], debug = TRUE) |>
    suppressMessages() |>
    suppressWarnings()
})

p1 <- INLAvaan:::visual_debug(fit, params = "y1~~y1", logscale = FALSE, points = TRUE, raw = FALSE) +
  annotate(
    "text",
    x = -4,
    y = 1 + 0.01,
    label = "Density",
    size = 3.5,
    hjust = -0.025,
    vjust = 0.15
  ) + coord_cartesian(ylim = c(0,1.05))
p2 <- INLAvaan:::visual_debug(fit, params = "y1~~y1", logscale = TRUE, points = TRUE, raw = FALSE) +
  annotate(
    "text",
    x = -4,
    y = 0.1,
    label = "Log density",
    size = 3.5,
    hjust = -0.025,
    vjust = 0.15
  ) + coord_cartesian(ylim = c(0.2,-4.1))
p2 / p1 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")
