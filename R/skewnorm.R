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
p_skewnorm <- p2 / p1 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")

## ----- Points only, before the fit -------------------------------------------
# The evaluated log-density ordinates alone. The density panel and legend stay
# as invisible placeholders so this frame lines up with the full figure above.
p2_points <- p2
p2_points$layers <- discard(p2$layers, \(l) inherits(l$geom, "GeomLine"))
p2_points$data <- filter(p2$data, type == "Evaluated")
p1_ghost <- p1
p1_ghost$layers <- list()
p1_ghost <- p1_ghost +
  geom_blank(aes(shape = type), show.legend = FALSE) +  # trains shape scale
  theme(
    axis.text = element_text(colour = "transparent"),
    strip.text = element_text(colour = "transparent"),
    panel.grid = element_line(colour = "transparent")
  )
p_skewnorm_points <- p2_points / p1_ghost + plot_layout(guides = "collect") &
  guides(shape = guide_legend(override.aes = list(alpha = 0))) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(colour = "transparent")
  )
