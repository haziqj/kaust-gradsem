library(tidyverse)
library(ggridges)
library(mvtnorm)
kaust_cols <- c(
  # purple_sat = "#9C6FAE",
  # blue_sat = "#5284C4",
  turquoise_sat = "#00A6AA",
  green_sat = "#CDCE00",
  yellow_sat = "#F0B500",
  orange_sat = "#F18F00"
)

# 1. Setup Data
zmax <- 3.5
x1_slices <- seq(-zmax, zmax, by = 0.25)
x2_path <- seq(-zmax, zmax, length.out = 500)
rho <- 0.5
scale_factor <- 30 # Define scale once so we can reuse it for the line

plot_data <- expand.grid(x1 = x1_slices, x2 = x2_path)

# 2. Calculate Density
# We remove the factor conversion for x1 to keep plotting coordinates numeric
plot_data <- plot_data %>%
  mutate(
    density = map2_dbl(
      x1,
      x2,
      ~ dmvnorm(c(.x, .y), mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), 2))
    )
  )

# 3. Create the "Spine" Data (The Peak Trace)
# The peak of a slice occurs exactly at the Conditional Mean: E[x2 | x1] = rho * x1
peak_data <- tibble(x1 = x1_slices) %>%
  mutate(
    x2 = rho * x1, # The location of the peak on the x-axis
    # Calculate height at the peak
    peak_height = map2_dbl(
      x1,
      x2,
      ~ dmvnorm(c(.x, .y), mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), 2))
    ),
    # Calculate the visual Y position: Base + (Height * Scale)
    y_plot = x1 + (peak_height * scale_factor)
  )

# 4. Plot
ggplot() +

  # Layer 1: The Ridges
  # We use 'group = x1' because x1 is now continuous
  geom_ridgeline_gradient(
    data = plot_data,
    aes(x = x2, y = x1, height = density, fill = density, group = x1),
    scale = scale_factor,
    color = "white",
    alpha = 0.3
  ) +

  # Layer 2: The Trace Line (The Spine)
  geom_line(
    data = peak_data,
    aes(x = x2, y = y_plot),
    color = "#b10f2e",
    size = 1,
    # linetype = "dashed"
  ) +

  # # Layer 3: The Peak Points
  # geom_point(
  #   data = peak_data,
  #   aes(x = x2, y = y_plot),
  #   color = "black",
  #   size = 1
  # ) +

  # # Layer 4: Height indicator at y = -1
  # geom_segment(
  #   data = peak_data %>% filter(x1 == -1),
  #   aes(x = x2, xend = x2, y = x1, yend = y_plot),
  #   color = "#b10f2e",
  #   size = 1,
  #   arrow = arrow(ends = "both", length = unit(0.15, "cm"))
  # ) +
  #
  # # Layer 5: SD indicator at y = -1 (conditional SD is sqrt(1 - rho^2))
  # geom_segment(
  #   data = peak_data %>% filter(x1 == -1) %>%
  #     mutate(
  #       cond_sd = sqrt(1 - rho^2),
  #       y_peak = y_plot
  #     ),
  #   aes(x = x2 - cond_sd, xend = x2 + cond_sd, y = y_peak, yend = y_peak),
  #   color = "#1f77b4",
  #   size = 1,
  #   arrow = arrow(ends = "both", length = unit(0.15, "cm"))
  # ) +

  scale_fill_gradientn(colors = kaust_cols) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none")
