library(lavaan)
library(INLAvaan)
# devtools::install_github("haziqj/brlavaan")
library(tidyverse)
library(patchwork)
library(sn)
library(furrr)
future::plan(multisession, workers = future::availableCores() - 1)
source("R/kidney.R")

# Data generation function with true parameter values
gen_data <- function(n = 100) {
  true_model <- "
    # -- Measurement model --
    # GlyCon: Glycemic control; HbA1c (%), FPG (mmol/L), insulin (µU/mL)
    GlyCon =~ 0.87*y1 + 1.1*y2 + 5*y3

    # KdnHlt: Kidney dysfunction; PCr (µmol/L), ACR (mg/g), BUN (mmol/L)
    KdnHlt =~ 22*y4 + 20*y5 + 1.8*y6

    # -- Structural model --
    KdnHlt ~ 0.60*GlyCon   # standardised ~0.51; R^2 ~26%

    # -- Latent variances --
    # Both (residual) variances fixed to 1, matching the std.lv = TRUE fit, so
    # these values are the true parameters of the fitted model
    GlyCon ~~ 1*GlyCon
    KdnHlt ~~ 1*KdnHlt

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
  lavaan::simulateData(true_model, sample.nobs = n)
}

mod <- "
  GlyCon =~ y1 + y2 + y3
  KdnHlt =~ y4 + y5 + y6

  KdnHlt ~ GlyCon
"

## COMPARE MCMC
dat <- gen_data(n = 250)
library(blavaan)
fit_blav <- bsem(
  mod, dat, meanstructure = TRUE, std.lv = TRUE, 
  bcontrol = list(cores = 3)
)
fit_inlv <- asem(mod, dat, meanstructure = TRUE, std.lv = TRUE)

p_compare <- INLAvaan:::compare_mcmc(fit_blav, INLAvaan = fit_inlv)$p_compare
save(p_compare, file = "R/p_compare.RData")


# True parameter values (for comparison)
true_params <- tribble(
  ~parameter     , ~true_value ,
  # Intercepts
  "y1~1"         , 6.5         ,
  "y2~1"         , 6.1         ,
  "y3~1"         , 12          ,
  "y4~1"         , 88          ,
  "y5~1"         , 60          ,
  "y6~1"         , 7.9         ,
  # Factor loadings
  "GlyCon=~y1"   , 0.87        ,
  "GlyCon=~y2"   , 1.1         ,
  "GlyCon=~y3"   , 5           ,
  "KdnHlt=~y4"   , 22          ,
  "KdnHlt=~y5"   , 20          ,
  "KdnHlt=~y6"   , 1.8         ,
  # Regressions
  "KdnHlt~GlyCon", 0.60        ,
  # Residual variances
  "y1~~y1"       , 0.33        ,
  "y2~~y2"       , 0.70        ,
  "y3~~y3"       , 11          ,
  "y4~~y4"       , 477         ,
  "y5~~y5"       , 225         ,
  "y6~~y6"       , 3.1
)

# Function to run one simulation
run_simulation <- function(n, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  dat <- gen_data(n = n)
  tmp <- capture.output(
    fit <<- suppressMessages(
      asem(mod, data = dat, std.lv = TRUE, meanstructure = TRUE, debug = TRUE)
    )
  )
  m <- length(coef(fit))

  estimates <- fit$summary |>
    rownames_to_column("parameter") |>
    select(parameter, Mean, SD, `2.5%`, `25%`, `75%`, `97.5%`)
  
  # Extract pdf_data for plotting
  pdf_data <- fit$pdf_data |>
    map_dfr(~ as.data.frame(.x), .id = "parameter")

  # Calculate log score (posterior density at true value)
  log_score_data <-
    pdf_data |>
    left_join(true_params, by = "parameter") |>
    group_by(parameter) |>
    summarise(
      log_score = {
        true_val <- true_value[1]
        approx_val <- approx(x, y, xout = true_val)$y
        ifelse(is.na(approx_val), NA_real_, log(pmax(approx_val, 1e-10)))
      },
      .groups = "drop"
    ) |>
    select(parameter, log_score)

  # Combine with true values
  results <- estimates |>
    left_join(true_params, by = "parameter") |>
    left_join(log_score_data, by = "parameter") |>
    mutate(
      n = n,
      seed = seed,
      bias = Mean - true_value,
      rel_bias = bias / true_value,
      coverage_95 = (true_value >= `2.5%`) & (true_value <= `97.5%`),
      coverage_50 = (true_value >= `25%`) & (true_value <= `75%`),
      ci_95_length = `97.5%` - `2.5%`,
      ci_50_length = `75%` - `25%`
    )

  # Store pdf_data as a list column
  results <- results |>
    left_join(
      pdf_data |>
        group_by(parameter) |>
        summarise(pdf = list(data.frame(x = x, y = y)), .groups = "drop"),
      by = "parameter"
    )

  return(results)
}
run_simulation_possibly <- possibly(run_simulation, otherwise = NULL)

# Run simulations across different sample sizes
sample_sizes <- c(50, 100, 250, 500, 1000)
n_replications <- 250 # Number of replications per sample size

# Run simulation study
cat("Running simulation study...\n")
tictoc::tic()
simulation_results <-
  map_df(sample_sizes, function(n) {
    cat(sprintf("\nSample size n = %d\n", n))
    future_map_dfr(
      1:n_replications,
      function(rep) run_simulation_possibly(n = n, seed = 2026 + rep),
      .id = "replication",
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
  })
tictoc::toc()

# Save results
# save(simulation_results, file = "R/A-simulation_results.RData")

# Summary statistics
summary_stats <-
  simulation_results |>
  group_by(n, parameter) |>
  summarise(
    mean_estimate = mean(Mean, na.rm = TRUE),
    sd_estimate = sd(Mean, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    rmse = sqrt(mean(bias^2, na.rm = TRUE)),
    coverage_rate_95 = mean(coverage_95, na.rm = TRUE),
    coverage_rate_50 = mean(coverage_50, na.rm = TRUE),
    mean_ci_95_length = mean(ci_95_length, na.rm = TRUE),
    mean_ci_50_length = mean(ci_50_length, na.rm = TRUE),
    mean_log_score = mean(log_score, na.rm = TRUE),
    true_value = first(true_value),
    .groups = "drop"
  )

# Print summary for key parameters
key_params <- c(
  "y3~1", # insulin intercept
  "GlyCon=~y3", # fasting insulin loading
  "KdnHlt~GlyCon", # structural path
  "y4~~y4", # creatinine residual variance
  "y2~~y2", # fasting glucose residual variance
  "y5~1" # albumin-creatinine ratio intercept
)

cat("\n\nSummary for key parameters:\n")
summary_stats |>
  filter(parameter %in% key_params) |>
  arrange(parameter, n)

# Visualization function: Plot posteriors across replications using pdf_data
plot_posterior_recovery <- function(
  params = NULL,
  results = simulation_results,
  show_facet_labels = FALSE
) {
  results_n <- results
  # Filter for specific parameters if provided
  if (!is.null(params)) {
    results_n <- results_n |> filter(parameter %in% params)
  } else {
    results_n <- filter(results_n, parameter %in% key_params)
  }

  # Extract and unnest pdf_data
  density_data <-
    results_n |>
    select(
      replication,
      parameter,
      n,
      true_value,
      log_score,
      coverage_95,
      coverage_50,
      pdf
    ) |>
    unnest(pdf) |>
    mutate(
      n_label = factor(
        paste0("n = ", n),
        levels = paste0("n = ", sample_sizes)
      ),
      parameter_label = factor(
        parameter,
        levels = unique(parameter)
      )
    )

  # Create the plot
  p <-
    density_data |>
    ggplot(aes(x = x, y = y, group = replication)) +
    geom_line(alpha = 0.11, linewidth = 0.35, color = "#00A6AA") +
    geom_vline(
      aes(xintercept = true_value),
      color = "#131516",
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    # Add mean log score annotations
    geom_text(
      data = density_data |>
        group_by(parameter, parameter_label, n, n_label) |>
        summarise(
          mean_log_score = mean(log_score, na.rm = TRUE),
          coverage_95 = mean(coverage_95, na.rm = TRUE) * 100,
          coverage_50 = mean(coverage_50, na.rm = TRUE) * 100,
          .groups = "drop"
        ) |>
        mutate(
          label = paste0(
            # "MLS = ",
            # sprintf("%.2f", mean_log_score),
            # "\n",
            "CR95 = ",
            sprintf("%.0f", coverage_95),
            "%\n",
            "CR50 = ",
            sprintf("%.0f", coverage_50),
            "%"
          )
        ),
      aes(
        x = Inf,
        y = Inf,
        label = label
      ),
      hjust = 1.05,
      vjust = 1.3,
      size = 2.5,
      color = "black",
      inherit.aes = FALSE
    ) +
    facet_grid(
      n_label ~ parameter_label,
      scales = "free"
    ) +
    theme_minimal() +
    labs(
      x = "Parameter value",
      y = "Posterior density"
    ) +
    theme(
      strip.text.x = element_text(size = 9),
      strip.text.y = if (show_facet_labels) {
        element_text(size = 9)
      } else {
        element_blank()
      },
      axis.text = element_text(size = 7)
    )

  return(p)
}

p1 <- plot_posterior_recovery("y3~1") + coord_cartesian(xlim = c(6, 22))
p2 <- plot_posterior_recovery("GlyCon=~y3") + coord_cartesian(xlim = c(1.2, 12))
p3 <- plot_posterior_recovery("KdnHlt=~y6") + coord_cartesian(xlim = c(0, 4.5))
p4 <- plot_posterior_recovery("KdnHlt~GlyCon") + coord_cartesian(xlim = c(-0.5, 3))
p5 <- plot_posterior_recovery("y3~~y3", show_facet_labels = TRUE) + coord_cartesian(xlim = c(0, 40))


p1 + p2 + p3 + p4 + p5 + 
  plot_layout(axes = "collect", guides = "collect", nrow = 1)

save(p1, p2, p3, p4, p5, file = "R/p_simulations.RData")
