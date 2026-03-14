# Meropenem in CRRT patients
# two-compartment model with linear elimination (iv)
# C++ analytical solution (no rxode2 dependency)

## ----- Loading library ----
library(pacman)
pacman::p_load(admr, Rcpp, tidyverse, MASS, data.table, pracma,
               writexl, scales, future, future.apply, progressr, parallel)

# Enable progress handlers globally
handlers(global = TRUE)
handlers("progress")

# Sets up parallel processing — cap at 8 workers to avoid memory/compilation issues
n_workers <- min(parallel::detectCores() - 1, 40)
plan(multisession, workers = n_workers)

## ----- Build & install C++ as a package ----
pkg_dir <- file.path(getwd(), "twocompiv")
# Only build if not already installed
if (!requireNamespace("twocompiv", quietly = TRUE)) {
  # Create package skeleton
  if (dir.exists(pkg_dir)) unlink(pkg_dir, recursive = TRUE)
  Rcpp::Rcpp.package.skeleton(
    name = "twocompiv",
    path = getwd(),
    cpp_files = character(0),
    example_code = FALSE
  )
  # Write the C++ source
  writeLines('
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix Cpp_twocomp_iv(NumericVector cl, NumericVector v1,
                             NumericVector v2, NumericVector q,
                             NumericVector ti) {
  int n_ind  = cl.size();
  int n_time = ti.size();
  NumericMatrix out(n_ind, n_time);
  for (int iter = 0; iter < n_ind; ++iter) {
    double k   = cl[iter] / v1[iter];
    double k12 = q[iter]  / v1[iter];
    double k21 = q[iter]  / v2[iter];
    double sum_rates = k + k12 + k21;
    double sqrt_disc = sqrt(sum_rates * sum_rates - 4.0 * k * k21);
    double alpha = 0.5 * (sum_rates + sqrt_disc);
    double beta  = 0.5 * (sum_rates - sqrt_disc);
    double A = (alpha - k21) / (alpha - beta) / v1[iter];
    double B = (k21 - beta)  / (alpha - beta) / v1[iter];
    for (int i = 0; i < n_time; ++i) {
      out(iter, i) = A * exp(-alpha * ti[i]) + B * exp(-beta * ti[i]);
    }
  }
  return out;
}
', file.path(pkg_dir, "src", "twocomp_iv.cpp"))
  # Clean up auto-generated files
  unlink(file.path(pkg_dir, "src", "rcpp_hello_world.cpp"), force = TRUE)
  unlink(file.path(pkg_dir, "src", "*.o"), force = TRUE)
  unlink(file.path(pkg_dir, "src", "*.so"), force = TRUE)
  unlink(file.path(pkg_dir, "man"), recursive = TRUE, force = TRUE)
  dir.create(file.path(pkg_dir, "man"), showWarnings = FALSE)
  # Generate Rcpp bindings (creates RcppExports.cpp and RcppExports.R)
  Rcpp::compileAttributes(pkg_dir)
  # Write NAMESPACE with explicit export
  writeLines(c(
    'useDynLib(twocompiv, .registration=TRUE)',
    'importFrom(Rcpp, evalCpp)',
    'export(Cpp_twocomp_iv)'
  ), file.path(pkg_dir, "NAMESPACE"))
  # Install the package
  install.packages(pkg_dir, repos = NULL, type = "source")
}
library(twocompiv)

## ----- Define Simulation Conditions ----

n_subjects   <- c(20, 40, 75, 100, 250, 500, 750, 1000, 2500, 10000)
n_timepoints <- 12
n_reps       <- 50

get_timepoints <- function(n) {
  switch(as.character(n),
         "4"  = c(0.1, 1, 6, 12),
         "5"  = c(0.1, 0.75, 2, 6, 12),
         "6"  = c(0.1, 0.5, 2, 4, 8, 12),
         "8"  = c(0.1, 0.5, 1, 2, 4, 6, 8, 12),
         "12" = c(0.08, 0.25, 0.5, 1, 1.5, 2, 3, 5, 7, 9, 10, 12),
         stop("Unsupported number of timepoints"))
}

times <- get_timepoints(n_timepoints)

## ----- Define Global True Parameters  ----

TRUE_PARAMS <- list(
  theta      = c(cl = 3.03, v1 = 14.2, v2 = 14.6, q = 15.9),
  iiv_sd     = c(cl = 0.52, v1 = 0.33, v2 = 0.33, q = 0.26),
  Sigma_prop = 0.153,
  dose       = 1000
)

TRUE_OMEGA <- diag(TRUE_PARAMS$iiv_sd^2)

## ----- Prediction wrapper ----

# Wrapper matching admr predder interface: returns matrix [n_ind x n_time]
pk_predder <- function(time, theta_i, dose = TRUE_PARAMS$dose) {
  theta_i <- as.matrix(theta_i)
  n_ind   <- nrow(theta_i)
  if (n_ind == 0) {
    theta_i <- matrix(theta_i, nrow = 1)
    n_ind   <- 1
  }
  doses <- rep(dose, n_ind)
  pred  <- twocompiv::Cpp_twocomp_iv(theta_i[, 1], theta_i[, 2],
                                     theta_i[, 3], theta_i[, 4], time)
  pred * doses
}

## ----- Compute TRUE_VALS using C++ ----

fine_times <- seq(0.001, 12, 0.01)
true_pk    <- as.numeric(
  pk_predder(fine_times, matrix(TRUE_PARAMS$theta, nrow = 1))
)

omega_names     <- c("omega_cl", "omega_v1", "omega_v2", "omega_q")
TRUE_OMEGA_DIAG <- setNames(diag(TRUE_OMEGA), omega_names)

TRUE_VALS <- c(
  TRUE_PARAMS$theta,
  auc            = trapz(fine_times, true_pk),
  cmax           = max(true_pk),
  residual_error = TRUE_PARAMS$Sigma_prop^2,
  TRUE_OMEGA_DIAG
)


## ----- Main Simulation (Parallel Execution) -----

extract_results <- function(fit, true_vals) {
  tryCatch({
    tp <- fit$transformed_params

    pk_params <- c("cl", "v1", "v2", "q")
    est_beta  <- tp$beta[pk_params]

    # AUC and Cmax from estimated beta using C++
    out_fine <- as.numeric(
      pk_predder(fine_times, matrix(est_beta, nrow = 1))
    )
    auc  <- trapz(fine_times, out_fine)
    cmax <- max(out_fine)

    re_est     <- as.numeric(tp$Sigma_prop)
    omega_diag <- setNames(diag(tp$Omega)[1:4], omega_names)

    res_vec <- c(est_beta, auc = auc, cmax = cmax,
                 residual_error = re_est, omega_diag)
    res_vec[names(true_vals)]

  }, error = function(e) NULL)
}

run_single_iteration <- function(n_sub, i) {
  tryCatch({
    set.seed(n_sub * 1000 + i)

    ## ---- Generate data ----
    mv <- MASS::mvrnorm(n_sub, rep(0, nrow(TRUE_OMEGA)), TRUE_OMEGA)
    params_all <- data.table(
      ID = 1:n_sub,
      cl = TRUE_PARAMS$theta['cl'] * exp(mv[, 1]),
      v1 = TRUE_PARAMS$theta['v1'] * exp(mv[, 2]),
      v2 = TRUE_PARAMS$theta['v2'] * exp(mv[, 3]),
      q  = TRUE_PARAMS$theta['q']  * exp(mv[, 4])
    )

    # Simulate individual concentrations using C++
    theta_mat <- as.matrix(params_all[, .(cl, v1, v2, q)])
    sim_mat   <- pk_predder(times, theta_mat)  # [n_sub x n_time]

    # Apply floor and use as DV
    sim_mat <- pmax(sim_mat, 0.01)

    # sim_mat is already [n_sub x n_time] — use directly as ind_wide
    ind_wide <- sim_mat

    agg_full <- admr::meancov(ind_wide)
    diag(agg_full$V) <- diag(agg_full$V) + (TRUE_PARAMS$Sigma_prop^2) * (agg_full$E^2)

    agg_var   <- agg_full
    agg_var$V <- diag(diag(agg_var$V))

    ## ----- Fitting model -----
    predder <- function(time, theta_i) {
      theta_i <- as.matrix(theta_i)
      n_ind   <- nrow(theta_i)
      if (n_ind == 0) {
        theta_i <- matrix(theta_i, nrow = 1)
        n_ind   <- 1
      }
      twocompiv::Cpp_twocomp_iv(theta_i[, 1], theta_i[, 2],
                                theta_i[, 3], theta_i[, 4], time) * TRUE_PARAMS$dose
    }

    base_opts <- list(
      time = times,
      p    = list(beta = TRUE_PARAMS$theta, Omega = TRUE_OMEGA,
                  Sigma_prop = TRUE_PARAMS$Sigma_prop),
      nsim = 5000, n = n_sub, fo_appr = FALSE,
      omega_expansion = 1, f = predder
    )

    opts_covar <- do.call(admr::genopts, c(base_opts, list(no_cov = FALSE)))
    opts_var   <- do.call(admr::genopts, c(base_opts, list(no_cov = TRUE)))

    fit_covar <- tryCatch(
      admr::fitIRMC(opts = opts_covar, obs = agg_full,
                    chains = 1, maxiter = 8000, use_grad = TRUE),
      error = function(e) NULL)

    fit_var <- tryCatch(
      admr::fitIRMC(opts = opts_var, obs = agg_var,
                    chains = 1, maxiter = 8000, use_grad = TRUE),
      error = function(e) NULL)

    if (is.null(fit_covar) || is.null(fit_var)) return(NULL)

    cov_res <- extract_results(fit_covar, TRUE_VALS)
    var_res <- extract_results(fit_var,   TRUE_VALS)

    if (is.null(cov_res) || is.null(var_res)) return(NULL)

    data.frame(
      N_Subject   = n_sub,
      Rep         = i,
      Metric      = names(TRUE_VALS),
      True_Val    = as.numeric(TRUE_VALS),
      Covar_Est   = as.numeric(cov_res[names(TRUE_VALS)]),
      VarOnly_Est = as.numeric(var_res[names(TRUE_VALS)])
    ) %>%
      mutate(
        Bias_Covar   = (Covar_Est   - True_Val) / True_Val * 100,
        Bias_VarOnly = (VarOnly_Est - True_Val) / True_Val * 100
      )

  }, error = function(e) NULL)
}


## ----- Run parallel grid -----

task_grid <- expand.grid(n_sub = n_subjects, rep = 1:n_reps)

start_time <- Sys.time()
cat("Simulation started at:", format(start_time, "%H:%M:%S"), "\n")

with_progress({
  p <- progressor(nrow(task_grid))
  results_list <- future_lapply(1:nrow(task_grid), function(idx) {
    res <- run_single_iteration(task_grid$n_sub[idx], task_grid$rep[idx])
    p()
    res
  }, future.seed = TRUE)
})

end_time <- Sys.time()
cat("Total Runtime:",
    round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

full_iteration_data <- dplyr::bind_rows(results_list)


## ----- Process and save results -----

summary_data <- full_iteration_data %>%
  group_by(N_Subject, Metric) %>%
  summarise(across(starts_with("Bias"),
                   list(Mean = mean, Var = var), na.rm = TRUE),
            .groups = 'drop')

output_file <- "Subjects_Results.xlsx"
write_xlsx(list(
  "All_Iterations" = full_iteration_data,
  "Summary_Stats"  = summary_data
), output_file)

cat("\nSimulation complete! Results saved to:", output_file)


## ----- Visualization -----

pk_levels    <- c("cl", "v1", "v2", "q")
pk_labels    <- c("CL (L/h)", "V1 (L)", "V2 (L)", "Q (L/h)")
aux_levels   <- c("auc", "cmax", "residual_error")
aux_labels   <- c("AUC", "Cmax", "Residual Error")
omega_levels <- c("omega_cl", "omega_v1", "omega_v2", "omega_q")
omega_labels <- c("omega(CL)", "omega(V1)", "omega(V2)", "omega(Q)")

method_colors <- c("Full Covariance" = "#2C7BB6", "Variance Only" = "#D7191C")
method_colors_shaded <- c("Full Covariance" = "#0072B2", "Variance Only" = "#D7191C")

# Helper: pivot to long bias format for a given metric subset
make_plot_data <- function(data, levels, labels) {
  data %>%
    dplyr::filter(Metric %in% levels) %>%
    dplyr::select(N_Subject, Metric, Bias_Covar, Bias_VarOnly) %>%
    tidyr::pivot_longer(
      cols = c(Bias_Covar, Bias_VarOnly),
      names_to = "Method", values_to = "Relative_Bias",
      names_prefix = "Bias_"
    ) %>%
    mutate(
      N_Label = factor(paste0("N=", N_Subject),
                       levels = paste0("N=", sort(unique(N_Subject)))),
      Method  = ifelse(Method == "Covar", "Full Covariance", "Variance Only"),
      Metric  = factor(Metric, levels = levels, labels = labels)
    )
}

# Helper: summarise to mean/SD per N, metric, method
make_summary_data <- function(data, levels, labels) {
  data %>%
    dplyr::filter(Metric %in% levels) %>%
    tidyr::pivot_longer(
      cols = c(Bias_Covar, Bias_VarOnly),
      names_to = "Method", values_to = "Relative_Bias",
      names_prefix = "Bias_"
    ) %>%
    group_by(N_Subject, Metric, Method) %>%
    summarise(
      Mean_Bias = mean(Relative_Bias, na.rm = TRUE),
      SD_Bias   = sd(Relative_Bias,   na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      Method = ifelse(Method == "Covar", "Full Covariance", "Variance Only"),
      Metric = factor(Metric, levels = levels, labels = labels)
    )
}

# Helper: boxplot
bias_boxplot <- function(data, title) {
  ggplot(data, aes(x = N_Label, y = Relative_Bias, fill = Method)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5, width = 0.7,
                 position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
    facet_wrap(~Metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = method_colors) +
    labs(title = title, x = "Sample Size", y = "Relative Bias (%)", fill = "Method",
         caption = "Red dashed line = zero bias") +
    theme_bw(base_size = 12) +
    theme(legend.position   = "right",
          strip.background  = element_rect(fill = "grey95"),
          strip.text        = element_text(face = "bold"),
          panel.grid.minor  = element_blank(),
          plot.title        = element_text(face = "bold", size = 14),
          axis.title        = element_text(face = "bold", size = 10))
}

# Helper: shaded convergence plot
bias_shaded <- function(data, title) {
  ggplot(data, aes(x = N_Subject, y = Mean_Bias, color = Method, fill = Method)) +
    geom_ribbon(aes(ymin = Mean_Bias - SD_Bias, ymax = Mean_Bias + SD_Bias),
                alpha = 0.2, linetype = 0) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    facet_wrap(~Metric, scales = "free_y", ncol = 2) +
    scale_x_continuous(breaks = n_subjects) +
    scale_color_manual(values = method_colors_shaded) +
    scale_fill_manual(values  = method_colors_shaded) +
    labs(title    = title,
         subtitle = "Solid lines = Mean Relative Bias; Shaded areas = \u00b11 SD",
         x = "Sample Size (N)", y = "Relative Bias (%)",
         color = "Estimation Method", fill = "Estimation Method") +
    theme_minimal(base_size = 12) +
    theme(legend.position  = "right",
          strip.text       = element_text(face = "bold", size = 11),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold"))
}


## --- PK parameters + auxiliary ---

p_pk <- bias_boxplot(
  make_plot_data(full_iteration_data, c(pk_levels, aux_levels), c(pk_labels, aux_labels)),
  "Parameter Estimation Bias: PK Parameters"
)

p_pk_shaded <- bias_shaded(
  make_summary_data(full_iteration_data, c(pk_levels, aux_levels), c(pk_labels, aux_labels)),
  "Mean Bias Convergence: PK Parameters"
)

## --- Omega (IIV variances) ---

p_omega <- bias_boxplot(
  make_plot_data(full_iteration_data, omega_levels, omega_labels),
  "Parameter Estimation Bias: IIV (Omega variances)"
)

p_omega_shaded <- bias_shaded(
  make_summary_data(full_iteration_data, omega_levels, omega_labels),
  "Mean Bias Convergence: IIV (Omega variances)"
)

print(p_pk)
print(p_pk_shaded)
print(p_omega)
print(p_omega_shaded)
