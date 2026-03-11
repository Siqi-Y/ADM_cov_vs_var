# Meropenem in CRRT patients
# two-compartment model with linear elimination (iv)

## ----- Loading library ----
library(pacman)
pacman::p_load(admr, rxode2, tidyverse, MASS, data.table, pracma, 
               writexl, scales, future, future.apply, progressr, parallel)

# Enable progress handlers globally
handlers(global = TRUE)
handlers("progress")

# Sets up parallel processing 
plan(multisession, workers = parallel::detectCores() - 1)

## ----- Define Simulation Conditions ----

n_subjects = c(20, 50, 80, 100, 150, 200) # Different sample size scenarios
n_timepoints = 8
n_reps <- 50

get_timepoints <- function(n) {
  switch(as.character(n),
         "4" = c(0.1, 1, 6, 12), 
         "5" = c(0.1, 0.75, 2, 6, 12),
         "6" = c(0.1, 0.5, 2, 4, 8, 12),
         "8" = c(0.1, 0.5, 1, 2, 4, 6, 8, 12),
         "12" = c(0.08, 0.25, 0.5, 1, 1.5, 2, 3, 5, 7, 9, 10, 12),
         stop("Unsupported number of timepoints"))
}

times <- get_timepoints(n_timepoints)

## ----- Define Global True Parameters  ----

TRUE_PARAMS <- list(
  theta = c(cl = 3.03, v1 = 14.2, v2 = 14.6, q = 15.9),
  iiv_sd = c(cl = 0.52, v1 = 0.33, v2 = 0.33, q = 0.26),
  sigma_prop = 0.153,
  dose = 1000
)

TRUE_OMEGA <- diag(TRUE_PARAMS$iiv_sd^2)

# Compile Model 
rxModel <- rxode2({
  ke = cl / v1
  k12 = q / v1
  k21 = q / v2
  d/dt(central) = -ke*central - k12*central + k21*peripheral
  d/dt(peripheral) = k12*central - k21*peripheral
  cp = central / v1
})


# Pre-create PK evaluation event table
pk_ev <- et() %>%
  et(amt = TRUE_PARAMS$dose, time = 0) %>%
  add.sampling(seq(0, 12, 0.01))

# Initial calculation for TRUE_VALS
TRUE_PK <- rxSolve(rxModel, TRUE_PARAMS$theta, pk_ev)
TRUE_VALS <- c(TRUE_PARAMS$theta,
               auc = trapz(TRUE_PK$time, TRUE_PK$cp),
               cmax = max(TRUE_PK$cp))


## ----- Main Simulation (Parallel Execution) -----

# Extract results funciton
extract_results <- function(fit, true_vals) {
  tryCatch({
    params <- c("cl", "v1", "v2", "q")
    
    p_df <- fit$param_df[fit$param_df$Parameter %in% params, ]
    idx <- match(params, p_df$Parameter)
    est_str <- p_df$`Back-transformed(95%CI)`[idx]
    est <- as.numeric(sub(" \\(.*", "", est_str))
    names(est) <- params
    
    local_model <- rxode2({
      ke = cl / v1
      k12 = q / v1
      k21 = q / v2
      d/dt(central) = -ke*central - k12*central + k21*peripheral
      d/dt(peripheral) = k12*central - k21*peripheral
      cp = central / v1
    })
    
    local_ev <- et() %>%
      et(amt = TRUE_PARAMS$dose, time = 0) %>%
      add.sampling(seq(0, 12, 0.01))
    
    # Solve model and calculate PK metrics
    out <- rxSolve(local_model, params = est, events = local_ev, 
                   inits = c(central=0, peripheral=0))
    auc <- trapz(out$time, out$cp)
    cmax <- max(out$cp)
    
    res_vec <- c(est, auc = auc, cmax = cmax)
    res_vec[names(true_vals)]
    
  }, error = function(e) NULL)
}

run_single_iteration <- function(n_sub, i){
  tryCatch({
    set.seed(n_sub * 1000 + i)
    
    ## ---- Generate data ----
    
    # Create individual parameters
    mv <- MASS::mvrnorm(n_sub, rep(0, nrow(TRUE_OMEGA)), TRUE_OMEGA)
    params_all <- data.table(
      ID = 1:n_sub,
      cl = TRUE_PARAMS$theta['cl'] * exp(mv[, 1]),
      v1 = TRUE_PARAMS$theta['v1'] * exp(mv[, 2]),
      v2 = TRUE_PARAMS$theta['v2'] * exp(mv[, 3]),
      q  = TRUE_PARAMS$theta['q']  * exp(mv[, 4])
    )
    
    # Create event table
    ev <- et() %>%
      et(amt = TRUE_PARAMS$dose, time = 0) %>%
      et(time = times) %>%
      et(ID = 1:n_sub)
    
    # Solve the model
    sim <- rxSolve(rxModel, params_all, ev, cores = 1) %>%
      as.data.frame() %>%
      dplyr::filter(time > 0) %>%  
      # dplyr::mutate(DV = cp * (1 + rnorm(dplyr::n(), 0, TRUE_PARAMS$sigma_prop)) + 
      #                 rnorm(dplyr::n(), 0, TRUE_PARAMS$sigma_add)) %>%
      dplyr::mutate(DV = pmax(cp, 0.01))  
    
    ind_wide <- sim %>%
      dplyr::select(id, time, DV) %>%
      tidyr::pivot_wider(names_from = time, values_from = DV) %>%
      dplyr::select(-id) %>%
      as.matrix()
    
    agg_full <- admr::meancov(ind_wide)
    diag(agg_full$V) <- diag(agg_full$V) + (TRUE_PARAMS$sigma_prop^2) * (agg_full$E^2)
    
    agg_var  <- agg_full
    agg_var$V <- diag(diag(agg_var$V)) # Variance-only matrix
    
    
    ## ----- Fitting model -----
    
    # Define predder inside worker 
    predder <- function(time, theta_i) {
      n_ind <- nrow(theta_i)
      if (is.null(n_ind) || n_ind == 0) {
        theta_i <- as.data.frame(t(theta_i))
        n_ind <- 1
      }
      
      # Create event table
      ev <- eventTable(amount.units="mg", time.units="hours")
      ev$add.dosing(dose = TRUE_PARAMS$dose, start.time = 0)
      ev$add.sampling(time)
      
      # Solve model
      out <- rxSolve(rxModel, theta_i, ev, cores = 1)
      
      # Format output
      cp_matrix <- matrix(out$cp, nrow = n_ind, ncol = length(time), 
                          byrow = TRUE)
      
      return(cp_matrix)
    }
    
    ## Set Options
    base_opts <- list(
      time = times,
      p = list(beta = TRUE_PARAMS$theta, Omega = TRUE_OMEGA,
               Sigma_prop = TRUE_PARAMS$sigma_prop^2),
      nsim = 2000,
      n = n_sub,
      fo_appr = FALSE,
      omega_expansion = 1,
      f = predder
    )
    
    opts_covar <- do.call(admr::genopts, c(base_opts, list(no_cov = FALSE)))
    opts_var <- do.call(admr::genopts, c(base_opts, list(no_cov = TRUE)))
    
    fit_covar <- tryCatch(
      admr::fitIRMC(
        opts = opts_covar,
        obs = agg_full,
        chains = 2,  # Number of chains
        maxiter = 800,  # Maximum iterations
        use_grad = T
      ), error = function(e) NULL)
    
    fit_var <- tryCatch(
      admr::fitIRMC(
        opts = opts_var,
        obs = agg_var,
        chains = 2,  # Number of chains
        maxiter = 800,  # Maximum iterations
        use_grad = T
      ),error = function(e) NULL)
    
    ##----- Result Handling -----
    
    if (is.null(fit_covar) || is.null(fit_var)) return(NULL)
    
    cov_res <- extract_results(fit_covar, TRUE_VALS)
    var_res <- extract_results(fit_var, TRUE_VALS)
    
    if (is.null(cov_res) || is.null(var_res)) return(NULL)
    
    # Output
    data.frame(
      N_Subject = n_sub,
      Rep = i,
      Metric = names(TRUE_VALS),
      True_Val = as.numeric(TRUE_VALS),
      Covar_Est = as.numeric(cov_res[names(TRUE_VALS)]),
      VarOnly_Est = as.numeric(var_res[names(TRUE_VALS)])
    ) %>%
      mutate(
        Bias_Covar = (Covar_Est - True_Val) / True_Val * 100,
        Bias_VarOnly = (VarOnly_Est - True_Val) / True_Val * 100
      )
    
  }, error = function(e) NULL)
}


# Run the task grid in parallel
task_grid <- expand.grid(n_sub = n_subjects, rep = 1:n_reps)

start_time <- Sys.time()
cat("Simulation started at:", format(start_time, "%H:%M:%S"), "\n")

# Tracks progress with a progress bar
with_progress({
  p <- progressor(nrow(task_grid))
  results_list <- future_lapply(1:nrow(task_grid), function(idx) {
    res <- run_single_iteration(task_grid$n_sub[idx], task_grid$rep[idx])
    p()
    res
  }, future.seed = T)
})

# Stop Timing 
end_time <- Sys.time()
cat("Total Runtime:", 
    round(difftime(end_time, start_time, units="mins"), 2), "minutes\n")

full_iteration_data <- dplyr::bind_rows(results_list)


## ----- Process and save results -----

# Calculate summary statistics (Mean & Variance) 
summary_data <- full_iteration_data %>%
  group_by(N_Subject, Metric) %>%
  summarise(across(starts_with("Bias"), 
                   list(Mean = mean, Var = var), na.rm = TRUE),
            .groups = 'drop')

# Save to Excel 
output_file <- "Subjects_Results.xlsx"
write_xlsx(list(
  "All_Iterations" = full_iteration_data,
  "Summary_Stats" = summary_data
), output_file)

cat("\nSimulation complete! Results saved to:", output_file)


## ----- Visualization -----

## Bias Comparision boxplot 

# Prepare the data for plotting
plot_data <- full_iteration_data %>%
  dplyr::select(N_Subject, Metric, Bias_Covar, Bias_VarOnly) %>%
  tidyr::pivot_longer(
    cols = c(Bias_Covar, Bias_VarOnly),
    names_to = "Method",
    values_to = "Relative_Bias",
    names_prefix = "Bias_"
  ) %>%
  mutate(
    N_Label = factor(paste0("N=", N_Subject), 
                     levels = paste0("N=", sort(unique(N_Subject)))),
    Method = ifelse(Method == "Covar", "Full Covariance", "Variance Only"),
    Metric = factor(Metric, 
                    levels = c("cl", "v1", "v2", "q", "auc", "cmax"),
                    labels = c("CL (L/h)", "V1 (L)", "V2 (L)", "Q (L/h)", "AUC", "Cmax"))
  )

#  Generate plot
p <- ggplot(plot_data, aes(x = N_Label, y = Relative_Bias, fill = Method)) +
  # Use a dodged boxplot to show distribution
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, width = 0.7, 
               position = position_dodge(width = 0.8)) +
  
  # Add a reference line at 0 (the 'True Value' target)
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.6) +
  
  # Facet by metric to compare parameters separately
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  
  scale_fill_manual(values = c("Full Covariance" = "#2C7BB6", "Variance Only" = "#D7191C")) +
  
  labs(
    title = "Comparison of Parameter Estimation Bias",
    x = "Sample Size (Number of Subjects)",
    y = "Relative Bias (%)",
    fill = "Method",
    caption = "Red dashed line indicates zero bias (true value)."
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 10)
  )


##  Mean Bias  with SD Shading 

summary_plot_data <- full_iteration_data %>%
  # Pivot to long format first
  tidyr::pivot_longer(
    cols = c(Bias_Covar, Bias_VarOnly),
    names_to = "Method",
    values_to = "Relative_Bias",
    names_prefix = "Bias_"
  ) %>%
  # Group by N, Metric, and Method to calculate statistics
  group_by(N_Subject, Metric, Method) %>%
  summarise(
    Mean_Bias = mean(Relative_Bias, na.rm = TRUE),
    SD_Bias = sd(Relative_Bias, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Method = ifelse(Method == "Covar", "Full Covariance", "Variance Only"),
    Metric = factor(Metric, 
                    levels = c("cl", "v1", "v2", "q", "auc", "cmax"),
                    labels = c("CL", "V1", "V2", "Q", "AUC", "Cmax"))
  )

#  Create the Convergence Plot
p_shaded <- ggplot(summary_plot_data, aes(x = N_Subject, y = Mean_Bias, color = Method, fill = Method)) +
  # Add the Shaded Ribbon (Standard Deviation)
  # alpha controls transparency; linetype = 0 removes the ribbon border
  geom_ribbon(aes(ymin = Mean_Bias - SD_Bias, ymax = Mean_Bias + SD_Bias), 
              alpha = 0.2, linetype = 0) +
  
  geom_line(linewidth = 0.8) +   # Add the Mean Line
  geom_point(size = 1) +    # Add Points for each N value
  
  # Horizontal reference line for zero bias
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +    # Facet by Metric
  
  # Professional Scales and Colors
  scale_x_continuous(breaks = n_subjects) + 
  scale_color_manual(values = c("Full Covariance" = "#0072B2", "Variance Only" = "#D7191C")) +
  scale_fill_manual(values = c("Full Covariance" = "#0072B2", "Variance Only" = "#D7191C")) +
  
  labs(
    title = "Mean Bias Convergence with SD Shading",
    subtitle = "Solid lines represent Mean Relative Bias; Shaded areas represent ±1 Standard Deviation",
    x = "Sample Size (N)",
    y = "Relative Bias (%)",
    fill = "Estimation Method",
    color = "Estimation Method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# Display
print(p)
print(p_shaded)

