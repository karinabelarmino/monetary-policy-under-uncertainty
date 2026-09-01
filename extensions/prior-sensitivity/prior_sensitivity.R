# Monetary Policy Under Uncertainty
# Prior-sensitivity exercise.
# Run the sections in order from the repository root.

# PACKAGES -----------------------------------------------------------------------------

required_packages <- c("BVAR", "coda")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the required packages before running the exercise: ",
    paste(missing_packages, collapse = ", "),
    "."
  )
}

library(BVAR)
library(coda)

if (as.character(utils::packageVersion("BVAR")) != "1.0.5") {
  warning("This exercise was prepared with BVAR 1.0.5.")
}

# USER SETTINGS ------------------------------------------------------------------------

# Keep FALSE for the complete 36-month evaluation used in the release.
# Change to TRUE only to check locally whether every section runs.
quick_mode <- FALSE

# Sign-restricted IRFs are computationally intensive. They are skipped in quick mode.
run_structural_comparison <- !quick_mode

set.seed(1)

if (quick_mode) {
  holdout <- 6
  evaluation_draws <- 1200
  evaluation_burn <- 400
  full_draws <- 4000
  full_burn <- 1000
  structural_draws <- 200
} else {
  holdout <- 36
  evaluation_draws <- 6000
  evaluation_burn <- 2000
  full_draws <- 60000
  full_burn <- 15000
  structural_draws <- 500
}

lags <- 2
output_dir <- "extensions/prior-sensitivity/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# LOADING DATABASE ---------------------------------------------------------------------

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the repository root.")
}

load("data/monetary_policy_data.RData")

model_variables <- c("igi", "ci", "gha", "exp", "selic")

if (!exists("data") || !all(c("date", model_variables) %in% names(data))) {
  stop("The frozen data file does not contain the expected variables.")
}

dates <- as.Date(data$date)
Y <- as.matrix(data[, model_variables])
storage.mode(Y) <- "double"

if (anyNA(dates) || anyNA(Y) || any(!is.finite(Y))) {
  stop("Dates and model variables must be complete and finite.")
}

# COMMON FUNCTIONS ---------------------------------------------------------------------

crps_sample <- function(draws, observation) {
  draws <- sort(as.numeric(draws))
  n_draws <- length(draws)
  first_term <- mean(abs(draws - observation))
  weights <- 2 * seq_len(n_draws) - n_draws - 1
  pairwise_term <- 2 * sum(weights * draws) / n_draws^2
  first_term - 0.5 * pairwise_term
}

make_var_matrices <- function(series, lags) {
  n_obs <- nrow(series)
  X <- matrix(1, nrow = n_obs - lags, ncol = 1 + ncol(series) * lags)

  for (lag_id in seq_len(lags)) {
    columns <- (2 + (lag_id - 1) * ncol(series)):(1 + lag_id * ncol(series))
    X[, columns] <- series[(lags + 1 - lag_id):(n_obs - lag_id), , drop = FALSE]
  }

  list(
    Y = series[(lags + 1):n_obs, , drop = FALSE],
    X = X
  )
}

posterior_stability <- function(beta_draws, n_vars, lags, max_draws = 2000) {
  draw_ids <- unique(round(seq(
    1,
    dim(beta_draws)[1],
    length.out = min(max_draws, dim(beta_draws)[1])
  )))

  stable <- vapply(draw_ids, function(draw_id) {
    top_block <- t(beta_draws[draw_id, -1, , drop = FALSE][1, , ])
    companion <- rbind(
      top_block,
      cbind(diag(n_vars * (lags - 1)), matrix(0, n_vars * (lags - 1), n_vars))
    )
    max(Mod(eigen(companion, only.values = TRUE)$values)) < 1
  }, logical(1))

  mean(stable)
}

forecast_from_draws <- function(beta_draws, sigma_draws, training_data, lags) {
  last_values <- unlist(
    lapply(seq_len(lags), function(lag_id) training_data[nrow(training_data) + 1 - lag_id, ])
  )
  predictor <- c(1, last_values)
  forecasts <- matrix(NA_real_, dim(beta_draws)[1], ncol(training_data))

  for (draw_id in seq_len(dim(beta_draws)[1])) {
    forecast_mean <- as.numeric(predictor %*% beta_draws[draw_id, , ])
    forecast_error <- t(chol(sigma_draws[draw_id, , ])) %*% rnorm(ncol(training_data))
    forecasts[draw_id, ] <- forecast_mean + as.numeric(forecast_error)
  }

  colnames(forecasts) <- colnames(training_data)
  forecasts
}

# DIAGNOSTIC 1: MINNESOTA SPECIFICATIONS ----------------------------------------------

add_soc <- function(Y, lags, par) {
  soc <- if (lags == 1) {
    diag(Y[1, ]) / par
  } else {
    diag(colMeans(Y[seq_len(lags), , drop = FALSE])) / par
  }

  Y_soc <- soc
  X_soc <- cbind(
    rep(0, ncol(Y)),
    matrix(rep(soc, lags), nrow = ncol(Y))
  )

  return(list(Y = Y_soc, X = X_soc))
}

add_sur <- function(Y, lags, par) {
  sur <- if (lags == 1) {
    Y[1, ] / par
  } else {
    colMeans(Y[seq_len(lags), , drop = FALSE]) / par
  }

  Y_sur <- sur
  X_sur <- c(1 / par, rep(sur, lags))

  return(list(Y = Y_sur, X = X_sur))
}

make_minnesota_prior <- function(prior_id) {
  mn <- bv_minnesota(
    lambda = bv_lambda(mode = 0.5, sd = 0.4, min = 0.0001, max = 5),
    alpha = bv_alpha(mode = 1),
    var = 1e07,
    b = if (prior_id == "mn_zero_mean") 0 else 1
  )

  prior_arguments <- list(
    hyper = c("lambda", "alpha"),
    mn = mn
  )

  if (prior_id %in% c("mn_soc", "published_baseline")) {
    prior_arguments$soc <- bv_dummy(
      mode = 1,
      sd = 1,
      min = 0.0001,
      max = 50,
      fun = add_soc
    )
  }

  if (prior_id %in% c("mn_sur", "published_baseline")) {
    prior_arguments$sur <- bv_dummy(
      mode = 1,
      sd = 1,
      min = 0.0001,
      max = 50,
      fun = add_sur
    )
  }

  do.call(bv_priors, prior_arguments)
}

mh <- bv_metropolis(adjust_acc = TRUE, adjust_burn = 0.50)

# DIAGNOSTIC 2: NON-MINNESOTA PRIORS --------------------------------------------------

# Both alternatives below retain the same Gaussian VAR likelihood and constant
# covariance matrix. The variables are standardized inside each estimation window so
# that the non-Minnesota prior scales have the same interpretation across equations.

draw_non_minnesota <- function(series, prior_id, n_draw, n_burn, seed) {
  set.seed(seed)

  series_mean <- colMeans(series)
  series_sd <- apply(series, 2, sd)
  standardized <- sweep(sweep(series, 2, series_mean, "-"), 2, series_sd, "/")
  matrices <- make_var_matrices(standardized, lags)
  y <- matrices$Y
  x <- matrices$X

  n_obs <- nrow(y)
  n_vars <- ncol(y)
  n_coefficients <- ncol(x)
  n_parameters <- n_coefficients * n_vars
  n_save <- n_draw - n_burn

  xtx <- crossprod(x)
  xtx_inverse <- solve(xtx)
  beta_ols <- solve(xtx, crossprod(x, y))
  residuals_ols <- y - x %*% beta_ols
  sigma <- crossprod(residuals_ols) / n_obs
  beta <- as.vector(beta_ols)

  prior_mean <- rep(0, n_parameters)
  intercepts <- 1 + (0:(n_vars - 1)) * n_coefficients

  if (prior_id == "normal_wishart") {
    prior_variance <- rep(10^2, n_parameters)
    prior_variance[intercepts] <- 100^2
  }

  if (prior_id == "ssvs") {
    equation_variances <- diag(sigma)
    ols_se <- unlist(lapply(seq_len(n_vars), function(equation_id) {
      sqrt(diag(xtx_inverse) * equation_variances[equation_id])
    }))
    ols_se <- pmax(ols_se, 1e-04)
    tau0 <- 0.1 * ols_se
    tau1 <- 10 * ols_se
    tau1[intercepts] <- 100
    inclusion <- rep(1, n_parameters)
    prior_variance <- tau1^2
  }

  sigma_df_prior <- n_vars + 2
  sigma_scale_prior <- diag(1, n_vars)

  beta_draws <- array(NA_real_, dim = c(n_save, n_coefficients, n_vars))
  sigma_draws <- array(NA_real_, dim = c(n_save, n_vars, n_vars))
  inclusion_sum <- rep(0, n_parameters)

  for (draw_id in seq_len(n_draw)) {
    beta_matrix <- matrix(beta, nrow = n_coefficients, ncol = n_vars)
    residuals <- y - x %*% beta_matrix
    sigma_scale_post <- sigma_scale_prior + crossprod(residuals)
    sigma_precision <- rWishart(
      1,
      sigma_df_prior + n_obs,
      solve(sigma_scale_post)
    )[, , 1]
    sigma <- solve(sigma_precision)

    prior_precision <- 1 / prior_variance
    beta_precision <- kronecker(sigma_precision, xtx) + diag(prior_precision)
    beta_rhs <- as.vector(crossprod(x, y) %*% sigma_precision) +
      prior_precision * prior_mean
    beta_mean <- solve(beta_precision, beta_rhs)
    beta_chol <- chol(beta_precision)
    beta <- beta_mean + backsolve(beta_chol, rnorm(n_parameters))

    if (prior_id == "ssvs") {
      log_spike <- dnorm(beta, 0, tau0, log = TRUE)
      log_slab <- dnorm(beta, 0, tau1, log = TRUE)
      inclusion_probability <- 1 / (1 + exp(log_spike - log_slab))
      inclusion <- rbinom(n_parameters, 1, inclusion_probability)
      inclusion[intercepts] <- 1
      prior_variance <- ifelse(inclusion == 1, tau1^2, tau0^2)
    }

    if (draw_id > n_burn) {
      save_id <- draw_id - n_burn
      beta_standardized <- matrix(beta, nrow = n_coefficients, ncol = n_vars)
      beta_original <- matrix(NA_real_, nrow = n_coefficients, ncol = n_vars)
      lag_matrices <- vector("list", lags)
      scale_matrix <- diag(series_sd)
      inverse_scale <- diag(1 / series_sd)

      for (lag_id in seq_len(lags)) {
        rows <- (2 + (lag_id - 1) * n_vars):(1 + lag_id * n_vars)
        lag_matrices[[lag_id]] <- scale_matrix %*%
          t(beta_standardized[rows, , drop = FALSE]) %*%
          inverse_scale
        beta_original[rows, ] <- t(lag_matrices[[lag_id]])
      }

      beta_original[1, ] <- series_mean +
        series_sd * beta_standardized[1, ] -
        Reduce(`+`, lapply(lag_matrices, function(a) a %*% series_mean))

      beta_draws[save_id, , ] <- beta_original
      sigma_draws[save_id, , ] <- scale_matrix %*% sigma %*% scale_matrix

      if (prior_id == "ssvs") {
        inclusion_sum <- inclusion_sum + inclusion
      }
    }
  }

  list(
    beta = beta_draws,
    sigma = sigma_draws,
    inclusion_probability = if (prior_id == "ssvs") {
      inclusion_sum / n_save
    } else {
      rep(NA_real_, n_parameters)
    }
  )
}

as_bvar_object <- function(draws, series, n_draw, n_burn) {
  matrices <- make_var_matrices(series, lags)
  explanatory_names <- c(
    "constant",
    unlist(lapply(seq_len(lags), function(lag_id) {
      paste0(model_variables, "-lag", lag_id)
    }))
  )

  structure(
    list(
      beta = draws$beta,
      sigma = draws$sigma,
      hyper = matrix(numeric(0), nrow = dim(draws$beta)[1], ncol = 0),
      ml = rep(NA_real_, dim(draws$beta)[1]),
      call = match.call(),
      variables = model_variables,
      explanatories = explanatory_names,
      meta = list(
        accepted = NA_real_,
        Y = matrices$Y,
        X = matrices$X,
        N = nrow(matrices$Y),
        K = ncol(matrices$X),
        M = ncol(series),
        lags = lags,
        n_draw = n_draw,
        n_burn = n_burn,
        n_save = dim(draws$beta)[1],
        n_thin = 1
      )
    ),
    class = "bvar"
  )
}

# MODEL CATALOG ------------------------------------------------------------------------

catalog <- data.frame(
  prior_id = c(
    "published_baseline",
    "mn_random_walk",
    "mn_zero_mean",
    "mn_soc",
    "mn_sur",
    "normal_wishart",
    "ssvs"
  ),
  diagnostic = c(
    "Reference",
    rep("1 - Minnesota specification", 4),
    rep("2 - Alternative prior families", 2)
  ),
  label = c(
    "Published baseline: MN + SOC + SUR",
    "Minnesota: random-walk mean",
    "Minnesota: zero mean",
    "Minnesota + SOC",
    "Minnesota + SUR",
    "Independent Normal-Wishart",
    "SSVS spike-and-slab"
  ),
  engine = c(rep("BVAR", 5), rep("custom_gibbs", 2)),
  stringsAsFactors = FALSE
)

# EXPANDING ONE-STEP FORECASTS ---------------------------------------------------------

holdout_start <- nrow(Y) - holdout + 1
origins <- (holdout_start - 1):(nrow(Y) - 1)
reference_scale <- apply(Y[seq_len(holdout_start - 1), , drop = FALSE], 2, sd)
forecast_results <- vector("list", nrow(catalog) * length(origins))
result_id <- 1

message(
  "Estimating ", nrow(catalog), " specifications over ", holdout,
  " expanding one-step forecasts."
)

for (model_id in seq_len(nrow(catalog))) {
  prior_id <- catalog$prior_id[model_id]
  message("Evaluating: ", catalog$label[model_id])

  for (origin_id in seq_along(origins)) {
    origin <- origins[origin_id]
    training_data <- Y[seq_len(origin), , drop = FALSE]
    observation <- Y[origin + 1, ]

    if (catalog$engine[model_id] == "BVAR") {
      set.seed(100000 + 1000 * model_id + origin_id)
      fit <- bvar(
        training_data,
        lags = lags,
        n_draw = evaluation_draws,
        n_burn = evaluation_burn,
        n_thin = 1,
        priors = make_minnesota_prior(prior_id),
        mh = mh,
        verbose = FALSE
      )

      set.seed(200000 + 1000 * model_id + origin_id)
      forecast <- predict(fit, horizon = 1, n_thin = 2)
      forecast_draws <- forecast$fcast[, 1, , drop = FALSE]
      dim(forecast_draws) <- c(dim(forecast$fcast)[1], length(model_variables))
    } else {
      draws <- draw_non_minnesota(
        training_data,
        prior_id,
        evaluation_draws,
        evaluation_burn,
        100000 + 1000 * model_id + origin_id
      )
      set.seed(200000 + 1000 * model_id + origin_id)
      forecast_draws <- forecast_from_draws(
        draws$beta,
        draws$sigma,
        training_data,
        lags
      )
    }

    forecast_mean <- colMeans(forecast_draws)
    lower_68 <- apply(forecast_draws, 2, quantile, probs = 0.16)
    upper_68 <- apply(forecast_draws, 2, quantile, probs = 0.84)
    lower_90 <- apply(forecast_draws, 2, quantile, probs = 0.05)
    upper_90 <- apply(forecast_draws, 2, quantile, probs = 0.95)

    forecast_results[[result_id]] <- data.frame(
      prior_id = prior_id,
      date = dates[origin + 1],
      variable = model_variables,
      observation = as.numeric(observation),
      forecast_mean = forecast_mean,
      error = forecast_mean - as.numeric(observation),
      absolute_error = abs(forecast_mean - as.numeric(observation)),
      squared_error = (forecast_mean - as.numeric(observation))^2,
      crps = vapply(seq_along(model_variables), function(variable_id) {
        crps_sample(forecast_draws[, variable_id], observation[variable_id])
      }, numeric(1)),
      coverage_68 = observation >= lower_68 & observation <= upper_68,
      coverage_90 = observation >= lower_90 & observation <= upper_90,
      scale = as.numeric(reference_scale[model_variables]),
      stringsAsFactors = FALSE
    )

    result_id <- result_id + 1
  }
}

forecast_scores <- do.call(rbind, forecast_results)
rownames(forecast_scores) <- NULL

# FORECAST METRICS ---------------------------------------------------------------------

prior_comparison <- do.call(rbind, lapply(split(forecast_scores, forecast_scores$prior_id), function(x) {
  data.frame(
    prior_id = x$prior_id[1],
    standardized_rmse = sqrt(mean(x$squared_error / x$scale^2)),
    standardized_mae = mean(x$absolute_error / x$scale),
    standardized_crps = mean(x$crps / x$scale),
    coverage_68 = mean(x$coverage_68),
    coverage_90 = mean(x$coverage_90),
    stringsAsFactors = FALSE
  )
}))

prior_comparison <- merge(catalog, prior_comparison, by = "prior_id", sort = FALSE)
prior_comparison$rank_rmse <- rank(prior_comparison$standardized_rmse, ties.method = "min")
prior_comparison$rank_mae <- rank(prior_comparison$standardized_mae, ties.method = "min")
prior_comparison$rank_crps <- rank(prior_comparison$standardized_crps, ties.method = "min")
prior_comparison$mean_rank <- rowMeans(prior_comparison[, c("rank_rmse", "rank_mae", "rank_crps")])
prior_comparison <- prior_comparison[order(prior_comparison$mean_rank, prior_comparison$standardized_crps), ]
prior_comparison$overall_rank <- seq_len(nrow(prior_comparison))
rownames(prior_comparison) <- NULL

diagnostic_1 <- prior_comparison[
  prior_comparison$prior_id %in% c(
    "published_baseline", "mn_random_walk", "mn_zero_mean", "mn_soc", "mn_sur"
  ),
]
diagnostic_2 <- prior_comparison[
  prior_comparison$prior_id %in% c("published_baseline", "normal_wishart", "ssvs"),
]

rank_diagnostic <- function(results, diagnostic_name) {
  results$rank_rmse_diagnostic <- rank(results$standardized_rmse, ties.method = "min")
  results$rank_mae_diagnostic <- rank(results$standardized_mae, ties.method = "min")
  results$rank_crps_diagnostic <- rank(results$standardized_crps, ties.method = "min")
  results$mean_rank_diagnostic <- rowMeans(results[, c(
    "rank_rmse_diagnostic", "rank_mae_diagnostic", "rank_crps_diagnostic"
  )])
  results <- results[order(results$mean_rank_diagnostic, results$standardized_crps), ]
  results$diagnostic_rank <- seq_len(nrow(results))
  results$diagnostic <- diagnostic_name
  results
}

diagnostic_comparison <- rbind(
  rank_diagnostic(diagnostic_1, "1 - Minnesota specification"),
  rank_diagnostic(diagnostic_2, "2 - Alternative prior families")
)

best_minnesota_id <- diagnostic_comparison$prior_id[
  diagnostic_comparison$diagnostic == "1 - Minnesota specification" &
    diagnostic_comparison$diagnostic_rank == 1
]
best_non_minnesota_id <- prior_comparison$prior_id[
  prior_comparison$prior_id %in% c("normal_wishart", "ssvs")
][which.min(prior_comparison$mean_rank[
  prior_comparison$prior_id %in% c("normal_wishart", "ssvs")
])]

utils::write.csv(
  forecast_scores,
  file.path(output_dir, "forecast_scores.csv"),
  row.names = FALSE
)
utils::write.csv(
  prior_comparison,
  file.path(output_dir, "prior_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  diagnostic_comparison,
  file.path(output_dir, "diagnostic_comparison.csv"),
  row.names = FALSE
)

# FULL-SAMPLE DIAGNOSTICS ---------------------------------------------------------------

message("Estimating the seven specifications on the full sample.")

full_fits <- setNames(vector("list", nrow(catalog)), catalog$prior_id)
posterior_diagnostics <- vector("list", nrow(catalog))

for (model_id in seq_len(nrow(catalog))) {
  prior_id <- catalog$prior_id[model_id]
  message("Full sample: ", catalog$label[model_id])

  if (catalog$engine[model_id] == "BVAR") {
    set.seed(300000 + model_id)
    fit <- bvar(
      Y,
      lags = lags,
      n_draw = full_draws,
      n_burn = full_burn,
      n_thin = 1,
      priors = make_minnesota_prior(prior_id),
      mh = mh,
      verbose = FALSE
    )
    monitoring_draws <- as.mcmc(fit, vars = c("lambda", "alpha"))
    inclusion_probability <- NA_real_
  } else {
    draws <- draw_non_minnesota(
      Y,
      prior_id,
      full_draws,
      full_burn,
      300000 + model_id
    )
    fit <- as_bvar_object(draws, Y, full_draws, full_burn)
    monitoring_draws <- mcmc(cbind(
      own_first_lag = fit$beta[, 2, 1],
      residual_variance = fit$sigma[, 1, 1]
    ))
    inclusion_probability <- if (prior_id == "ssvs") {
      coefficient_count <- 1 + ncol(Y) * lags
      intercept_ids <- 1 + (0:(ncol(Y) - 1)) * coefficient_count
      mean(draws$inclusion_probability[-intercept_ids], na.rm = TRUE)
    } else {
      NA_real_
    }
  }

  geweke <- geweke.diag(monitoring_draws)$z
  effective_size <- effectiveSize(monitoring_draws)
  posterior_diagnostics[[model_id]] <- data.frame(
    prior_id = prior_id,
    minimum_effective_sample_size = min(effective_size),
    maximum_absolute_geweke_z = max(abs(geweke)),
    stable_draw_share = posterior_stability(fit$beta, ncol(Y), lags),
    mean_inclusion_probability = inclusion_probability,
    stringsAsFactors = FALSE
  )
  full_fits[[prior_id]] <- fit
}

posterior_diagnostics <- do.call(rbind, posterior_diagnostics)
utils::write.csv(
  posterior_diagnostics,
  file.path(output_dir, "posterior_diagnostics.csv"),
  row.names = FALSE
)

flagged_diagnostics <- posterior_diagnostics$prior_id[
  posterior_diagnostics$maximum_absolute_geweke_z > 1.96
]

if (length(flagged_diagnostics) > 0) {
  warning(
    "Inspect the Geweke diagnostic before release for: ",
    paste(flagged_diagnostics, collapse = ", "),
    "."
  )
}

# FORECAST FIGURE ----------------------------------------------------------------------

catalog_lookup <- setNames(catalog$label, catalog$prior_id)

grDevices::png(
  file.path(output_dir, "prior-performance.png"),
  width = 2200,
  height = 1050,
  res = 180
)
old_par <- par(mfrow = c(1, 2), mar = c(5, 10, 4, 3), las = 1)

short_labels <- c(
  published_baseline = "Published baseline",
  mn_random_walk = "MN: random-walk",
  mn_zero_mean = "MN: zero mean",
  mn_soc = "MN + SOC",
  mn_sur = "MN + SUR",
  normal_wishart = "Normal-Wishart",
  ssvs = "SSVS spike-and-slab"
)

short_titles <- c(
  "1 - Minnesota specification" = "Diagnostic 1: Minnesota",
  "2 - Alternative prior families" = "Diagnostic 2: other prior families"
)

for (diagnostic_name in c(
  "1 - Minnesota specification",
  "2 - Alternative prior families"
)) {
  plot_data <- diagnostic_comparison[diagnostic_comparison$diagnostic == diagnostic_name, ]
  plot_data <- plot_data[order(plot_data$standardized_crps, decreasing = TRUE), ]
  colors <- ifelse(plot_data$diagnostic_rank == 1, "#2f8f5b", "#687078")
  midpoints <- barplot(
    plot_data$standardized_crps,
    names.arg = short_labels[plot_data$prior_id],
    horiz = TRUE,
    col = colors,
    border = NA,
    xlab = "Standardized CRPS (lower is better)",
    main = short_titles[diagnostic_name],
    cex.names = 0.75,
    xlim = c(0, max(plot_data$standardized_crps) * 1.25)
  )
  text(
    plot_data$standardized_crps,
    midpoints,
    labels = sprintf(" %.3f", plot_data$standardized_crps),
    pos = 4,
    cex = 0.8
  )
}

par(old_par)
dev.off()

# STRUCTURAL SENSITIVITY ---------------------------------------------------------------

if (run_structural_comparison) {
  message("Computing sign-restricted IRFs for the baseline and both diagnostic leaders.")

  restrictions <- matrix(
    c(
      1, NA, -1,  1,  0,
      NA,  1,  1, -1,  0,
      NA, NA, NA, NA, NA,
      NA, NA, NA, NA, NA,
      NA, NA, -1, -1,  1
    ),
    ncol = 5,
    nrow = 5
  )

  irf_options <- bv_irf(
    sign_restr = restrictions,
    horizon = 24,
    fevd = FALSE,
    identification = TRUE,
    sign_lim = 60000
  )

  structural_ids <- unique(c(
    "published_baseline",
    best_minnesota_id,
    best_non_minnesota_id
  ))
  structural_irfs <- setNames(vector("list", length(structural_ids)), structural_ids)

  for (structural_id in seq_along(structural_ids)) {
    prior_id <- structural_ids[structural_id]
    retained_draws <- full_fits[[prior_id]]$meta$n_save
    irf_thin <- max(1, floor(retained_draws / structural_draws))
    set.seed(400000 + structural_id)
    structural_irfs[[prior_id]] <- BVAR::irf(
      full_fits[[prior_id]],
      irf_options,
      conf_bands = c(0.16),
      n_thin = irf_thin,
      verbose = FALSE
    )
  }

  response_paths <- data.frame(
    impulse = c("igi", "igi", "ci", "selic"),
    response = c("exp", "gha", "exp", "exp"),
    title = c(
      "Uncertainty shock -> inflation expectations",
      "Uncertainty shock -> output gap",
      "Credibility shock -> inflation expectations",
      "Selic shock -> inflation expectations"
    ),
    stringsAsFactors = FALSE
  )

  extract_irf <- function(irf_object, impulse, response) {
    response_id <- match(response, irf_object$variables)
    impulse_id <- match(impulse, irf_object$variables)
    values <- irf_object$quants[, response_id, , impulse_id, drop = FALSE]
    dim(values) <- c(dim(irf_object$quants)[1], dim(irf_object$quants)[3])
    values
  }

  line_colors <- c("#4f565c", "#2f8f5b", "#b25f36")
  names(line_colors) <- structural_ids

  grDevices::png(
    file.path(output_dir, "irf-sensitivity.png"),
    width = 1800,
    height = 1300,
    res = 170
  )
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 4, 1), oma = c(1, 1, 3, 1))
  horizons <- 0:23

  for (path_id in seq_len(nrow(response_paths))) {
    values <- lapply(structural_irfs, function(irf_result) {
      extract_irf(
        irf_result,
        response_paths$impulse[path_id],
        response_paths$response[path_id]
      )
    })
    y_range <- range(unlist(lapply(values, function(x) x[c(1, 3), ])), 0, finite = TRUE)

    plot(
      horizons,
      values[[1]][2, ],
      type = "n",
      ylim = y_range,
      xlab = "Months",
      ylab = "Response",
      main = response_paths$title[path_id]
    )
    abline(h = 0, col = "#a8adb2", lty = 2)

    for (prior_id in structural_ids) {
      polygon(
        c(horizons, rev(horizons)),
        c(values[[prior_id]][1, ], rev(values[[prior_id]][3, ])),
        border = NA,
        col = grDevices::adjustcolor(line_colors[prior_id], alpha.f = 0.08)
      )
      lines(horizons, values[[prior_id]][2, ], col = line_colors[prior_id], lwd = 2)
    }

    if (path_id == 1) {
      legend(
        "topright",
        legend = catalog_lookup[structural_ids],
        col = line_colors[structural_ids],
        lwd = 3,
        cex = 0.7,
        bty = "n"
      )
    }
  }

  mtext("Structural sensitivity under the published restrictions", outer = TRUE, cex = 1.25, font = 2)
  par(old_par)
  dev.off()
}

# SESSION INFORMATION ------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "session-info.txt")
)

message("Prior-sensitivity exercise completed.")
print(prior_comparison[, c(
  "overall_rank",
  "label",
  "standardized_rmse",
  "standardized_mae",
  "standardized_crps"
)], row.names = FALSE)
