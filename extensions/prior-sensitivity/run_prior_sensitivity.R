# Monetary Policy Under Uncertainty
# Prior-sensitivity extension for the published five-variable BVAR.
# Run from the repository root with:
# Rscript extensions/prior-sensitivity/run_prior_sensitivity.R

options(stringsAsFactors = FALSE)

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the repository root.")
}

required_packages <- c("BVAR", "coda")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the required packages before running the exercise: ",
    paste(missing_packages, collapse = ", "),
    "."
  )
}

if (as.character(utils::packageVersion("BVAR")) != "1.0.5") {
  stop("This release is pinned to BVAR 1.0.5.")
}

source("extensions/prior-sensitivity/R/prior_helpers.R")

output_dir <- "extensions/prior-sensitivity/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load("data/monetary_policy_data.RData")

if (!exists("data") || !is.data.frame(data)) {
  stop("The frozen RData file must create one data frame named `data`.")
}

model_variables <- c("igi", "ci", "gha", "exp", "selic")
required_columns <- c("date", model_variables)

if (!all(required_columns %in% names(data))) {
  stop("The frozen data file does not contain the expected baseline variables.")
}

dates <- as.Date(data$date)
Y <- as.matrix(data[, model_variables])
storage.mode(Y) <- "double"

if (anyNA(dates) || anyNA(Y) || any(!is.finite(Y))) {
  stop("Dates and model variables must be complete and finite.")
}

if (anyDuplicated(dates) || !identical(order(dates), seq_along(dates))) {
  stop("Dates must be unique and ordered.")
}

mode <- match.arg(
  tolower(Sys.getenv("PRIOR_SENSITIVITY_MODE", unset = "release")),
  c("release", "quick")
)

settings <- if (mode == "release") {
  list(
    holdout = 36L,
    evaluation_draws = 6000L,
    evaluation_burn = 2000L,
    prediction_thin = 4L,
    full_draws = 60000L,
    full_burn = 15000L,
    structural_draws = 1000L
  )
} else {
  list(
    holdout = 12L,
    evaluation_draws = 3000L,
    evaluation_burn = 1000L,
    prediction_thin = 2L,
    full_draws = 12000L,
    full_burn = 4000L,
    structural_draws = 300L
  )
}

settings$holdout <- get_env_integer(
  "PRIOR_SENSITIVITY_HOLDOUT",
  settings$holdout,
  minimum = 6L
)
settings$evaluation_draws <- get_env_integer(
  "PRIOR_SENSITIVITY_EVAL_DRAWS",
  settings$evaluation_draws,
  minimum = 2000L
)
settings$evaluation_burn <- get_env_integer(
  "PRIOR_SENSITIVITY_EVAL_BURN",
  settings$evaluation_burn,
  minimum = 500L
)
settings$full_draws <- get_env_integer(
  "PRIOR_SENSITIVITY_FULL_DRAWS",
  settings$full_draws,
  minimum = 5000L
)
settings$full_burn <- get_env_integer(
  "PRIOR_SENSITIVITY_FULL_BURN",
  settings$full_burn,
  minimum = 1000L
)

if (settings$evaluation_burn >= settings$evaluation_draws) {
  stop("Evaluation burn-in must be smaller than the number of draws.")
}
if (settings$full_burn >= settings$full_draws) {
  stop("Full-sample burn-in must be smaller than the number of draws.")
}
if (settings$holdout >= nrow(Y) - 24L) {
  stop("The holdout is too large for this sample and lag specification.")
}

run_structural <- get_env_flag(
  "PRIOR_SENSITIVITY_STRUCTURAL",
  default = identical(mode, "release")
)

available_cores <- parallel::detectCores(logical = FALSE)
if (is.na(available_cores)) available_cores <- 1L
default_cores <- if (.Platform$OS.type == "windows") 1L else min(5L, max(1L, available_cores - 1L))
evaluation_cores <- get_env_integer(
  "PRIOR_SENSITIVITY_CORES",
  default_cores,
  minimum = 1L
)

catalog <- prior_catalog()
model_ids <- catalog$prior_id
names(model_ids) <- catalog$label

lags <- 2L
holdout_start <- nrow(Y) - settings$holdout + 1L
origins <- (holdout_start - 1L):(nrow(Y) - 1L)
reference_scale <- apply(Y[seq_len(holdout_start - 1L), , drop = FALSE], 2L, stats::sd)

if (any(!is.finite(reference_scale)) || any(reference_scale <= 0)) {
  stop("Reference scales must be finite and strictly positive.")
}

mh <- BVAR::bv_metropolis(adjust_acc = TRUE, adjust_burn = 0.50)

evaluate_prior <- function(prior_id, prior_index) {
  prior <- make_prior(prior_id)

  result <- lapply(seq_along(origins), function(origin_index) {
    origin <- origins[origin_index]
    training_data <- Y[seq_len(origin), , drop = FALSE]
    observation <- Y[origin + 1L, ]

    set.seed(100000L + 1000L * prior_index + origin_index)
    fit <- BVAR::bvar(
      training_data,
      lags = lags,
      n_draw = settings$evaluation_draws,
      n_burn = settings$evaluation_burn,
      n_thin = 1L,
      priors = prior,
      mh = mh,
      verbose = FALSE
    )

    set.seed(200000L + 1000L * prior_index + origin_index)
    forecast <- stats::predict(
      fit,
      horizon = 1L,
      n_thin = settings$prediction_thin
    )
    draws <- forecast$fcast[, 1L, , drop = FALSE]
    dim(draws) <- c(dim(forecast$fcast)[1L], length(model_variables))

    forecast_mean <- colMeans(draws)
    lower_68 <- apply(draws, 2L, stats::quantile, probs = 0.16, names = FALSE)
    upper_68 <- apply(draws, 2L, stats::quantile, probs = 0.84, names = FALSE)
    lower_90 <- apply(draws, 2L, stats::quantile, probs = 0.05, names = FALSE)
    upper_90 <- apply(draws, 2L, stats::quantile, probs = 0.95, names = FALSE)

    data.frame(
      prior_id = prior_id,
      date = dates[origin + 1L],
      variable = model_variables,
      observation = as.numeric(observation),
      forecast_mean = forecast_mean,
      error = forecast_mean - as.numeric(observation),
      absolute_error = abs(forecast_mean - as.numeric(observation)),
      squared_error = (forecast_mean - as.numeric(observation))^2,
      crps = vapply(
        seq_along(model_variables),
        function(variable_id) crps_sample(draws[, variable_id], observation[variable_id]),
        numeric(1)
      ),
      coverage_68 = observation >= lower_68 & observation <= upper_68,
      coverage_90 = observation >= lower_90 & observation <= upper_90,
      scale = as.numeric(reference_scale[model_variables]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, result)
}

message(
  "Evaluating ", length(model_ids), " prior specifications over ",
  settings$holdout, " expanding one-step forecasts."
)

evaluation_jobs <- Map(
  function(prior_id, prior_index) list(prior_id = prior_id, prior_index = prior_index),
  unname(model_ids),
  seq_along(model_ids)
)

evaluation_results <- if (evaluation_cores > 1L && .Platform$OS.type != "windows") {
  parallel::mclapply(
    evaluation_jobs,
    function(job) evaluate_prior(job$prior_id, job$prior_index),
    mc.cores = min(evaluation_cores, length(evaluation_jobs)),
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  lapply(evaluation_jobs, function(job) evaluate_prior(job$prior_id, job$prior_index))
}

forecast_scores <- do.call(rbind, evaluation_results)
rownames(forecast_scores) <- NULL

if (anyNA(forecast_scores$crps)) {
  stop("At least one predictive score could not be computed.")
}

variable_metrics <- do.call(
  rbind,
  lapply(split(forecast_scores, list(forecast_scores$prior_id, forecast_scores$variable)), function(x) {
    data.frame(
      prior_id = x$prior_id[1L],
      variable = x$variable[1L],
      rmse = sqrt(mean(x$squared_error)),
      mae = mean(x$absolute_error),
      mean_crps = mean(x$crps),
      coverage_68 = mean(x$coverage_68),
      coverage_90 = mean(x$coverage_90),
      stringsAsFactors = FALSE
    )
  })
)
rownames(variable_metrics) <- NULL

overall_metrics <- do.call(
  rbind,
  lapply(split(forecast_scores, forecast_scores$prior_id), function(x) {
    data.frame(
      prior_id = x$prior_id[1L],
      standardized_rmse = sqrt(mean(x$squared_error / x$scale^2)),
      standardized_mae = mean(x$absolute_error / x$scale),
      standardized_crps = mean(x$crps / x$scale),
      coverage_68 = mean(x$coverage_68),
      coverage_90 = mean(x$coverage_90),
      stringsAsFactors = FALSE
    )
  })
)
rownames(overall_metrics) <- NULL

overall_metrics$rank_rmse <- rank(overall_metrics$standardized_rmse, ties.method = "min")
overall_metrics$rank_mae <- rank(overall_metrics$standardized_mae, ties.method = "min")
overall_metrics$rank_crps <- rank(overall_metrics$standardized_crps, ties.method = "min")
overall_metrics$mean_rank <- rowMeans(overall_metrics[, c("rank_rmse", "rank_mae", "rank_crps")])
overall_metrics <- merge(catalog, overall_metrics, by = "prior_id", sort = FALSE)
overall_metrics <- overall_metrics[order(overall_metrics$mean_rank, overall_metrics$standardized_crps), ]
overall_metrics$performance_rank <- seq_len(nrow(overall_metrics))
rownames(overall_metrics) <- NULL

winner_id <- overall_metrics$prior_id[1L]
comparison_id <- if (winner_id == "published_baseline") {
  overall_metrics$prior_id[overall_metrics$prior_id != "published_baseline"][1L]
} else {
  winner_id
}

utils::write.csv(
  forecast_scores,
  file.path(output_dir, "forecast_scores.csv"),
  row.names = FALSE
)
utils::write.csv(
  variable_metrics,
  file.path(output_dir, "forecast_metrics_by_variable.csv"),
  row.names = FALSE
)
utils::write.csv(
  overall_metrics,
  file.path(output_dir, "prior_comparison.csv"),
  row.names = FALSE
)

message("Fitting each prior specification on the full sample for diagnostics.")

full_fits <- vector("list", length(model_ids))
names(full_fits) <- model_ids
diagnostic_rows <- vector("list", length(model_ids))

for (model_index in seq_along(model_ids)) {
  prior_id <- unname(model_ids[model_index])
  fit_model <- function(n_draw, n_burn, seed) {
    set.seed(seed)
    BVAR::bvar(
      Y,
      lags = lags,
      n_draw = n_draw,
      n_burn = n_burn,
      n_thin = 1L,
      priors = make_prior(prior_id),
      mh = mh,
      verbose = FALSE
    )
  }

  diagnose_fit <- function(fit) {
    hyper_draws <- coda::as.mcmc(fit, vars = c("lambda", "alpha"))
    geweke <- coda::geweke.diag(hyper_draws)$z
    list(
      hyper_draws = hyper_draws,
      geweke = geweke,
      geweke_p = stats::pnorm(abs(geweke), lower.tail = FALSE) * 2,
      effective_size = coda::effectiveSize(hyper_draws)
    )
  }

  fitted_draws <- settings$full_draws
  fitted_burn <- settings$full_burn
  fit <- fit_model(fitted_draws, fitted_burn, 300000L + model_index)
  diagnostics <- diagnose_fit(fit)

  if (mode == "release" && any(diagnostics$geweke_p < 0.05)) {
    fitted_draws <- 2L * settings$full_draws
    fitted_burn <- 2L * settings$full_burn
    message(
      "Geweke diagnostic flagged ", catalog$label[catalog$prior_id == prior_id],
      "; refitting with ", fitted_draws, " draws."
    )
    fit <- fit_model(fitted_draws, fitted_burn, 300000L + model_index)
    diagnostics <- diagnose_fit(fit)
  }

  hyper_draws <- diagnostics$hyper_draws
  geweke <- diagnostics$geweke
  geweke_p <- diagnostics$geweke_p
  effective_size <- diagnostics$effective_size

  diagnostic_rows[[model_index]] <- data.frame(
    prior_id = prior_id,
    hyperparameter = colnames(hyper_draws),
    posterior_median = unname(apply(hyper_draws, 2L, stats::median)),
    effective_sample_size = unname(as.numeric(effective_size[colnames(hyper_draws)])),
    geweke_z = unname(as.numeric(geweke[colnames(hyper_draws)])),
    geweke_p_value = unname(as.numeric(geweke_p[colnames(hyper_draws)])),
    n_draw = fitted_draws,
    n_burn = fitted_burn,
    acceptance_rate = fit$meta$accepted / fit$meta$n_save,
    stable_draw_share = posterior_stability_rate(fit),
    median_log_likelihood = stats::median(fit$ml),
    stringsAsFactors = FALSE
  )

  if (prior_id %in% unique(c("published_baseline", comparison_id))) {
    full_fits[[prior_id]] <- fit
  }
}

posterior_diagnostics <- do.call(rbind, diagnostic_rows)
rownames(posterior_diagnostics) <- NULL
utils::write.csv(
  posterior_diagnostics,
  file.path(output_dir, "posterior_diagnostics.csv"),
  row.names = FALSE
)

catalog_lookup <- setNames(catalog$label, catalog$prior_id)

grDevices::png(
  file.path(output_dir, "prior-performance.png"),
  width = 1600,
  height = 1000,
  res = 160
)
old_par <- graphics::par(mar = c(5, 16, 4, 2), las = 1)
plot_order <- rev(seq_len(nrow(overall_metrics)))
bar_colors <- ifelse(
  overall_metrics$prior_id[plot_order] == winner_id,
  "#2f8f5b",
  "#5f666d"
)
midpoints <- graphics::barplot(
  overall_metrics$standardized_crps[plot_order],
  names.arg = overall_metrics$label[plot_order],
  horiz = TRUE,
  col = bar_colors,
  border = NA,
  xlab = "Standardized CRPS (lower is better)",
  main = "Prior sensitivity: out-of-sample performance",
  cex.names = 0.9,
  xlim = c(0, max(overall_metrics$standardized_crps) * 1.13)
)
graphics::text(
  overall_metrics$standardized_crps[plot_order],
  midpoints,
  labels = sprintf(" %.3f", overall_metrics$standardized_crps[plot_order]),
  pos = 4,
  cex = 0.9
)
graphics::par(old_par)
grDevices::dev.off()

if (run_structural) {
  message(
    "Computing sign-restricted IRFs for the published baseline and ",
    catalog_lookup[[comparison_id]], "."
  )

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

  irf_options <- BVAR::bv_irf(
    sign_restr = restrictions,
    horizon = 24L,
    fevd = FALSE,
    identification = TRUE,
    sign_lim = 60000L
  )

  structural_ids <- unique(c("published_baseline", comparison_id))
  structural_irfs <- setNames(vector("list", length(structural_ids)), structural_ids)

  for (structural_index in seq_along(structural_ids)) {
    prior_id <- structural_ids[structural_index]
    retained_draws <- full_fits[[prior_id]]$meta$n_save
    irf_thin <- max(1L, floor(retained_draws / settings$structural_draws))
    set.seed(400000L + structural_index)
    structural_irfs[[prior_id]] <- BVAR::irf(
      full_fits[[prior_id]],
      irf_options,
      conf_bands = c(0.05, 0.16),
      n_thin = irf_thin,
      verbose = FALSE
    )
  }

  saveRDS(
    structural_irfs,
    file.path(output_dir, "structural_irfs.rds"),
    compress = "xz"
  )

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
    dim(values) <- c(dim(irf_object$quants)[1L], dim(irf_object$quants)[3L])
    values
  }

  grDevices::png(
    file.path(output_dir, "irf-sensitivity.png"),
    width = 1800,
    height = 1300,
    res = 170
  )
  old_par <- graphics::par(mfrow = c(2, 2), mar = c(4, 4, 4, 1), oma = c(1, 1, 3, 1))
  horizons <- 0:23

  for (path_id in seq_len(nrow(response_paths))) {
    baseline_values <- extract_irf(
      structural_irfs[["published_baseline"]],
      response_paths$impulse[path_id],
      response_paths$response[path_id]
    )
    comparison_values <- extract_irf(
      structural_irfs[[comparison_id]],
      response_paths$impulse[path_id],
      response_paths$response[path_id]
    )
    y_range <- range(
      baseline_values[c(1L, 5L), ],
      comparison_values[c(1L, 5L), ],
      0,
      finite = TRUE
    )

    graphics::plot(
      horizons,
      baseline_values[3L, ],
      type = "n",
      ylim = y_range,
      xlab = "Months",
      ylab = "Response",
      main = response_paths$title[path_id]
    )
    graphics::abline(h = 0, col = "#a8adb2", lty = 2)
    graphics::polygon(
      c(horizons, rev(horizons)),
      c(baseline_values[2L, ], rev(baseline_values[4L, ])),
      border = NA,
      col = grDevices::adjustcolor("#5f666d", alpha.f = 0.16)
    )
    graphics::polygon(
      c(horizons, rev(horizons)),
      c(comparison_values[2L, ], rev(comparison_values[4L, ])),
      border = NA,
      col = grDevices::adjustcolor("#2f8f5b", alpha.f = 0.18)
    )
    graphics::lines(horizons, baseline_values[3L, ], col = "#5f666d", lwd = 2)
    graphics::lines(horizons, comparison_values[3L, ], col = "#2f8f5b", lwd = 2)

    if (path_id == 1L) {
      graphics::legend(
        "topright",
        legend = c("Published baseline", "MN: random-walk mean"),
        col = c("#5f666d", "#2f8f5b"),
        lwd = 3,
        cex = 0.8,
        bty = "n"
      )
    }
  }

  graphics::mtext("Structural sensitivity under the published sign restrictions", outer = TRUE, cex = 1.25, font = 2)
  graphics::par(old_par)
  grDevices::dev.off()
}

session_lines <- capture.output(utils::sessionInfo())
writeLines(session_lines, file.path(output_dir, "session-info.txt"))

summary_lines <- c(
  paste0("Mode: ", mode),
  paste0("Holdout: ", format(dates[holdout_start]), " to ", format(tail(dates, 1L))),
  paste0("One-step forecasts: ", settings$holdout),
  paste0("Winner: ", catalog_lookup[[winner_id]]),
  paste0("Winner standardized CRPS: ", sprintf("%.6f", overall_metrics$standardized_crps[1L])),
  paste0("Published baseline rank: ", overall_metrics$performance_rank[overall_metrics$prior_id == "published_baseline"]),
  paste0("Structural comparison generated: ", run_structural)
)
writeLines(summary_lines, file.path(output_dir, "release-summary.txt"))

message("Prior-sensitivity exercise completed.")
message("Winner: ", catalog_lookup[[winner_id]])
print(overall_metrics[, c(
  "performance_rank",
  "label",
  "standardized_rmse",
  "standardized_mae",
  "standardized_crps",
  "coverage_68",
  "coverage_90"
)], row.names = FALSE)
