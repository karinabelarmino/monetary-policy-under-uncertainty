add_soc <- function(Y, lags, par) {
  soc <- if (lags == 1L) {
    diag(Y[1L, ]) / par
  } else {
    diag(colMeans(Y[seq_len(lags), , drop = FALSE])) / par
  }

  list(
    Y = soc,
    X = cbind(
      rep(0, ncol(Y)),
      matrix(rep(soc, lags), nrow = ncol(Y))
    )
  )
}

add_sur <- function(Y, lags, par) {
  sur <- if (lags == 1L) {
    Y[1L, ] / par
  } else {
    colMeans(Y[seq_len(lags), , drop = FALSE]) / par
  }

  list(Y = sur, X = c(1 / par, rep(sur, lags)))
}

make_minnesota <- function(prior_mean = 1) {
  BVAR::bv_minnesota(
    lambda = BVAR::bv_lambda(
      mode = 0.5,
      sd = 0.4,
      min = 0.0001,
      max = 5
    ),
    alpha = BVAR::bv_alpha(mode = 1),
    var = 1e07,
    b = prior_mean
  )
}

make_soc <- function() {
  BVAR::bv_dummy(
    mode = 1,
    sd = 1,
    min = 0.0001,
    max = 50,
    fun = add_soc
  )
}

make_sur <- function() {
  BVAR::bv_dummy(
    mode = 1,
    sd = 1,
    min = 0.0001,
    max = 50,
    fun = add_sur
  )
}

make_prior <- function(prior_id) {
  prior_args <- list(
    hyper = c("lambda", "alpha"),
    mn = make_minnesota(prior_mean = if (prior_id == "mn_zero_mean") 0 else 1)
  )

  if (prior_id %in% c("mn_soc", "published_baseline")) {
    prior_args$soc <- make_soc()
  }

  if (prior_id %in% c("mn_sur", "published_baseline")) {
    prior_args$sur <- make_sur()
  }

  do.call(BVAR::bv_priors, prior_args)
}

prior_catalog <- function() {
  data.frame(
    prior_id = c(
      "published_baseline",
      "mn_random_walk",
      "mn_zero_mean",
      "mn_soc",
      "mn_sur"
    ),
    label = c(
      "Published baseline: MN + SOC + SUR",
      "Minnesota: random-walk mean",
      "Minnesota: zero mean",
      "Minnesota + SOC",
      "Minnesota + SUR"
    ),
    prior_mean = c(1, 1, 0, 1, 1),
    soc = c(TRUE, FALSE, FALSE, TRUE, FALSE),
    sur = c(TRUE, FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
}

crps_sample <- function(draws, observation) {
  draws <- sort(as.numeric(draws))
  n_draws <- length(draws)

  if (n_draws < 2L || !is.finite(observation) || any(!is.finite(draws))) {
    return(NA_real_)
  }

  first_term <- mean(abs(draws - observation))
  weights <- 2 * seq_len(n_draws) - n_draws - 1
  pairwise_term <- 2 * sum(weights * draws) / (n_draws^2)

  first_term - 0.5 * pairwise_term
}

posterior_stability_rate <- function(fit, max_draws = 2000L) {
  beta <- fit$beta
  n_available <- dim(beta)[1L]
  draw_ids <- unique(round(seq(1, n_available, length.out = min(max_draws, n_available))))
  n_vars <- fit$meta$M
  lags <- fit$meta$lags

  stable <- vapply(draw_ids, function(draw_id) {
    coefficients <- beta[draw_id, -1L, , drop = FALSE]
    dim(coefficients) <- c(n_vars * lags, n_vars)
    top_block <- do.call(
      cbind,
      lapply(seq_len(lags), function(lag_id) {
        rows <- ((lag_id - 1L) * n_vars + 1L):(lag_id * n_vars)
        t(coefficients[rows, , drop = FALSE])
      })
    )

    companion <- if (lags == 1L) {
      top_block
    } else {
      rbind(
        top_block,
        cbind(diag(n_vars * (lags - 1L)), matrix(0, n_vars * (lags - 1L), n_vars))
      )
    }

    max(Mod(eigen(companion, only.values = TRUE)$values)) < 1
  }, logical(1))

  mean(stable)
}

get_env_integer <- function(name, default, minimum = 1L) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }

  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < minimum) {
    stop(name, " must be an integer greater than or equal to ", minimum, ".")
  }

  parsed
}

get_env_flag <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  if (!value %in% c("true", "false")) {
    stop(name, " must be either true or false.")
  }

  identical(value, "true")
}
