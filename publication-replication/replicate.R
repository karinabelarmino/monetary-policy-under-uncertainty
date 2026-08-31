# Monetary Policy Under Uncertainty
# Central Bank Review publication replication with sign and zero restrictions.
# Run the sections in order.

# PACKAGES -----------------------------------------------------------------------------

required_packages <- c("vars", "BVAR", "coda")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the required packages before running the replication: ",
    paste(missing_packages, collapse = ", "),
    "."
  )
}

library(vars)
library(BVAR)
library(coda)

if (as.character(utils::packageVersion("BVAR")) != "1.0.5") {
  warning("This replication was prepared with BVAR 1.0.5.")
}

# LOADING DATABASE ---------------------------------------------------------------------

# Run the script from the repository root so the relative path remains valid.
load("data/monetary_policy_data.RData")

# PRIOR DEFINITIONS --------------------------------------------------------------------

set.seed(1)

mn <- bv_minnesota(
  lambda = bv_lambda(mode = 0.5, sd = 0.4, min = 0.0001, max = 5),
  alpha = bv_alpha(mode = 1),
  var = 1e07
)

# This function creates the sum-of-coefficients dummy prior.
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

soc <- bv_dummy(
  mode = 1,
  sd = 1,
  min = 0.0001,
  max = 50,
  fun = add_soc
)

# This function creates the single-unit-root dummy prior.
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

sur <- bv_dummy(
  mode = 1,
  sd = 1,
  min = 0.0001,
  max = 50,
  fun = add_sur
)

priors <- bv_priors(
  hyper = c("lambda", "alpha"),
  mn = mn,
  sur = sur,
  soc = soc
)

mh <- bv_metropolis(adjust_acc = TRUE, adjust_burn = 0.50)

opt_irf <- bv_irf(
  horizon = 12,
  fevd = TRUE,
  identification = TRUE
)

# VAR SELECTION ------------------------------------------------------------------------

ord <- as.matrix(data[, c("igi", "ci", "gha", "exp", "selic")])

lag_selection <- VARselect(ord, lag.max = 8, type = "const")

message("This print shows the lag order selected by each information criterion.")
print(lag_selection$selection)

# BVAR ESTIMATION ----------------------------------------------------------------------

IC <- bvar(
  ord,
  lags = 2,
  n_draw = 60000,
  n_burn = 15000,
  n_thin = 1,
  priors = priors,
  mh = mh,
  verbose = TRUE
)

plot(IC)

# CONVERGENCE DIAGNOSTICS --------------------------------------------------------------

MCIC <- as.mcmc(IC, vars = c("lambda", "alpha"))
GDMCIC <- geweke.diag(MCIC)

message("This print reports the Geweke convergence diagnostic for lambda and alpha.")
print(GDMCIC)

geweke_p_values <- pnorm(abs(GDMCIC$z), lower.tail = FALSE) * 2

message("This print reports the two-sided p-values associated with the Geweke diagnostic.")
print(geweke_p_values)

# SIGN AND ZERO RESTRICTIONS -----------------------------------------------------------

rs <- matrix(
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

message("This print shows the contemporaneous sign and zero restrictions used for identification.")
print(rs)

sinal <- bv_irf(
  sign_restr = rs,
  horizon = 24,
  fevd = TRUE,
  identification = TRUE,
  sign_lim = 60000
)

message("This print summarizes the structural identification and IRF settings.")
print(sinal)

# BASELINE MODEL -----------------------------------------------------------------------

irfs <- BVAR::irf(
  IC,
  sinal,
  conf_bands = c(0.05, 0.16),
  verbose = TRUE
)

# RESPONSES TO AN UNCERTAINTY SHOCK ----------------------------------------------------

plot(irfs, area = TRUE, vars_impulse = c("igi"), vars_response = c("igi"))
plot(irfs, area = TRUE, vars_impulse = c("igi"), vars_response = c("ci"))
plot(irfs, area = TRUE, vars_impulse = c("igi"), vars_response = c("gha"))
plot(irfs, area = TRUE, vars_impulse = c("igi"), vars_response = c("exp"))
plot(irfs, area = TRUE, vars_impulse = c("igi"), vars_response = c("selic"))

# RESPONSES TO A CREDIBILITY SHOCK -----------------------------------------------------

plot(irfs, area = TRUE, vars_impulse = c("ci"), vars_response = c("ci"))
plot(irfs, area = TRUE, vars_impulse = c("ci"), vars_response = c("igi"))
plot(irfs, area = TRUE, vars_impulse = c("ci"), vars_response = c("gha"))
plot(irfs, area = TRUE, vars_impulse = c("ci"), vars_response = c("exp"))
plot(irfs, area = TRUE, vars_impulse = c("ci"), vars_response = c("selic"))

# RESPONSES TO AN OUTPUT-GAP SHOCK -----------------------------------------------------

plot(irfs, area = TRUE, vars_impulse = c("gha"), vars_response = c("gha"))
plot(irfs, area = TRUE, vars_impulse = c("gha"), vars_response = c("igi"))
plot(irfs, area = TRUE, vars_impulse = c("gha"), vars_response = c("ci"))
plot(irfs, area = TRUE, vars_impulse = c("gha"), vars_response = c("exp"))
plot(irfs, area = TRUE, vars_impulse = c("gha"), vars_response = c("selic"))

# RESPONSES TO AN INFLATION-EXPECTATIONS SHOCK ----------------------------------------

plot(irfs, area = TRUE, vars_impulse = c("exp"), vars_response = c("exp"))
plot(irfs, area = TRUE, vars_impulse = c("exp"), vars_response = c("igi"))
plot(irfs, area = TRUE, vars_impulse = c("exp"), vars_response = c("ci"))
plot(irfs, area = TRUE, vars_impulse = c("exp"), vars_response = c("gha"))
plot(irfs, area = TRUE, vars_impulse = c("exp"), vars_response = c("selic"))

# RESPONSES TO A MONETARY-POLICY SHOCK -------------------------------------------------

plot(irfs, area = TRUE, vars_impulse = c("selic"), vars_response = c("selic"))
plot(irfs, area = TRUE, vars_impulse = c("selic"), vars_response = c("igi"))
plot(irfs, area = TRUE, vars_impulse = c("selic"), vars_response = c("ci"))
plot(irfs, area = TRUE, vars_impulse = c("selic"), vars_response = c("gha"))
plot(irfs, area = TRUE, vars_impulse = c("selic"), vars_response = c("exp"))
