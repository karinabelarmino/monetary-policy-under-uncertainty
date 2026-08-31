required_packages <- c("BVAR")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

if (!file.exists("DESCRIPTION")) {
  stop("Run the smoke test from the repository root.")
}

source("extensions/prior-sensitivity/R/prior_helpers.R")

catalog <- prior_catalog()
stopifnot(nrow(catalog) == 5L)
stopifnot(length(unique(catalog$prior_id)) == nrow(catalog))

constructed_priors <- lapply(catalog$prior_id, make_prior)
stopifnot(all(vapply(constructed_priors, inherits, logical(1), what = "bv_priors")))

known_crps <- crps_sample(c(-1, 0, 1), 0)
stopifnot(isTRUE(all.equal(known_crps, 2 / 9, tolerance = 1e-12)))

load("data/monetary_policy_data.RData")
stopifnot(exists("data"), is.data.frame(data), nrow(data) == 234L)
stopifnot(all(c("date", "igi", "ci", "gha", "exp", "selic") %in% names(data)))
stopifnot(!anyNA(data[, c("igi", "ci", "gha", "exp", "selic")]))

message("Prior-sensitivity smoke test passed.")
