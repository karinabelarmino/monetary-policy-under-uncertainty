# Backward-compatible entry point for the publication-linked replication.
# The maintained script now lives in publication-replication/replicate.R.

if (!file.exists("publication-replication/replicate.R")) {
  stop("Run this entry point from the repository root.")
}

source("publication-replication/replicate.R", chdir = FALSE)
