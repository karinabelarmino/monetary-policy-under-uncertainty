# Monetary Policy Under Uncertainty

This repository documents a reproducible research agenda on uncertainty, monetary-policy credibility and macroeconomic dynamics in Brazil. It separates the workflow linked to the published article from post-publication methodological extensions while keeping the same frozen dataset and model context.

The starting point is [*Uncertainty, credibility and monetary policy in Brazil: A BVAR approach with sign restrictions*](https://doi.org/10.1016/j.cbrev.2026.100257), published in *Central Bank Review* in 2026.

## Project map

| Module | Analytical question | Main evidence | Status |
| --- | --- | --- | :---: |
| [Central Bank Review publication replication](publication-replication/) | Can the published five-variable BVAR be transparently reconstructed with the frozen data and `BVAR` 1.0.5? | Baseline estimation, convergence diagnostics and sign-restricted IRFs | Released · Aug 25, 2026 |
| [Prior sensitivity](extensions/prior-sensitivity/) | Does the baseline remain competitive under alternative Minnesota and dummy-prior specifications? | 36 expanding one-step forecasts, full-sample diagnostics and structural comparison | Released · Aug 31, 2026 |
| `bsvarSIGNs` assessment | How portable are the model and identification scheme to another sign-restricted BVAR implementation? | Cross-package specification and output comparison | Planned |

## What the new extension adds

The published script uses a hierarchical Minnesota prior together with sum-of-coefficients (SOC) and single-unit-root (SUR) dummy priors. The sensitivity exercise preserves the five variables, two lags, sample, likelihood and sign restrictions, then changes only the prior specification.

Five configurations are evaluated over an expanding holdout from July 2019 to June 2022. The Minnesota prior with a random-walk mean and no additional dummy priors provides the best overall out-of-sample performance: it ranks first on standardized MAE and CRPS. The zero-mean Minnesota alternative records the lowest standardized RMSE, so the result is reported as a multi-metric ranking rather than a claim that one prior dominates every criterion.

| Specification | Standardized RMSE | Standardized MAE | Standardized CRPS | Overall rank |
| --- | ---: | ---: | ---: | ---: |
| Minnesota: random-walk mean | 0.432 | 0.241 | 0.207 | 1 |
| Minnesota: zero mean | 0.427 | 0.243 | 0.207 | 2 |
| Minnesota + SUR | 0.438 | 0.242 | 0.209 | 3 |
| Minnesota + SOC | 0.451 | 0.249 | 0.214 | 4 |
| Published baseline: Minnesota + SOC + SUR | 0.455 | 0.249 | 0.215 | 5 |

The leading specification reduces standardized CRPS by 3.8% relative to the published baseline over this holdout. This is evidence about predictive performance in a demanding period that includes the COVID-19 shock and the beginning of a rapid monetary-tightening cycle; it is not evidence that the published prior was universally incorrect.

Under the published sign restrictions, the baseline and the leading alternative produce closely aligned median IRFs with substantially overlapping credible regions. The prior choice matters more for the forecast ranking than for the central structural interpretation in this application.

## Shared model and data

The frozen analytical file contains 234 monthly observations from January 2003 to June 2022. The baseline system orders five endogenous variables as IGI, CI, GHA, EXP and SELIC and uses two lags. Public and derived series come from the Central Bank of Brazil, Ipeadata, IBRE/FGV, FecomercioSP and the Brazilian Economic Policy Uncertainty project.

See [`data/DATA_SOURCES.md`](data/DATA_SOURCES.md) for the data dictionary and source notes.

## Repository structure

```text
.
├── assets/                              Shared explanatory figures
├── data/                                Frozen analytical data and source notes
├── extensions/
│   └── prior-sensitivity/               Post-publication prior comparison
├── publication-replication/             Article-linked reconstruction
├── social/                              Release copy for professional communication
├── CITATION.cff                         Citation metadata
├── DESCRIPTION                          R dependencies
├── LICENSE                              MIT License
├── replicate.R                          Backward-compatible replication entry point
└── README.md                            Project-level navigation
```

The root-level `replicate.R` is retained only as a backward-compatible entry point for links from the first release.

## Reproducibility

Requirements:

- R 4.3 or newer
- `BVAR` 1.0.5
- `coda`
- `vars`

Run the published workflow:

```bash
Rscript publication-replication/replicate.R
```

Run the complete prior-sensitivity release:

```bash
Rscript extensions/prior-sensitivity/run_prior_sensitivity.R
```

For a short execution check, use `PRIOR_SENSITIVITY_MODE=quick`. The committed results were generated in `release` mode and should not be replaced with quick-mode output.

## Credits and funding

The published paper was coauthored by Karina Oliveira Belarmino de Almeida, Wilson Luiz Rotatori Corrêa and Luckas Sabioni Lopes.

The research received support from FAPEMIG through a Researcher in Science, Technology and Innovation Development scholarship, BDCTI-IV.

## License

The code is available under the MIT License. Original data providers are credited in the source documentation and in the published article.
