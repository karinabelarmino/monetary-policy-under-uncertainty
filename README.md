# Monetary Policy Under Uncertainty

This repository documents a reproducible research agenda on uncertainty, monetary-policy credibility and macroeconomic dynamics in Brazil. The project is organized by analytical purpose: the workflow tied to the published article is kept intact, while post-publication exercises live under `extensions/` and share the same frozen dataset.

The starting point is [*Uncertainty, credibility and monetary policy in Brazil: A BVAR approach with sign restrictions*](https://doi.org/10.1016/j.cbrev.2026.100257), published in *Central Bank Review* in 2026.

## Project map

| Module | Analytical question | Main evidence | Status |
| --- | --- | --- | :---: |
| [Central Bank Review publication replication](replication/) | Can the published five-variable BVAR be transparently reconstructed with the frozen data and `BVAR` 1.0.5? | Baseline estimation, convergence diagnostics and sign-restricted IRFs | Published workflow |
| [Prior sensitivity](extensions/prior-sensitivity/) | How sensitive is the baseline to Minnesota specifications and to non-Minnesota prior families? | Two diagnostics, 36 one-step forecasts, posterior checks and structural sensitivity | Local validation candidate |
| `bsvarSIGNs` assessment | How portable are the model and identification scheme to another sign-restricted BVAR implementation? | Cross-package specification and output comparison | Planned |

## Prior-sensitivity extension

The extension separates two questions that should not be conflated:

1. **Within-Minnesota sensitivity:** prior mean and SOC/SUR dummy-prior choices.
2. **Across-family sensitivity:** independent Normal-Wishart and SSVS spike-and-slab priors.

The data, five-variable system, two lags, Gaussian likelihood and published sign and zero restrictions remain fixed. Across 36 expanding one-step forecasts from July 2019 to June 2022, the zero-mean Minnesota is the strongest Minnesota specification, while SSVS provides the best aggregate performance among all tested alternatives.

Relative to the published baseline, SSVS reduces standardized CRPS by 14.4%, MAE by 7.4% and RMSE by 6.2%. Its predictive intervals also have substantially better empirical coverage. These are conditional forecasting results, not evidence that the prior used for the structural publication was incorrect.

The main impulse-response interpretation remains qualitatively stable across the published baseline, the zero-mean Minnesota and SSVS, with broadly overlapping 68% credible regions.

## Shared model and data

The frozen analytical file contains 234 monthly observations from January 2003 to June 2022. The baseline system orders five endogenous variables as IGI, CI, GHA, EXP and SELIC and uses two lags. Public and derived series come from the Central Bank of Brazil, Ipeadata, IBRE/FGV, FecomercioSP and the Brazilian Economic Policy Uncertainty project.

See [`data/DATA_SOURCES.md`](data/DATA_SOURCES.md) for the data dictionary and source notes.

## Repository structure

```text
.
├── assets/                              Shared explanatory figures
├── data/                                Frozen analytical data and source notes
├── extensions/
│   └── prior-sensitivity/
│       ├── outputs/                     Generated tables and figures
│       ├── prior_sensitivity.R          Both prior diagnostics in one script
│       └── README.md                    Methods, results and interpretation
├── replication/                         Article-linked reconstruction
├── social/                              Draft professional communication
├── CITATION.cff                         Citation metadata
├── DESCRIPTION                          R dependencies
├── LICENSE                              MIT License
├── replicate.R                          Backward-compatible replication entry point
└── README.md                            Project-level navigation
```

This structure is preferable to placing every file in a long publication-named folder: it keeps the root concise, gives the article workflow a stable home, and leaves a predictable `extensions/<exercise>/` path for future work. The root-level `replicate.R` is retained as a backward-compatible entry point.

## Reproducibility

Requirements:

- R 4.3 or newer
- `BVAR` 1.0.5
- `coda`
- `vars`

Run the publication workflow:

```bash
Rscript replication/replicate.R
```

Run both prior diagnostics:

```bash
Rscript extensions/prior-sensitivity/prior_sensitivity.R
```

The complete sensitivity exercise is computationally intensive. For a short local code-path check, set `quick_mode <- TRUE` inside `prior_sensitivity.R`; quick-mode scores are not release results.

## Credits and funding

The published paper was coauthored by Karina Oliveira Belarmino de Almeida, Wilson Luiz Rotatori Corrêa and Luckas Sabioni Lopes.

The research received support from FAPEMIG through a Researcher in Science, Technology and Innovation Development scholarship, BDCTI-IV.

## License

The code is available under the MIT License. Original data providers are credited in the source documentation and in the published article.
