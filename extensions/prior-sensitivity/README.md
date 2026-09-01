# Prior Sensitivity

This post-publication exercise asks two different questions about the exact five-variable, two-lag BVAR used in the publication workflow:

1. How sensitive is the model to changes within the Minnesota-prior specification?
2. What changes when the Minnesota prior is replaced by a different prior family?

Keeping these questions separate avoids treating changes to the Minnesota mean or dummy priors as if they were evidence about non-Minnesota priors.

## Why priors matter

A prior is a probability distribution that records how plausible different parameter values are before the current sample is used. The posterior distribution combines that prior information with the likelihood supplied by the data. A prior does not mechanically determine the result by itself, but it can materially affect regularization and posterior uncertainty when a VAR contains many coefficients relative to the available time observations.

The published baseline combines a hierarchical Minnesota prior with sum-of-coefficients (SOC) and single-unit-root (SUR) dummy priors. The data, five-variable ordering, two lags, Gaussian likelihood, constant covariance matrix and structural restrictions are held fixed throughout this exercise.

## Figure 1: how to read the Minnesota-prior illustration

<p align="center">
  <img
    src="../../assets/minnesota-prior-explainer.png"
    width="500"
    alt="Intuitive explanation of Minnesota-prior shrinkage"
  >
</p>

<p align="center">
  <strong>Figure 1.</strong> Schematic intuition behind Minnesota-prior shrinkage in a BVAR.<br>
  <em>Source: Author's elaboration.</em>
</p>

Each row represents one VAR equation and each group of columns represents coefficients associated with a lag. Within the first-lag block, the green diagonal cells are the first own lags: the first lag of variable 1 in the equation for variable 1, the first lag of variable 2 in its own equation, and so on. The Minnesota prior gives these coefficients more prior freedom. In a random-walk specification they are centered on one; in the zero-mean specification tested here they are centered on zero.

The gray off-diagonal cells are cross-variable effects and are shrunk toward zero. The progressively smaller cells for later lags represent lag decay: distant lags receive stronger shrinkage because they are assumed less likely to carry useful signal.

The cell size and color show **prior freedom**, not estimated coefficient magnitude, posterior significance or coefficient sign. The figure also does not represent the SOC and SUR dummy priors, the Normal-Wishart prior or SSVS; it is only an intuition for the Minnesota component.

## Diagnostic 1: changes within the Minnesota prior

| Specification | Prior mean | SOC | SUR | Interpretation |
| --- | :---: | :---: | :---: | --- |
| Published baseline | 1 | Yes | Yes | Prior configuration used in the publication workflow |
| Minnesota: random-walk mean | 1 | No | No | Own first lags are centered on persistence; other coefficients on zero |
| Minnesota: zero mean | 0 | No | No | All autoregressive coefficients are centered on zero |
| Minnesota + SOC | 1 | Yes | No | Adds information about preservation of long-run levels |
| Minnesota + SUR | 1 | No | Yes | Adds information associated with unit-root-like behavior |

All five specifications retain the hierarchical `lambda` and `alpha` hyperparameters used in the baseline.

## Diagnostic 2: other prior families

Two genuinely non-Minnesota alternatives are estimated with a direct Gibbs sampler in the same script:

| Prior | Main assumption | Role in the comparison |
| --- | --- | --- |
| Independent Normal-Wishart | Diffuse zero-centered Normal priors for standardized coefficients, independent of an inverse-Wishart prior for the covariance matrix | Weak regularization without the Minnesota own-lag and lag-decay structure |
| SSVS spike-and-slab | Every autoregressive coefficient can switch probabilistically between a tightly concentrated spike near zero and a diffuse slab | Data-informed sparsity rather than Minnesota-style lag decay |

The semiautomatic SSVS scale follows the approach documented in the [`bvartools` SSVS vignette](https://cran.r-project.org/web/packages/bvartools/vignettes/ssvs.html), based on [George, Sun and Ni (2008)](https://doi.org/10.1016/j.jeconom.2007.08.017). The implementation is written directly in `prior_sensitivity.R` so the whole exercise remains in one sequential script and does not introduce another model specification or a hidden helper workflow.

The `BVAR` package exposes Minnesota and dummy-observation priors through its public estimation interface. For that reason, the non-Minnesota posterior is sampled directly and then converted to the same reduced-form layout used by `BVAR` for the structural-identification stage.

## Validation design

The ranking uses 36 expanding one-step-ahead forecasts from July 2019 to June 2022. Every forecast origin uses only information available up to that date.

The three scores are standardized by the pre-holdout standard deviation of each variable before aggregation:

- RMSE emphasizes larger point-forecast errors.
- MAE is less sensitive to individual large errors.
- CRPS evaluates the complete posterior predictive distribution; lower values are better.

The ranking is defined separately for each diagnostic as the mean rank across RMSE, MAE and CRPS. Prior selection is completed before the structural IRFs are inspected.

## Results: Diagnostic 1

| Specification | Standardized RMSE | Standardized MAE | Standardized CRPS | Diagnostic rank |
| --- | ---: | ---: | ---: | ---: |
| **Minnesota: zero mean** | **0.427** | 0.242 | **0.207** | **1** |
| Minnesota: random-walk mean | 0.432 | **0.241** | 0.207 | 2 |
| Minnesota + SUR | 0.438 | 0.242 | 0.209 | 3 |
| Minnesota + SOC | 0.452 | 0.249 | 0.214 | 4 |
| Published baseline: MN + SOC + SUR | 0.455 | 0.249 | 0.215 | 5 |

There is no winner in every point-forecast criterion within the Minnesota family. The zero-mean Minnesota has the lowest RMSE and CRPS and ranks first overall, while the random-walk mean has the lowest MAE. Relative to the published baseline, the zero-mean alternative reduces standardized CRPS by 3.9% and standardized RMSE by 6.3%.

## Results: Diagnostic 2

| Specification | Standardized RMSE | Standardized MAE | Standardized CRPS | Diagnostic rank |
| --- | ---: | ---: | ---: | ---: |
| **SSVS spike-and-slab** | **0.427** | **0.231** | **0.184** | **1** |
| Independent Normal-Wishart | 0.430 | 0.242 | 0.188 | 2 |
| Published baseline: MN + SOC + SUR | 0.455 | 0.249 | 0.215 | 3 |

The SSVS prior ranks first on all three aggregate metrics. Relative to the baseline, it reduces standardized CRPS by 14.4%, MAE by 7.4% and RMSE by 6.2%. The independent Normal-Wishart also improves standardized CRPS by 12.5%.

This does not mean that SSVS dominates for every variable or objective. Its aggregate gain is concentrated in better density forecasts for uncertainty, inflation expectations and the output gap, while the baseline remains competitive for some series. The result supports SSVS for this forecasting exercise, not a universal replacement of the published prior.

Predictive coverage also changes materially. The baseline covers 26.7% of observations with its 68% intervals and 41.7% with its 90% intervals. SSVS raises those shares to 71.1% and 83.3%, respectively. The 90% coverage remains below nominal, which is relevant because the holdout includes the COVID-19 shock and the beginning of a rapid monetary-tightening cycle.

<p align="center">
  <img
    src="outputs/prior-performance.png"
    width="850"
    alt="Out-of-sample CRPS results for the two prior diagnostics"
  >
</p>

<p align="center">
  <strong>Figure 2.</strong> Standardized CRPS by diagnostic; lower values indicate better probabilistic forecasts.<br>
  <em>Source: Author's elaboration based on the frozen analytical sample.</em>
</p>

## Structural sensitivity

The published sign and zero restrictions are applied to the baseline, the zero-mean Minnesota and SSVS posterior draws. The principal responses retain the same qualitative interpretation: uncertainty raises inflation expectations and contracts the output gap, credibility lowers inflation expectations, and a positive Selic shock lowers inflation expectations on impact.

The posterior medians differ in magnitude and persistence, especially for the output-gap and credibility responses, but their 68% credible regions overlap broadly. In this application, changing the prior improves predictive performance more than it changes the central structural interpretation.

<p align="center">
  <img
    src="outputs/irf-sensitivity.png"
    width="850"
    alt="Impulse-response sensitivity under three prior specifications"
  >
</p>

<p align="center">
  <strong>Figure 3.</strong> Posterior median IRFs and 68% credible regions under the published baseline and both diagnostic leaders.<br>
  <em>Source: Author's elaboration based on the frozen analytical sample.</em>
</p>

## How to run

Open `prior_sensitivity.R` and run the sections in order from the repository root, as in the publication replication script. Alternatively:

```bash
Rscript extensions/prior-sensitivity/prior_sensitivity.R
```

The default `quick_mode <- FALSE` runs the complete 36-month evaluation. To check the code path locally before a long run, change only this line at the beginning of the script:

```r
quick_mode <- TRUE
```

Quick mode uses six forecast origins, fewer posterior draws and skips the structural comparison. Its scores must not be interpreted as release results.

## Generated outputs

| File | Content |
| --- | --- |
| `outputs/diagnostic_comparison.csv` | Separate rankings for the two diagnostic questions |
| `outputs/prior_comparison.csv` | All seven specifications in one audit table |
| `outputs/forecast_scores.csv` | Forecast-origin and variable-level predictions and scores |
| `outputs/posterior_diagnostics.csv` | Effective sample size, Geweke and stability summaries |
| `outputs/prior-performance.png` | CRPS comparison for both diagnostics |
| `outputs/irf-sensitivity.png` | Structural comparison under the published restrictions |
| `outputs/session-info.txt` | R and package environment used for the run |

## Scope and limitations

- The ranking is conditional on this sample, variable ordering, two-lag model and 2019–2022 holdout.
- Standardizing the non-Minnesota systems is an invertible reparameterization used to make prior scales comparable across equations; forecasts and structural draws are transformed back to the original units.
- The package can emit a non-fatal warning when its automatic `psi` AR fit is estimated in a small number of expanding windows. The warning is left visible.
- At 60,000 full-sample draws, the random-walk Minnesota has a maximum absolute Geweke statistic of 2.76. It is not the structural-comparison leader, but this chain should be inspected or rerun for longer before publication if its full-sample posterior is reported directly.
- SSVS inclusion indicators introduce a different form of regularization, not merely a different Minnesota hyperparameter.
- Forecast performance and structural plausibility answer different questions. A forecasting gain does not make the published structural choice “wrong.”
- Sign restrictions identify a set of admissible models rather than one unique structural model.
