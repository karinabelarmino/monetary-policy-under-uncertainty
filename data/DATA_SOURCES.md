# Data sources

`monetary_policy_data.RData` is the updated frozen analytical file supplied for this project. Loading the file creates one data frame named `dados`, with 234 monthly observations from January 2003 to June 2022. The monthly reference dates are stored directly as values of class `Date`.

The file contains eight columns: `date`, `exp`, `ghp`, `gha`, `gap`, `selic`, `igi` and `ci`. It was converted directly from the updated analytical workbook and validated for column order, monthly coverage, duplicate dates, missing values and non-finite numeric observations. The original Excel workbook is not required to run the replication.

The underlying series are publicly available from the Central Bank of Brazil, Ipeadata, IBRE/FGV, FecomercioSP and the Brazilian Economic Policy Uncertainty project. The General Macroeconomic Uncertainty Indicator is derived from EPU, IIE-BR, consumer uncertainty and EMBI+ Brazil. Definitions, transformations and source codes are reported in the published article.

Article: https://doi.org/10.1016/j.cbrev.2026.100257
