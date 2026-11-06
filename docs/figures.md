# Figures

All figures are generated in R from the reproduced result tables. Each is supplied as a vector PDF and a PNG. The figures show deterministic probability sums or model-based expectations.

| File stem | What it shows |
|---|---|
| `gain_by_MTD` | PCS change at each true MTD location for the two 90% tables |
| `gain_and_changes` | PCS gain and frequency of changed recommendations on the independent curves for the full, 80% and 90% tables |
| `prior_sensitivity` | Performance of the shorter tables under central, high and equal location priors |
| `changed_recommendations` | Corrected selections, newly incorrect selections and other changes per thousand modeled trials |
| `expansion_consequences` | Predicted reduction in expansion DLTs under the MN and DD historical sources |

PCS is the probability of selecting the true MTD. The 80% and 90% labels refer to the retained fraction of the full table's planning gain. They are not confidence levels. In the clinical figure, expansion fills to eighteen total patients at the selected dose and retains the specified safety checks. MN denotes the metronomic source and DD the dose-dense source. Both use historical likelihood power 0.5 in this figure.
