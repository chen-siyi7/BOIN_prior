# Reproducing the calculations

## Inputs and saved rules

`data/core.rds` contains the six-dose model, planning priors, ordered partition, endpoint graph and saved rules. `data/clinical.rds` contains the three-dose model, historical priors, preferences, risk specifications and six clinical rules. `data/supporting.rds` supplies the three additional trial designs. `data/evaluation.rds` contains the independent evaluation curves. Gain-ceiling prices and candidate sets are in `data/gain_ceilings.rds`. `data/secondary.rds` supplies the robust and matched-method rules, conditional-family metadata, robust prior and risk prices, logistic model parameters, and supporting partitions, rules and evaluation libraries.

The RDS files preserve the numerical inputs used in the study. The replay regenerates all five terminal catalogs from the BOIN conduct rules. It also regenerates every evaluation library used here, the supporting logistic planning curves, the clinical planning sample and all seven historical priors. CSV evaluation inputs are checked within their decimal round-trip precision. Binary RDS inputs are compared directly. The original six-dose planning support remains a supplied input. The script checks that all supplied ordered partitions cover their starting cells. A discarded child is accepted only when it contains no ordered curve. The eight-dose partition starts from five intervals per coordinate, while the others use eight.

## Verification

`scripts/verify.R` checks the file inventory, parses the R sources and tests dose lookup. It requires only `jsonlite` and `digest`.

`scripts/reproduce.R --mode all` performs the scientific replay. It uses `Rcpp` and `gmp` in addition to the packages above.

| Calculation | Implementation |
|---|---|
| Uniform risk bounds | Independent outward-rounded Q63 arithmetic in Rcpp |
| Endpoint flow certificate | Integer capacity checks and maximum flow in Rcpp |
| Planning gains and table-size checks | Exact rational arithmetic in R with `gmp` |
| Clinical score gains | Exact rational arithmetic over the saved finite priors |
| Nominal gain ceilings | Supplied nonnegative prices, independently enclosed likelihoods and rational weak-duality calculations |
| Robust-prior gains and ceilings | Conditional-family probability enclosures, exact prior normalization, membership checks and rational weak-duality calculations |
| Matched logistic proposals | Refit of the calibrated priors, Gaussian posterior quadrature at 181 and 241 nodes per dimension, action-stability and independent marginal-CDF checks |
| Supporting trial designs | Full saved rules, partition coverage, integer risk checks, planning gains and terminal-outcome evaluation |
| Generator intervals | Stratified paired standard errors across four independent curve families at each MTD position |
| Risk-bound comparisons | Global adverse, ordered adverse, coordinatewise signed and ordered signed calculations for the same fixed rules |
| CSE comparison | Both independent-sequence and bounded-stopping moment embeddings, signed and affine tilt bounds, original refinement replay and fresh benchmarking |
| Evaluation outcomes | Terminal probability sums on the fixed independent curves |
| Continuous historical priors | R Gauss-Legendre integration at three orders |
| Expansion predictions | R backward recursion, checked against a separate exact forward recursion |

Saved binary64 probabilities and weights are interpreted exactly in the rational checks. Clinical weights are normalized exactly for the gain audits. The integer risk checks reconstruct likelihood bounds and check the required label and toxicity-window masks. Floating-point calculations propose rules and report descriptive outcomes. A proposed rule must pass the independent risk check before it is returned.

The gain-ceiling check validates the candidate universe and the domain of each priced risk constraint. Nonnegative prices give a weak-duality bound without requiring an optimal numerical solve. The Rcpp likelihood enclosures can give slightly wider bounds than the archived high-precision calculation. Both values are retained in the output. These are bounds on the stated planning objectives and risk specifications.

The full replay includes the nominal bounds used in the clinical benefit thresholds and both robust-prior bounds. The supporting analyses are executable parts of the replay. Original optimizer histories and timing records describe the searches that produced the frozen rules. New fitting commands construct and verify rules from those supplied design specifications. A time-limited refit need not return the original optimizer path or the same table.

## Fitting

Rule fitting requires `highs` and `Matrix`. The R interface calls HiGHS directly. `scripts/fit_rules.R` accepts `--mode short`, `--mode full` or `--mode transport`. `--seconds` is the solver limit for each constraint round. `--retention` is the fraction of the supplied full table's planning gain sought by the shorter table.

`scripts/fit_clinical.R --rule NAME` accepts `primary_MN`, `primary_DD`, `expanded_MN`, `expanded_DD`, `uncertain_MN` or `uncertain_sources`. The last two fit the worst gain across the specified prior or preference vertices. The frozen clinical priors and risk specifications are shared with the verification code.

`scripts/fit_secondary.R` supports `--mode robust`, `--mode matched` and `--mode supporting`. Robust fitting uses one auxiliary value per MTD position and a further worst-case objective over location vertices. Matched fitting restricts departures to the chosen posterior proposal. Supporting fitting uses the corresponding trial catalog and partition. All modes perform a separate integer risk check before returning a rule. `--start saved` retains the supplied feasible incumbents. `--start boin` starts from BOIN alone.

## CSE and timing

The default CSE replay reconstructs the original ordered-likelihood and tilt bounds by repeating refinement to the original final leaf counts. It checks these values against the supplied original results. The current replay time is measured separately. `scripts/reproduce_secondary.R --mode benchmark` instead applies a fresh 60-second budget and 200,000-leaf cap to each method after initialization.

The CSE results use numerical moment and tilt calculations. They do not receive the integer certification used for the ordered-likelihood rule. The completed ordered-likelihood comparison receives its own partition and integer risk checks. Independent direct terminal sums also check the stopped-trial likelihood-ratio moment bounds.

Original phase boundaries, times and peak memory records are in `data/secondary/timing/`. They distinguish nominal construction from later analyses that start with frozen inputs. The new replay writes `timing.csv` and `summary.json`. Local hardware and the availability of resource measurements are recorded in `verification/`. The original workflow and the new replay have different scopes, so their durations should not be substituted for one another.

## Figures

All figures are drawn by `R/plots.R` through `scripts/plot_results.R`. They use only R graphics. No external plotting program is called. See [figure descriptions](figures.md).

## Local records

The supplied `verification/` folder records the local R version, package versions, full replay and example fits. GitHub Actions runs the inventory and lookup checks on pushes. The full scientific replay is also available as a manual workflow. The [result map](result_map.csv) identifies the files that reproduce each reported table. The [technical note](technical_note/technical_note.pdf) gives model specifications and the supporting comparisons.
