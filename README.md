# Prior information for final dose selection in BOIN trials

This folder accompanies *The value of prior information for final dose selection in Bayesian optimal interval trials* by Siyi Chen. R handles rule construction, evaluation, sensitivity analyses and figures. Rcpp provides likelihood calculations, endpoint checks and independent risk verification. Optimization uses the R interface to HiGHS.

Explanatory comments have been removed from the R and C++ sources. Required Rcpp export annotations are retained. Figures are generated in R.

## Start here

Use R 4.4 or later. Compiling the Rcpp sources requires a GCC or Clang toolchain with 128-bit integer support. Run these commands from this folder.

```console
Rscript scripts/install_packages.R
Rscript scripts/verify.R
Rscript scripts/reproduce.R --mode all --output ../boin-results
Rscript scripts/plot_results.R --input ../boin-results --output ../boin-figures
```

Each output directory must be new and outside this code folder. The full run rebuilds five trial catalogs, checks partition coverage, verifies the saved rules, recalculates gains and operating characteristics, and reproduces the clinical prior and expansion analyses. It also checks the robust-prior gain bounds, reconstructs the logistic posterior proposals, evaluates the matched methods and supporting designs, reports generator intervals, and replays the ordered-likelihood and CSE comparisons. Results are checked against the supplied reference tables. The plotting command creates five figures, each as a PDF and a PNG, using R graphics.

The supplied [results](results/) and [figures](figures/) show the completed local run. [Reproduction details](docs/reproduction.md) describe the calculations and their scope. The [technical note](docs/technical_note/technical_note.pdf) is a complete, readable document with its LaTeX source. The [result map](docs/result_map.csv) links every article table to its inputs and calculation.

## Look up a recommendation

```console
Rscript scripts/select_dose.R --table central_sparse90 --counts 12,12,3,3,0,0 --dlts 2,2,0,1,0,0
Rscript scripts/select_dose.R --table primary_MN --counts 0,6,12 --dlts 0,1,5
```

The first example uses the six-dose design. The second uses the three-dose clinical design. The program rejects count patterns outside the chosen protocol's terminal catalog. A recommendation of `-1` means no dose is selected. The tables require the conduct, stopping rules and prespecified planning assumptions used to construct them.

## Fit a rule

```console
Rscript scripts/fit_rules.R --prior high --mode short --retention 0.9 --output ../boin-fit
Rscript scripts/fit_rules.R --prior central --mode short --budget 0.01 --output ../boin-budget-fit
Rscript scripts/fit_clinical.R --rule primary_MN --output ../boin-clinical-fit
Rscript scripts/fit_secondary.R --mode robust --prior joint --output ../boin-family-fit
Rscript scripts/fit_secondary.R --mode matched --prior central --proposal quantile --seconds 8 --rounds 6 --output ../boin-logit-fit
Rscript scripts/fit_secondary.R --mode supporting --design N36 --prior high --seconds 8 --rounds 8 --output ../boin-supporting-fit
```

The fitter supports a full table, a shorter table and the endpoint flow construction. Risk budgets are 0.005, 0.01, 0.02 and 0.03. The clinical fitter supports all six clinical specifications. Secondary fitting supports high or joint location/shape uncertainty, matched posterior proposals and all three supporting designs. Each newly saved rule passes a separate Rcpp risk check. Solver versions, time limits and ties can produce different feasible tables. The distributed tables define the manuscript results.

The secondary fitter starts with the saved feasible rules by default. Use `--start boin` to start with BOIN alone. The robust search uses 60 seconds per round by default. The logistic command above uses the article's six eight-second limits. The supporting designs are `N36`, `J8` and `cohort2`. Matched proposals are `quantile`, `weak`, `mode` and `mean`.

## Run the supporting analyses separately

```console
Rscript scripts/reproduce_secondary.R --mode all --output ../boin-secondary-results
Rscript scripts/reproduce_secondary.R --mode benchmark --output ../boin-cse-benchmark
```

The first command verifies the saved supporting results. Its CSE replay uses the original refinement counts. The second repeats the 60-second, 200,000-leaf comparison on the current computer. Fresh wall-clock benchmarks can stop at different refinement counts. Original construction and benchmark timings are retained separately in `data/secondary/timing/`.

## Folder contents

| Folder | Contents |
|---|---|
| `R/` | BOIN conduct, input generators, fitting, evaluation, gain audits and plotting |
| `src/` | Rcpp likelihood calculations and independent integer checks |
| `scripts/` | R command line entry points |
| `data/` | Frozen planning inputs, decision tables and reference results |
| `results/` | Reproduced tables from the current R/Rcpp run |
| `figures/` | Figures generated by the R plotting script |
| `verification/` | Local check results and tested package versions |
| `tests/` | Reachable terminal-data examples |

The data files contain modeled curves and published aggregate counts. The clinical expansion results are planning predictions. The work has not been prospectively evaluated as a change to clinical practice.

## Citation and license

Use [CITATION.cff](CITATION.cff) to cite this code. The accompanying manuscript is unpublished. The code is distributed under the [MIT license](LICENSE).
