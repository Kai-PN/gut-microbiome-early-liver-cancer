# Gut bacterial and fungal signatures for early HCC detection

R Markdown pipeline for the final classification model in our study of cross-kingdom (bacterial + fungal) gut dysbiosis in early hepatocellular carcinoma (HCC), using shotgun metagenomic data. It trains HCC-vs-control classifiers on bacterial (`bac`), fungal (`fun`) and combined (`bacfun`) species profiles, validates them on three external cohorts (China, Thailand, Korea), and explains the final models with SHAP.

## Pipeline overview

All steps are in `Final_Model_AdaBoost_N20.Rmd`.

1. **Data preparation** 
– species-level aggregation of discovery and validation phyloseq objects; bacterial and fungal tables merged for the `bacfun` model.
2. **Filtering and batch correction** 
– taxa kept if present in ≥10% of samples; relative abundances batch-corrected with MMUPHin (`Batch` as batch variable, `Group` as covariate).
3. **Validation alignment** 
– external-cohort taxa matched to discovery taxa by species name; unmatched discovery taxa set to zero.
4. **Repeated training (20 repeats)** 
– each repeat: stratified 80/20 split, down-sampling of the majority class, near-zero-variance removal, AdaBoost.M1 (decision stumps, 50 iterations) for feature ranking, top 20 features passed to a random forest tuned by 5-fold × 20-repeat cross-validation (ROC metric).
5. **Evaluation** 
– AUC, accuracy, sensitivity, specificity and F1 on the hold-out set and each external cohort, summarized as mean with t-based 95% CI across repeats; mean ROC curves vertically averaged across repeats.
6. **Feature selection frequency** 
– how often each feature was selected across the 20 repeats, annotated with taxonomy.
7. **Final refit and SHAP** 
– random forest refit on the full discovery data using the 20 consensus features (ranked by selection frequency, then mean importance); per-sample SHAP values computed on the external cohorts with `fastshap` (500 Monte Carlo simulations).

## Requirements

R (v.4.5 recommended) with:

- CRAN: `dplyr`, `tidyr`, `tibble`, `caret`, `pROC`, `randomForest`, `ggplot2`, `scales`, `adabag`, `fastshap`, `ggbeeswarm`, `ggpubr`, `tidytext`
- Bioconductor: `phyloseq`, `MMUPHin`, `microbiome`
- GitHub: `speedyseq` (`remotes::install_github("mikemc/speedyseq")`)

```r
install.packages(c("dplyr", "tidyr", "tibble", "caret", "pROC", "randomForest",
                   "ggplot2", "scales", "adabag", "fastshap", "ggbeeswarm",
                   "ggpubr", "tidytext", "remotes", "BiocManager"))
BiocManager::install(c("phyloseq", "MMUPHin", "microbiome"))
remotes::install_github("mikemc/speedyseq")
```

## Data

Input data are **not included** in this repository because they contain participant-level information. The pipeline expects:

- Discovery cohort
- Validation cohorts (China, Thailand, Korea)

File paths are currently set in the setup and data chunks (`project_dir`, `validation_paths`, and the `readRDS()` / `read.csv()` calls in `prep_data()`). Edit these to match your local directory before running.

## Usage

Open `Final_Model_AdaBoost_N20.Rmd` in RStudio and run the chunks in order, or knit the whole document. Intermediate results are cached as `.rds` files in `checkpoints/` through `run_or_load()`, so re-running skips completed steps; pass `force = TRUE` to recompute a step.

Outputs:

- `performance/` – per-repeat metrics, summary table, raw predictions, mean ROC curves, feature selection frequency, and an `.RData` bundle of all results
- `figures/` – performance summary, ROC curves, feature frequency, SHAP beeswarm plots

Training all three models (20 repeats each, with repeated CV inside every repeat) is computationally heavy; expect a long first run.

## Citation

Manuscript in preparation. Citation details will be added on publication.

## Declaration of generative AI 
The authors used the assistance of Claude and Gemini to finish the code.