# metaCpG: Methylation-Based Overall Survival Prediction in Medulloblastoma

A point-and-click tool for estimating overall survival (OS) in medulloblastoma
from a single DNA methylation array. You provide raw IDAT files, choose a
molecular context, and the app returns a **metaCpG risk score**, the sample's
position within a reference cohort, and a **model-based predicted OS**.

> **Research use only.** This tool is intended for research and educational use.
> It is not a medical device and must not be used to guide individual patient
> management. See [Intended use and limitations](#intended-use-and-limitations).

---

## What it predicts

| Molecular context | Output | Covariates in the survival model |
|---|---|---|
| **SHH** | Predicted **5-year OS** | metaCpG score + *MYCN* amplification status |
| **Group 3/4** (early) | Predicted **5-year OS** | metaCpG score + M-stage + stratified by *MYC* amplification |
| **Group 3/4** (late) | Predicted **10-year mortality** (late events) | metaCpG score alone — the *late* metaCpG |

The **late metaCpG** for Group 3/4 is a distinct signature that captures the
risk of late mortality (i.e. probability of 10-year survival), and is reported
alongside the early 5-year estimate.

---

## How it works

1. **Read arrays.** Raw IDATs are imported with `minfi`. The tool auto-detects
   450K vs EPIC arrays by file size and, if a batch contains both, combines them
   onto the common 450K probe space (`combineArrays`).
2. **Preprocess.** `preprocessNoob` (single-dye), detection p-values, then
   `ratioConvert` → `mapToGenome` → beta values.
3. **Score.** Betas are converted to M-values, subset to the signature CpGs, and
   collapsed to a single weighted-mean **metaCpG score**, which is then z-scaled
   to the discovery cohort (stored centre/scale).
4. **Locate.** The score is placed on the reference-cohort distribution so you
   can see where a sample sits relative to previously analysed cases.
5. **Predict OS.** The score (with the relevant clinical covariates) is passed to
   a pre-fitted Cox proportional-hazards model to produce the predicted survival
   probability at 5 and/or 10 years.

---

## Repository layout

```
.
├── MB_appJun26_4.R                 # Shiny app entry point (UI + server)
├── MB_app_source_functionsJun26.R  # Package loading, processing & plotting functions
├── create_cox_models.r             # Builds the Cox OS models from phenotype CSVs
├── ExampleDataScript.R             # Headless, scriptable end-to-end example
├── introduction.md                 # Intro tab content
├── Tutorial/tutorial.md            # Tutorial tab content
├── www/MB_risk.jpg                 # Sidebar logo
└── mb_app/
    ├── SHHn49.model.rds            # SHH signature: CpGs, weights, scaling
    ├── G3_G4_subgroup.model.rds    # G3/4 (subgroup-aware) signature
    ├── G3_G4_no_subgroup.model.rds # G3/4 (combined) signature
    ├── *.dist.rds                  # Reference-cohort score distributions
    ├── shh49.train.fit.rds         # Cox model: SHH 5-year OS
    ├── g34early19.train.fit.rds    # Cox model: G3/4 early 5-year OS
    └── g34late49.train.fit.rds     # Cox model: G3/4 late 10-year mortality
```

> Some `.rds` files are large and may be tracked with
> [Git LFS](https://git-lfs.com/). After cloning, run `git lfs pull` if needed.

---

## Installation

Requires **R ≥ 4.2** and Bioconductor.

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "minfi", "minfiData", "sva", "bumphunter",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
))

# CRAN packages
install.packages(c(
  "shiny", "bslib", "shinyWidgets", "shinydashboard", "DT",
  "ggplot2", "ggpubr", "ggnewscale", "ggrepel", "pheatmap",
  "shinycssloaders", "shinybusy", "waiter",
  "foreach", "survival", "gridExtra", "fs", "future", "promises", "tinytex"
))
```

The source file will attempt to install any missing packages on first run, but a
clean install up front is faster and more reliable.

---

## Quick start (Shiny app)

```r
shiny::runApp("MB_appJun26_4.R")
```

Then in the browser:

1. **Upload IDATs.** Select the IDAT files for your sample(s). Each sample needs
   **both** the `_Grn.idat` and `_Red.idat` files; filenames must share the same
   sentrix basename (e.g. `200xxxx_R01C01_Grn.idat` and `..._Red.idat`).
2. **Select a model** (SHH / G3_G4_sub / G3_G4_no_sub).
3. **Run analysis.** Processing runs in the background; progress is shown on the
   bar. The Results and Download tabs appear when finished.
4. **Review** the risk table, distribution plot, and survival readout.
5. **Download** results as CSV or PDF.

---

## Quick start (headless / scripted)

For batch use or reproducible pipelines, adapt `ExampleDataScript.R`:

```r
source("./MB_app_source_functionsJun26.R")

# Folder containing the IDAT pairs
idats     <- "/path/to/idat/folder"
temp.base <- get_basenames(idats)
processed <- process_idats(temp.base)

# Choose a signature
meta <- SHH        # or G3_G4_sub / G3_G4_no_sub

# Compute the metaCpG score
res <- extract.metagene(
  as.character(meta[[1]]$genes),
  as.numeric(meta[[1]]$weights),
  beta2m(processed$betas),
  as.numeric(meta[[2]])
)
round(res, 3)
```

A built-in sanity check uses the `minfiData` example arrays; expected values for
the legacy reference signatures are documented inline in `ExampleDataScript.R`.

---

## Refitting the survival models

`create_cox_models.r` rebuilds the Cox OS models from phenotype tables. Each CSV
must contain `OS_Time`, `OS_Status`, the metagene score column (`fixedMG`), and
the relevant covariates (`M._versus_M.`, `ConsensusMYC`, `ConsensusMYCN`). This
is only needed if you are retraining; the pre-fitted `*.train.fit.rds` files ship
with the repository.

---

## Input requirements

- **Array types:** Illumina Infinium 450K or EPIC (v1). Mixed batches are
  combined onto the 450K overlap automatically.
- **File pairs:** every sample requires both `_Grn.idat` and `_Red.idat`.
- **Coverage:** all signature CpGs must be present on the array. EPIC and 450K
  both cover the signatures used here; very sparse or heavily failed arrays
  (high detection-p) may give unreliable scores.

---

## Intended use and limitations

- **Research use only.** Not validated as a diagnostic or prognostic medical
  device and not a substitute for clinical judgement, histopathology, or
  established risk stratification.
- **Population-level estimates.** Predicted OS is a model output describing the
  average outcome for cases with a similar profile in the training cohorts — not
  a deterministic prediction for an individual.
- **Cohort dependence.** Scores are z-scaled to, and percentiles are computed
  against, the discovery cohort. Estimates assume the input sample is drawn from
  a comparable population (subgroup, treatment era, assay).
- **Assay sensitivity.** Pre-analytic factors (low input, degraded DNA, batch
  effects, failed probes) will affect the methylation profile and therefore the score.
- **No subgroup classification.** The tool assumes the molecular context (SHH vs
  Group 3/4) is already known; it does not assign subgroup.

---

## Citation

If you use this tool, please cite:

> *[Schwalbe et al., 2026]*

---

## Acknowledgements & provenance

This application is repurposed from an earlier methylation risk-scoring tool
originally developed for malignant rhabdoid tumours (MRT/ATRT/ECRT). The
processing, scoring, and plotting framework is shared; the signatures, reference
distributions, and survival models here are specific to medulloblastoma.

## License

*[LICENSE — ]*

## Contact

*[CONTACT — Dan Williamson daniel.willamson@ncl.ac.uk 
Ed Schwalbe ed.schwalbe@northumrbia.ac.uk (corresponding authors) / 
James Hacking james.hacking@ncl.ac.uk (repository maintainer)]*
