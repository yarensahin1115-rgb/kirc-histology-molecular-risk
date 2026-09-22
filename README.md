# Histology-predicted molecular risk in ccRCC

Code, fitted models, and analysis utilities for predicting a prognostic
25-gene molecular risk signature from H&E whole-slide images in clear-cell
renal cell carcinoma (ccRCC).

The primary development and held-out validation analyses use TCGA-KIRC.
The locked molecular expression signature is additionally evaluated in the
independent CPTAC-CCRCC cohort.

## Contents

### Analysis

- `analysis.R` — prognostic signature development; held-out survival analysis;
  patient-level bootstrap confidence intervals; time-dependent AUC analysis;
  clinical and image-plus-clinical Cox models; nuclear-grade analyses;
  grade/stage-adjusted image–molecular correlation; repeated-split sensitivity
  analysis; slide-preparation batch-effect assessment; and direct
  histology-to-survival benchmarking

- `extract_features.py` — WSI tissue detection, patch sampling, Phikon feature
  extraction, and slide-level feature aggregation

- `extract_features2.py` — additional/revised feature-extraction utility
  used during pipeline development

### Fitted models and coefficients

- `imza_katsayilar.csv` — coefficients of the 25-gene LASSO-Cox prognostic
  signature

- `imza_up_down.csv` — signature genes with tumor-versus-normal log2 fold
  change, univariable survival hazard ratio, and prognostic direction

- `ridge_coef.csv` — coefficients of the histology-to-molecular-risk ridge model

- `signature_model.rds` — fitted LASSO-Cox model and associated signature genes

### Figure-generation scripts

- `fig_updown4.py` — visualization of signature-gene direction and associated
  molecular characteristics

- `fig_gene_morph.py` — exploratory morphology montages for selected
  image-predictable signature genes

## Study design

The matched TCGA-KIRC cohort comprised 528 patients and was divided at the
patient level into 369 development and 159 held-out test cases.

A 25-gene prognostic expression signature was derived in the development set
using variance filtering, univariable Cox-based gene ranking, and
LASSO-penalized Cox regression.

H&E whole-slide images were processed at approximately 20× magnification.
Up to 200 eligible 224 × 224-pixel tissue patches per slide were sampled and
encoded using the frozen Phikon histopathology foundation model. Patch
embeddings were mean pooled into 768-dimensional slide representations.

A ridge-regression model was trained in the development set to predict the
continuous expression-derived molecular risk score from histology.

The held-out test set was used for evaluation of image–molecular concordance,
survival discrimination, clinicopathologic adjustment, time-dependent AUC,
and sensitivity analyses.

A direct ridge-Cox histology-to-survival model was also evaluated as a
benchmark against the molecular-intermediate approach.

The locked 25-gene molecular signature was externally evaluated in the
independent CPTAC-CCRCC cohort without gene reselection, coefficient refitting,
or outcome-based recalibration.

## Reproducibility notes

- Patient-level development/test separation was used throughout supervised
  modeling.
- Repeated-split analyses reran gene screening, LASSO fitting,
  image-model fitting, clinical modeling, and held-out evaluation independently
  across 10 event-stratified patient-level splits.
- Penalized Cox models were fitted with `glmnet`.
- Cox tie handling for the direct histology-to-survival ridge model was fixed
  to Breslow for reproducibility.
- The direct histology-to-survival benchmark used 10-fold cross-validation and
  `lambda.min` for penalty selection.
- Principal held-out C-index confidence intervals were estimated using 2,000
  patient-level percentile bootstrap resamples.
- Time-dependent AUCs were evaluated at 1, 3, and 5 years.
- No explicit stain-normalization algorithm was applied before Phikon feature
  extraction.
- Slide-preparation sensitivity analyses evaluated FFPE and frozen tumor-slide
  composition.

## Data availability

Source TCGA-KIRC and CPTAC data are not redistributed in this repository.

The source transcriptomic, histologic, clinicopathologic, and survival data are
publicly available through the NCI Genomic Data Commons:

- TCGA-KIRC:
  https://portal.gdc.cancer.gov/projects/TCGA-KIRC

- CPTAC:
  https://portal.gdc.cancer.gov/projects/CPTAC-3

Patient-level source data are not included in this repository.

## Requirements

### R

- R 4.4.2
- DESeq2 1.52.0
- glmnet 5.0
- survival 3.8-6
- survminer 0.5.2
- TCGAbiolinks 2.40.0
- uwot 0.2.4
- pheatmap 1.0.13
- patchwork 1.3.2
- ggplot2 4.0.3
- EnhancedVolcano 1.30.0
- timeROC

### Python

- Python 3.12.10
- PyTorch 2.5.1+cu121
- Transformers 5.12.1
- openslide-python 1.4.6
- NumPy 2.4.4
- Pillow 12.2.0
- Matplotlib 3.11.0

## Model and annotation

Feature encoder: Phikon (`owkin/phikon`), used as a frozen histopathology
foundation model without task-specific fine-tuning.

Gene annotation: GENCODE v36 / GRCh38, consistent with the GDC STAR-Counts
pipeline.

## Citation

If you use this repository, please cite the associated manuscript once
published.
