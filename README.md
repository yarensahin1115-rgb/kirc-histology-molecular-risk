# Histology-predicted molecular risk in ccRCC

Code and fitted coefficients for predicting a 25-gene molecular prognostic 
signature from H&E whole-slide images in clear-cell renal cell carcinoma 
(TCGA-KIRC), with external validation in CPTAC-CCRCC.

## Contents

**Analysis**
- `analysis.R` — signature development, survival analysis, bootstrap CIs, 
  time-dependent AUC, nuclear-grade multivariable Cox, clinical-model 
  comparison, repeated-split sensitivity, batch-effect assessment
- `extract_features.py` — WSI tissue detection, patch sampling, and Phikon 
  feature extraction

**Fitted coefficients**
- `imza_katsayilar.csv` — 25-gene signature coefficients (LASSO Cox)
- `imza_up_down.csv` — signature genes with tumor–normal log2FC, hazard 
  ratios, and prognostic direction
- `ridge_coef.csv` — image-to-signature ridge model coefficients
- `signature_model.rds` — fitted LASSO Cox model and gene list

**Figures**
- `fig_updown4.py`, `fig_gene_morph.py` — morphology montages

## Data availability

Source TCGA and CPTAC data are not redistributed. They are publicly available 
from the NCI Genomic Data Commons (https://portal.gdc.cancer.gov). Patient-level 
data are not included in this repository.

## Requirements

Python 3.12 (PyTorch 2.5.1, transformers, openslide-python); 
R 4.4.2 (DESeq2, glmnet, survival, timeROC, TCGAbiolinks).

## Model

Feature encoder: Phikon (owkin/phikon), a ViT-B/16 histopathology foundation 
model. Gene annotation: GENCODE v36 (GDC STAR-Counts pipeline).
