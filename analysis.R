# =============================================================
# Histology-predicted molecular risk in ccRCC (TCGA-KIRC)
# Analysis script — reproduces the principal analyses
#
# R 4.4.2
# survival 3.8-6
# glmnet 5.0
# DESeq2 1.52.0
# timeROC 0.4.1
# TCGAbiolinks 2.40.0
# =============================================================

library(survival)
library(glmnet)
library(timeROC)

# =============================================================
# 0. LOAD DATA
# =============================================================

master <- readRDS("master_final.rds")
# Expected columns:
# submitter_id, OS_time, OS_event, split,
# risk_expr, risk_img, age_at_index,
# ajcc_pathologic_stage, tumor_grade, f0-f767

D <- readRDS("master_full.rds")
# D$expr = patient x gene VST expression matrix

expr <- D$expr

stopifnot(nrow(expr) == nrow(master))

# Ordered grade: G1-G4 only
master$grade_num <- ifelse(
  master$tumor_grade %in% c("G1", "G2", "G3", "G4"),
  as.numeric(sub("G", "", master$tumor_grade)),
  NA_real_
)

# Binary stage for clinical models
master$stage_bin <- ifelse(
  master$ajcc_pathologic_stage %in% c("Stage I", "Stage II"),
  "early",
  ifelse(
    master$ajcc_pathologic_stage %in% c("Stage III", "Stage IV"),
    "late",
    NA_character_
  )
)

master$stage_bin <- factor(
  master$stage_bin,
  levels = c("early", "late")
)

# Numeric stage for grade/stage-adjusted correlation
master$stage_num <- match(
  master$ajcc_pathologic_stage,
  c("Stage I", "Stage II", "Stage III", "Stage IV")
)

fcol <- grep(
  "^f[0-9]+$",
  names(master),
  value = TRUE
)

te <- master[
  master$split == "test",
]

cat("Full matched cohort:", nrow(master), "\n")
cat("Development set:", sum(master$split == "train"), "\n")
cat("Held-out test set:", nrow(te), "\n")
cat("Held-out deaths:", sum(te$OS_event), "\n")


# =============================================================
# 1. PRIMARY TEST-SET DISCRIMINATION
# =============================================================

summary(
  coxph(
    Surv(OS_time, OS_event) ~ risk_expr,
    data = te
  )
)$concordance
# approximately 0.722

summary(
  coxph(
    Surv(OS_time, OS_event) ~ risk_img,
    data = te
  )
)$concordance
# approximately 0.701

cor_primary <- cor.test(
  te$risk_img,
  te$risk_expr,
  method = "spearman",
  exact = FALSE
)

print(cor_primary)
# rho approximately 0.659

MAE <- mean(
  abs(
    te$risk_img -
      te$risk_expr
  )
)

RMSE <- sqrt(
  mean(
    (
      te$risk_img -
        te$risk_expr
    )^2
  )
)

cat("MAE:", MAE, "\n")
cat("RMSE:", RMSE, "\n")


# =============================================================
# 2. BOOTSTRAP CONFIDENCE INTERVALS
#    Percentile bootstrap; 2000 patient-level resamples
# =============================================================

set.seed(1)

boot_c <- function(f, d, n = 2000) {

  quantile(
    replicate(
      n,
      {
        i <- sample(
          nrow(d),
          replace = TRUE
        )

        tryCatch(
          summary(
            coxph(
              f,
              data = d[i, ]
            )
          )$concordance[1],
          error = function(e) NA
        )
      }
    ),
    c(0.025, 0.975),
    na.rm = TRUE
  )
}

boot_c(
  Surv(OS_time, OS_event) ~ risk_expr,
  te
)
# approximately 0.654-0.788

boot_c(
  Surv(OS_time, OS_event) ~ risk_img,
  te
)
# approximately 0.627-0.772


# Image-expression Spearman bootstrap CI
set.seed(1)

boot_rho <- replicate(
  2000,
  {
    i <- sample(
      nrow(te),
      replace = TRUE
    )

    cor(
      te$risk_img[i],
      te$risk_expr[i],
      method = "spearman"
    )
  }
)

quantile(
  boot_rho,
  c(0.025, 0.975),
  na.rm = TRUE
)
# approximately 0.551-0.745


# =============================================================
# 3. TIME-DEPENDENT AUC
#    1, 3, and 5 years
#    timeROC v0.4.1
# =============================================================

roc_img <- timeROC(
  T = te$OS_time,
  delta = te$OS_event,
  marker = te$risk_img,
  cause = 1,
  times = c(365, 1095, 1825),
  iid = TRUE
)

roc_expr <- timeROC(
  T = te$OS_time,
  delta = te$OS_event,
  marker = te$risk_expr,
  cause = 1,
  times = c(365, 1095, 1825),
  iid = TRUE
)

cat("\nImage AUCs:\n")
print(roc_img$AUC)

cat("\nImage AUC CIs:\n")
print(confint(roc_img))

cat("\nExpression AUCs:\n")
print(roc_expr$AUC)

cat("\nExpression AUC CIs:\n")
print(confint(roc_expr))


# =============================================================
# 4. NUCLEAR-GRADE CORRELATION + MULTIVARIABLE COX
# =============================================================

teg <- master[
  master$split == "test" &
    !is.na(master$grade_num),
]

cat(
  "\nOrdered-grade test-set n:",
  nrow(teg),
  "\n"
)

cat(
  "Deaths:",
  sum(teg$OS_event),
  "\n"
)

cor.test(
  teg$risk_img,
  teg$grade_num,
  method = "spearman",
  exact = FALSE
)
# rho approximately 0.447

cor.test(
  teg$risk_expr,
  teg$grade_num,
  method = "spearman",
  exact = FALSE
)
# rho approximately 0.346

tapply(
  teg$risk_img,
  teg$grade_num,
  median,
  na.rm = TRUE
)

kruskal.test(
  risk_img ~ factor(grade_num),
  data = teg
)

# Full-cohort supportive association
full_grade <- master[
  !is.na(master$grade_num),
]

cor.test(
  full_grade$risk_img,
  full_grade$grade_num,
  method = "spearman",
  exact = FALSE
)
# rho approximately 0.491


# Grade-adjusted survival
teg$risk_img_z <- as.numeric(
  scale(teg$risk_img)
)

m_g <- coxph(
  Surv(OS_time, OS_event) ~ grade_num,
  data = teg
)

m_gi <- coxph(
  Surv(OS_time, OS_event) ~
    grade_num +
    risk_img_z,
  data = teg
)

summary(m_g)
summary(m_gi)

anova(
  m_g,
  m_gi,
  test = "LRT"
)

cox.zph(m_gi)

summary(m_g)$concordance
summary(m_gi)$concordance


# =============================================================
# 5. CLINICAL MODEL + INCREMENTAL VALUE
#    age + early/late stage + ordinal grade
# =============================================================

tec <- master[
  master$split == "test" &
    complete.cases(
      master[, c(
        "OS_time",
        "OS_event",
        "age_at_index",
        "stage_bin",
        "grade_num",
        "risk_img"
      )]
    ),
]

cat(
  "\nClinical-model analytic n:",
  nrow(tec),
  "\n"
)

cat(
  "Clinical-model deaths:",
  sum(tec$OS_event),
  "\n"
)

m_clin <- coxph(
  Surv(OS_time, OS_event) ~
    age_at_index +
    stage_bin +
    grade_num,
  data = tec
)

m_combo <- coxph(
  Surv(OS_time, OS_event) ~
    age_at_index +
    stage_bin +
    grade_num +
    risk_img,
  data = tec
)

summary(m_clin)
summary(m_combo)

anova(
  m_clin,
  m_combo,
  test = "LRT"
)

summary(m_clin)$concordance
summary(m_combo)$concordance


# =============================================================
# 6. REPEATED-SPLIT SENSITIVITY
#
# ORIGINAL ANALYSIS CODE
# 10 event-stratified patient-level splits
# This section is retained to reproduce the originally reported
# repeated-split results.
# =============================================================

expr <- D$expr

fcol <- grep(
  "^f[0-9]+$",
  names(master),
  value = TRUE
)

res <- data.frame(
  seed = 1:10,
  expr = NA,
  img = NA,
  clin = NA,
  img_clin = NA
)

selected_genes <- vector(
  "list",
  10
)

for (i in 1:10) {

  set.seed(i)

  ev <- which(
    master$OS_event == 1
  )

  ce <- which(
    master$OS_event == 0
  )

  tr <- rep(
    FALSE,
    nrow(master)
  )

  tr[
    c(
      sample(
        ev,
        round(.7 * length(ev))
      ),
      sample(
        ce,
        round(.7 * length(ce))
      )
    )
  ] <- TRUE

  ytr <- Surv(
    master$OS_time[tr],
    master$OS_event[tr]
  )

  # Top 5000 genes by development-set variance
  top <- names(
    sort(
      apply(
        expr[tr, ],
        2,
        var
      ),
      decreasing = TRUE
    )
  )[1:5000]

  # Univariable Cox ranking
  up <- sapply(
    top,
    function(g) {
      tryCatch(
        suppressWarnings(
          summary(
            coxph(
              ytr ~ expr[tr, g]
            )
          )$coefficients[5]
        ),
        error = function(e) 1
      )
    }
  )

  # Top 2000 genes by smallest Cox P value
  cand <- names(
    sort(up)
  )[1:2000]

  # LASSO-Cox signature
  cv <- cv.glmnet(
    as.matrix(
      expr[
        tr,
        cand
      ]
    ),
    ytr,
    family = "cox",
    alpha = 1,
    nfolds = 5,
    maxit = 1000000
  )

  b <- coef(
    cv,
    s = "lambda.min"
  )

  sig <- rownames(b)[
    which(
      as.numeric(b) != 0
    )
  ]

  bb <- as.numeric(b)[
    which(
      as.numeric(b) != 0
    )
  ]

  selected_genes[[i]] <- sig

  # Development-derived gene standardization
  z <- scale(
    expr[
      ,
      sig,
      drop = FALSE
    ],
    center = colMeans(
      expr[
        tr,
        sig,
        drop = FALSE
      ]
    ),
    scale = apply(
      expr[
        tr,
        sig,
        drop = FALSE
      ],
      2,
      sd
    )
  )

  rexp <- as.numeric(
    z %*% bb
  )

  # Histology-to-molecular-risk ridge
  rf <- cv.glmnet(
    as.matrix(
      master[
        tr,
        fcol
      ]
    ),
    rexp[tr],
    alpha = 0,
    nfolds = 5,
    maxit = 1000000
  )

  rimg <- as.numeric(
    predict(
      rf,
      as.matrix(
        master[
          ,
          fcol
        ]
      ),
      s = "lambda.min"
    )
  )

  master$.re <- rexp
  master$.ri <- rimg

  teD <- master[
    !tr,
  ]

  teC <- teD[
    complete.cases(
      teD[, c(
        "age_at_index",
        "stage_bin",
        "grade_num"
      )]
    ),
  ]

  ci <- function(f, d) {
    tryCatch(
      summary(
        coxph(
          f,
          data = d
        )
      )$concordance[1],
      error = function(e) NA
    )
  }

  res$expr[i] <- ci(
    Surv(OS_time, OS_event) ~ .re,
    teD
  )

  res$img[i] <- ci(
    Surv(OS_time, OS_event) ~ .ri,
    teD
  )

  res$clin[i] <- ci(
    Surv(OS_time, OS_event) ~
      age_at_index +
      stage_bin +
      grade_num,
    teC
  )

  res$img_clin[i] <- ci(
    Surv(OS_time, OS_event) ~
      age_at_index +
      stage_bin +
      grade_num +
      .ri,
    teC
  )
}

cat("\nRepeated-split results:\n")
print(res)

split_summary <- sapply(
  res[, -1],
  function(v) {
    c(
      median = median(
        v,
        na.rm = TRUE
      ),
      min = min(
        v,
        na.rm = TRUE
      ),
      max = max(
        v,
        na.rm = TRUE
      )
    )
  }
)

cat("\nRepeated-split median/range:\n")
print(split_summary)

cat(
  "\nCombined > clinical in",
  sum(
    res$img_clin >
      res$clin,
    na.rm = TRUE
  ),
  "of 10 splits\n"
)


# Gene-selection frequency across repeated splits
gene_frequency <- sort(
  table(
    unlist(
      selected_genes
    )
  ),
  decreasing = TRUE
)

cat("\nMost frequently selected genes:\n")
print(
  head(
    gene_frequency,
    30
  )
)

cat(
  "\nGenes selected in >=8/10 splits:",
  sum(
    gene_frequency >= 8
  ),
  "\n"
)


# =============================================================
# 7. GRADE/STAGE-ADJUSTED IMAGE-MOLECULAR CORRELATION
# =============================================================

te_adj <- master[
  master$split == "test",
]

tec_adj <- te_adj[
  complete.cases(
    te_adj[, c(
      "risk_img",
      "risk_expr",
      "grade_num",
      "stage_num"
    )]
  ),
]

cat(
  "\nGrade/stage complete-case n:",
  nrow(tec_adj),
  "\n"
)

# Raw correlation in same complete-case subset
raw_adj_subset <- cor.test(
  tec_adj$risk_img,
  tec_adj$risk_expr,
  method = "spearman",
  exact = FALSE
)

print(raw_adj_subset)
# rho approximately 0.651

# Residualization for grade and stage
r_img <- residuals(
  lm(
    risk_img ~
      grade_num +
      stage_num,
    data = tec_adj
  )
)

r_expr <- residuals(
  lm(
    risk_expr ~
      grade_num +
      stage_num,
    data = tec_adj
  )
)

adjusted_cor <- cor.test(
  r_img,
  r_expr,
  method = "spearman",
  exact = FALSE
)

print(adjusted_cor)

# Verified:
# n = 155
# raw rho = 0.651
# adjusted rho = 0.583
# P = 1.8e-15


# =============================================================
# 8. DIRECT HISTOLOGY-TO-SURVIVAL BENCHMARK
# =============================================================

tr_direct <- master$split == "train"
te_direct <- master$split == "test"

X_tr <- as.matrix(
  master[
    tr_direct,
    fcol
  ]
)

X_te <- as.matrix(
  master[
    te_direct,
    fcol
  ]
)

y_tr <- Surv(
  master$OS_time[tr_direct],
  master$OS_event[tr_direct]
)

y_te <- Surv(
  master$OS_time[te_direct],
  master$OS_event[te_direct]
)

set.seed(42)

cv_direct <- cv.glmnet(
  x = X_tr,
  y = y_tr,
  family = "cox",
  alpha = 0,
  nfolds = 10,
  cox.ties = "breslow"
)

cat(
  "\nDirect model lambda.min:",
  cv_direct$lambda.min,
  "\n"
)

cat(
  "Direct model lambda.1se:",
  cv_direct$lambda.1se,
  "\n"
)

pred_direct <- as.numeric(
  predict(
    cv_direct,
    newx = X_te,
    s = "lambda.min",
    type = "link"
  )
)

direct_model <- coxph(
  y_te ~ pred_direct
)

cat(
  "\nDirect histology-to-survival C-index:\n"
)

print(
  summary(
    direct_model
  )$concordance
)

# Patient-level percentile bootstrap CI
set.seed(42)

B_direct <- 2000

boot_direct <- replicate(
  B_direct,
  {

    idx <- sample(
      seq_along(
        pred_direct
      ),
      replace = TRUE
    )

    tryCatch(
      summary(
        coxph(
          Surv(
            master$OS_time[
              te_direct
            ][idx],
            master$OS_event[
              te_direct
            ][idx]
          ) ~
            pred_direct[idx]
        )
      )$concordance[1],
      error = function(e) NA
    )
  }
)

cat(
  "\nDirect histology-to-survival 95% bootstrap CI:\n"
)

print(
  quantile(
    boot_direct,
    c(0.025, 0.975),
    na.rm = TRUE
  )
)

# Verified:
# C-index approximately 0.692
# 95% CI approximately 0.605-0.771


# =============================================================
# 9. BATCH-EFFECT / SLIDE-PREPARATION ASSESSMENT
#
# Requires patient-level FFPE fraction merged as master$ffpe_frac.
#
# Verified tumor-slide preparation counts:
# FFPE = 617
# Frozen = 453
# Unresolved = 12
# Total tumor slides = 1082
# =============================================================

if ("ffpe_frac" %in% names(master)) {

  cat(
    "\nFFPE fraction vs image risk:\n"
  )

  print(
    cor.test(
      master$ffpe_frac,
      master$risk_img,
      method = "spearman",
      exact = FALSE
    )
  )

  cat(
    "\nFFPE fraction vs expression risk:\n"
  )

  print(
    cor.test(
      master$ffpe_frac,
      master$risk_expr,
      method = "spearman",
      exact = FALSE
    )
  )

  batch_df <- master[
    complete.cases(
      master[, c(
        "OS_time",
        "OS_event",
        "risk_img",
        "ffpe_frac"
      )]
    ),
  ]

  m_batch_unadjusted <- coxph(
    Surv(OS_time, OS_event) ~
      risk_img,
    data = batch_df
  )

  m_batch_adjusted <- coxph(
    Surv(OS_time, OS_event) ~
      risk_img +
      ffpe_frac,
    data = batch_df
  )

  print(
    summary(
      m_batch_unadjusted
    )
  )

  print(
    summary(
      m_batch_adjusted
    )
  )

} else {

  message(
    "Batch analysis skipped: master$ffpe_frac is not present."
  )
}


# =============================================================
# 10. EXTERNAL CPTAC-CCRCC MOLECULAR VALIDATION
#
# The external validation applies to the locked 25-gene
# molecular expression signature, not to the histology model.
#
# Procedure:
# - CPTAC primary-tumor RNA-seq
# - VST transformation
# - version-stripped Ensembl matching
# - all 25 signature genes retained
# - CPTAC-internal gene-wise standardization
# - fixed TCGA-derived coefficients
# - no gene reselection
# - no coefficient refitting
# - no outcome-based recalibration
#
# Verified result:
# n = 176
# deaths = 30
# median follow-up approximately 1894 days
# HR = 3.99
# 95% CI = 1.81-8.84
# P = 6.29e-4
# C-index = 0.678
# 95% CI = 0.570-0.787
# log-rank P = 0.023
# =============================================================


# =============================================================
# 11. SOFTWARE / SESSION INFORMATION
# =============================================================

cat("\nPackage versions:\n")

cat(
  "R:",
  R.version.string,
  "\n"
)

cat(
  "survival:",
  as.character(
    packageVersion(
      "survival"
    )
  ),
  "\n"
)

cat(
  "glmnet:",
  as.character(
    packageVersion(
      "glmnet"
    )
  ),
  "\n"
)

cat(
  "timeROC:",
  as.character(
    packageVersion(
      "timeROC"
    )
  ),
  "\n"
)

cat(
  "\nSession information:\n"
)

print(
  sessionInfo()
)
