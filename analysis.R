
# =============================================================
# Histology-predicted molecular risk in ccRCC (TCGA-KIRC)
# Analysis script — reproduces the principal analyses
# R 4.4.2; packages: survival, glmnet, DESeq2, timeROC, TCGAbiolinks
# =============================================================

library(survival); library(glmnet)

# ---- Load ----
master <- readRDS("master_final.rds")   # submitter_id, OS_time, OS_event, split,
                                         # risk_expr, risk_img, age_at_index,
                                         # ajcc_pathologic_stage, tumor_grade, f0-f767
D      <- readRDS("master_full.rds")     # $expr = patient x gene (VST); $master
master$grade_num <- as.numeric(gsub("G","", master$tumor_grade))
master$stage_bin <- ifelse(master$ajcc_pathologic_stage %in% c("Stage I","Stage II"), "early",
                    ifelse(master$ajcc_pathologic_stage %in% c("Stage III","Stage IV"), "late", NA))

# =============================================================
# 1. PRIMARY TEST-SET DISCRIMINATION
# =============================================================
te <- master[master$split=="test", ]
summary(coxph(Surv(OS_time,OS_event) ~ risk_expr, data=te))$concordance   # 0.722
summary(coxph(Surv(OS_time,OS_event) ~ risk_img,  data=te))$concordance   # 0.701
cor.test(te$risk_img, te$risk_expr, method="spearman")                    # r=0.659
mean(abs(te$risk_img - te$risk_expr))                                     # MAE
sqrt(mean((te$risk_img - te$risk_expr)^2))                                # RMSE

# =============================================================
# 2. BOOTSTRAP CONFIDENCE INTERVALS (percentile, 2000 resamples)
# =============================================================
set.seed(1)
boot_c <- function(f, d, n=2000) quantile(replicate(n, {
  i <- sample(nrow(d), replace=TRUE)
  tryCatch(summary(coxph(f, data=d[i,]))$concordance[1], error=function(e) NA)
}), c(0.025,0.975), na.rm=TRUE)
boot_c(Surv(OS_time,OS_event)~risk_expr, te)   # 0.654-0.788
boot_c(Surv(OS_time,OS_event)~risk_img,  te)   # 0.627-0.772

# =============================================================
# 3. TIME-DEPENDENT AUC (1,3,5 years)
# =============================================================
library(timeROC)
timeROC(T=te$OS_time, delta=te$OS_event, marker=te$risk_img,  cause=1, times=c(365,1095,1825), iid=TRUE)
timeROC(T=te$OS_time, delta=te$OS_event, marker=te$risk_expr, cause=1, times=c(365,1095,1825), iid=TRUE)

# =============================================================
# 4. NUCLEAR-GRADE CORRELATION + MULTIVARIABLE COX
# =============================================================
teg <- master[master$split=="test" & !is.na(master$grade_num), ]
cor.test(teg$risk_img, teg$grade_num, method="spearman")                  # rho=0.447
teg$risk_img_z <- as.numeric(scale(teg$risk_img))
m_g  <- coxph(Surv(OS_time,OS_event) ~ grade_num, data=teg)
m_gi <- coxph(Surv(OS_time,OS_event) ~ grade_num + risk_img_z, data=teg)
summary(m_gi)$conf.int                                                    # image HR 1.48 (1.09-2.01)
anova(m_g, m_gi)


analysis_code <- r"(
# =============================================================
# Histology-predicted molecular risk in ccRCC (TCGA-KIRC)
# Analysis script — reproduces the principal analyses
# R 4.4.2; packages: survival, glmnet, DESeq2, timeROC, TCGAbiolinks
# =============================================================

library(survival); library(glmnet)

# ---- Load ----
master <- readRDS("master_final.rds")   # submitter_id, OS_time, OS_event, split,
                                         # risk_expr, risk_img, age_at_index,
                                         # ajcc_pathologic_stage, tumor_grade, f0-f767
D      <- readRDS("master_full.rds")     # $expr = patient x gene (VST); $master
master$grade_num <- as.numeric(gsub("G","", master$tumor_grade))
master$stage_bin <- ifelse(master$ajcc_pathologic_stage %in% c("Stage I","Stage II"), "early",
                    ifelse(master$ajcc_pathologic_stage %in% c("Stage III","Stage IV"), "late", NA))

# =============================================================
# 1. PRIMARY TEST-SET DISCRIMINATION
# =============================================================
te <- master[master$split=="test", ]
summary(coxph(Surv(OS_time,OS_event) ~ risk_expr, data=te))$concordance   # 0.722
summary(coxph(Surv(OS_time,OS_event) ~ risk_img,  data=te))$concordance   # 0.701
cor.test(te$risk_img, te$risk_expr, method="spearman")                    # r=0.659
mean(abs(te$risk_img - te$risk_expr))                                     # MAE
sqrt(mean((te$risk_img - te$risk_expr)^2))                                # RMSE

# =============================================================
# 2. BOOTSTRAP CONFIDENCE INTERVALS (percentile, 2000 resamples)
# =============================================================
set.seed(1)
boot_c <- function(f, d, n=2000) quantile(replicate(n, {
  i <- sample(nrow(d), replace=TRUE)
  tryCatch(summary(coxph(f, data=d[i,]))$concordance[1], error=function(e) NA)
}), c(0.025,0.975), na.rm=TRUE)
boot_c(Surv(OS_time,OS_event)~risk_expr, te)   # 0.654-0.788
boot_c(Surv(OS_time,OS_event)~risk_img,  te)   # 0.627-0.772

# =============================================================
# 3. TIME-DEPENDENT AUC (1,3,5 years)
# =============================================================
library(timeROC)
timeROC(T=te$OS_time, delta=te$OS_event, marker=te$risk_img,  cause=1, times=c(365,1095,1825), iid=TRUE)
timeROC(T=te$OS_time, delta=te$OS_event, marker=te$risk_expr, cause=1, times=c(365,1095,1825), iid=TRUE)

# =============================================================
# 4. NUCLEAR-GRADE CORRELATION + MULTIVARIABLE COX
# =============================================================
teg <- master[master$split=="test" & !is.na(master$grade_num), ]
cor.test(teg$risk_img, teg$grade_num, method="spearman")                  # rho=0.447
teg$risk_img_z <- as.numeric(scale(teg$risk_img))
m_g  <- coxph(Surv(OS_time,OS_event) ~ grade_num, data=teg)
m_gi <- coxph(Surv(OS_time,OS_event) ~ grade_num + risk_img_z, data=teg)
summary(m_gi)$conf.int                                                    # image HR 1.48 (1.09-2.01)
anova(m_g, m_gi)
View(plot_df)
analysis_code <- r"(
# =============================================================
# Histology-predicted molecular risk in ccRCC (TCGA-KIRC)
# Analysis script — reproduces the principal analyses
# R 4.4.2; packages: survival, glmnet, DESeq2, timeROC, TCGAbiolinks
# =============================================================

library(survival); library(glmnet)

# ---- Load ----
master <- readRDS("master_final.rds")   # submitter_id, OS_time, OS_event, split,
                                         # risk_expr, risk_img, age_at_index,
                                         # ajcc_pathologic_stage, tumor_grade, f0-f767
D      <- readRDS("master_full.rds")     # $expr = patient x gene (VST); $master
master$grade_num <- as.numeric(gsub("G","", master$tumor_grade))
master$stage_bin <- ifelse(master$ajcc_pathologic_stage %in% c("Stage I","Stage II"), "early",
                    ifelse(master$ajcc_pathologic_stage %in% c("Stage III","Stage IV"), "late", NA))

# =============================================================
# 1. PRIMARY TEST-SET DISCRIMINATION
# =============================================================
te <- master[master$split=="test", ]
summary(coxph(Surv(OS_time,OS_event) ~ risk_expr, data=te))$concordance   # 0.722
summary(coxph(Surv(OS_time,OS_event) ~ risk_img,  data=te))$concordance   # 0.701
cor.test(te$risk_img, te$risk_expr, method="spearman")                    # r=0.659
mean(abs(te$risk_img - te$risk_expr))                                     # MAE
sqrt(mean((te$risk_img - te$risk_expr)^2))                                # RMSE

# =============================================================
# 2. BOOTSTRAP CONFIDENCE INTERVALS (percentile, 2000 resamples)
# =============================================================
set.seed(1)
boot_c <- function(f, d, n=2000) quantile(replicate(n, {
  i <- sample(nrow(d), replace=TRUE)
  tryCatch(summary(coxph(f, data=d[i,]))$concordance[1], error=function(e) NA)
}), c(0.025,0.975), na.rm=TRUE)
boot_c(Surv(OS_time,OS_event)~risk_expr, te)   # 0.654-0.788
boot_c(Surv(OS_time,OS_event)~risk_img,  te)   # 0.627-0.772

# =============================================================
# 3. TIME-DEPENDENT AUC (1,3,5 years)
# =============================================================
library(timeROC)
timeROC(T=te$OS_time, delta=te$OS_event, marker=te$risk_img,  cause=1, times=c(365,1095,1825), iid=TRUE)
timeROC(T=te$OS_time, delta=te$OS_event, marker=te$risk_expr, cause=1, times=c(365,1095,1825), iid=TRUE)

# =============================================================
# 4. NUCLEAR-GRADE CORRELATION + MULTIVARIABLE COX
# =============================================================
teg <- master[master$split=="test" & !is.na(master$grade_num), ]
cor.test(teg$risk_img, teg$grade_num, method="spearman")                  # rho=0.447
teg$risk_img_z <- as.numeric(scale(teg$risk_img))
m_g  <- coxph(Surv(OS_time,OS_event) ~ grade_num, data=teg)
m_gi <- coxph(Surv(OS_time,OS_event) ~ grade_num + risk_img_z, data=teg)
summary(m_gi)$conf.int                                                    # image HR 1.48 (1.09-2.01)
anova(m_g, m_gi)                                                          # LR p=0.011
cox.zph(m_gi)                                                             # PH global p=0.087

# =============================================================
# 5. CLINICAL MODEL (collapsed stage) + INCREMENTAL VALUE
#    covariates: age (continuous, per year),
#                grade (ordinal 1-4), stage (binary early/late)
# =============================================================
tec <- master[master$split=="test" &
              complete.cases(master[,c("OS_time","OS_event","age_at_index","stage_bin","grade_num","risk_img")]) &
              master$split=="test", ]
tec <- tec[tec$split=="test",]
m_clin  <- coxph(Surv(OS_time,OS_event) ~ age_at_index + stage_bin + grade_num, data=tec)
m_combo <- coxph(Surv(OS_time,OS_event) ~ age_at_index + stage_bin + grade_num + risk_img, data=tec)
summary(m_clin)$conf.int; summary(m_clin)$concordance                     # C=0.803
summary(m_combo)$conf.int; summary(m_combo)$concordance                   # C=0.804
anova(m_clin, m_combo)                                                    # LR p=0.040

# =============================================================
# 6. REPEATED-SPLIT SENSITIVITY (10 splits, leakage-safe:
#    all standardization/fitting from training fold only)
# =============================================================
expr <- D$expr; fcol <- grep("^f[0-9]+$", names(master), value=TRUE)
res <- data.frame(seed=1:10, expr=NA, img=NA, clin=NA, img_clin=NA)
for(i in 1:10){
  set.seed(i)
  ev <- which(master$OS_event==1); ce <- which(master$OS_event==0)
  tr <- rep(FALSE,nrow(master)); tr[c(sample(ev,round(.7*length(ev))),sample(ce,round(.7*length(ce))))] <- TRUE
  ytr <- Surv(master$OS_time[tr], master$OS_event[tr])
  top <- names(sort(apply(expr[tr,],2,var),decreasing=TRUE))[1:5000]
  up  <- sapply(top, function(g) tryCatch(summary(coxph(ytr~expr[tr,g]))$coefficients[5],error=function(e)1))
  cand<- names(sort(up))[1:2000]
  cv  <- cv.glmnet(as.matrix(expr[tr,cand]), ytr, family="cox", alpha=1, nfolds=5)
  b   <- coef(cv,s="lambda.min"); sig <- rownames(b)[which(as.numeric(b)!=0)]; bb <- as.numeric(b)[which(as.numeric(b)!=0)]
  z   <- scale(expr[,sig,drop=FALSE], center=colMeans(expr[tr,sig,drop=FALSE]), scale=apply(expr[tr,sig,drop=FALSE],2,sd))
  rexp<- as.numeric(z %*% bb)
  rf  <- cv.glmnet(as.matrix(master[tr,fcol]), rexp[tr], alpha=0, nfolds=5)
  rimg<- as.numeric(predict(rf, as.matrix(master[,fcol]), s="lambda.min"))
  master$.re <- rexp; master$.ri <- rimg; teD <- master[!tr,]
  teC <- teD[complete.cases(teD[,c("age_at_index","stage_bin","grade_num")]),]
  ci <- function(f,d) tryCatch(summary(coxph(f,data=d))$concordance[1],error=function(e)NA)
  res$expr[i]<-ci(Surv(OS_time,OS_event)~.re,teD); res$img[i]<-ci(Surv(OS_time,OS_event)~.ri,teD)
  res$clin[i]<-ci(Surv(OS_time,OS_event)~age_at_index+stage_bin+grade_num,teC)
  res$img_clin[i]<-ci(Surv(OS_time,OS_event)~age_at_index+stage_bin+grade_num+.ri,teC)
}
sapply(res[,-1], function(v) c(median=median(v),min=min(v),max=max(v)))

# =============================================================
# 7. BATCH-EFFECT ASSESSMENT (FFPE fraction from slide barcodes)
#    slide type: BS=FFPE, TS=frozen (from GDC file_name)
# =============================================================
# per-patient FFPE fraction merged into master$ffpe_frac (see repository notes)
# cor.test(ffpe_frac, risk_img); adjusted Cox unchanged HR; subgroup C-index

# =============================================================
# 8. EXTERNAL VALIDATION (CPTAC-CCRCC) — frozen-model transfer
#    fixed development coefficients applied to CPTAC (no refit)
# =============================================================
# See external_validation.R:
#   CPTAC-3 primary-tumor RNA-seq -> VST -> 25 genes by Ensembl ID
#   -> z-standardize -> risk = sum(beta * z) -> Cox / KM
#   Result: HR 3.99 (1.81-8.84), C-index 0.678 (0.570-0.787), log-rank p=0.023

