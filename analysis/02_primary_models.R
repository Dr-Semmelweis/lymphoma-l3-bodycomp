# -----------------------------------------------------------------------------
# 02_primary_models.R
#
# Co-primary Cox proportional-hazards models.
#
# Reproduces:
#   - Table 3: both co-primary models (overall and progression-free survival),
#     all covariates, with concordance and Schoenfeld tests;
#   - the bootstrap confidence intervals for the primary exposure;
#   - the restricted-cubic-spline tests for non-linearity;
#   - the comparison models without the continuous age term, quoted in the
#     Discussion;
#   - the events-per-variable ratios quoted in Methods.
#
# Input:  data/clinical_dataset_n155.csv
#         data/imaging_L3_quantitative_n155.csv
# Output: printed to standard output; no files are written.
#
# Run from the package root:  Rscript code/02_primary_models.R
# Requires R >= 4.5 and the packages 'survival', 'boot', 'rms'.
# The bootstrap (2 x 5,000 refits) takes one to two minutes.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(survival)
  library(boot)
  library(rms)
})

clin <- read.csv(file.path("data", "clinical_dataset_n155.csv"), na.strings = c("NA", ""))
img  <- read.csv(file.path("data", "imaging_L3_quantitative_n155.csv"), na.strings = c("NA", ""))
stopifnot(nrow(clin) == 155, identical(clin$patient_id, img$patient_id))

# Analysis frame. Reference levels are fixed explicitly: male sex and low IPI.
# The exposure and age are standardised to mean 0, SD 1 within the cohort, so
# that hazard ratios are expressed per one standard deviation.
z <- function(x) (x - mean(x)) / sd(x)
df <- data.frame(
  os_months  = clin$os_months,  os_event  = clin$os_event,
  pfs_months = clin$pfs_months, pfs_event = clin$pfs_event,
  smd_z = z(img$smd_mean_hu),
  age_z = z(clin$age_years),
  sex   = factor(clin$sex, levels = c("M", "F")),
  ipi   = factor(clin$ipi_group, levels = c("Low", "L-Int", "H-Int", "High"))
)
cat(sprintf("n = %d | OS events = %d | PFS events = %d\n", nrow(df),
            sum(df$os_event), sum(df$pfs_event)))
cat(sprintf("1 SD of muscle density = %.3f HU | 1 SD of age = %.1f years\n\n",
            sd(img$smd_mean_hu), sd(clin$age_years)))

# Prints every coefficient of a fitted model, then concordance and the
# Schoenfeld proportional-hazards tests (global, and for the exposure).
report <- function(m, label) {
  s <- summary(m)
  cat("====", label, "====\n")
  for (v in rownames(s$coefficients)) {
    ci <- s$conf.int[v, ]
    cat(sprintf("  %-12s HR %.3f (%.3f-%.3f)  p=%.4f\n",
                v, ci[1], ci[3], ci[4], s$coefficients[v, 5]))
  }
  cat(sprintf("  concordance %.3f | events %d | n %d\n",
              s$concordance[1], m$nevent, m$n))
  zph <- cox.zph(m)
  cat(sprintf("  Schoenfeld: global p=%.2f | muscle density p=%.2f\n\n",
              zph$table["GLOBAL", "p"], zph$table["smd_z", "p"]))
}

# ---- Co-primary models ------------------------------------------------------
# Efron's method for tied event times is the coxph() default.
os_fit  <- coxph(Surv(os_months,  os_event)  ~ smd_z + sex + ipi + age_z, data = df)
pfs_fit <- coxph(Surv(pfs_months, pfs_event) ~ smd_z + sex + ipi + age_z, data = df)
report(os_fit,  "Overall survival (co-primary)")
report(pfs_fit, "Progression-free survival (co-primary)")

# ---- Comparison without the continuous age term -----------------------------
# Quantifies what continuous age adjustment adds beyond the dichotomised age
# term already inside the IPI (Discussion: OS 0.42 -> 0.51; PFS 0.66 -> 0.71).
hr1 <- function(m) {
  ci <- summary(m)$conf.int["smd_z", ]
  sprintf("HR %.3f (%.3f-%.3f)", ci[1], ci[3], ci[4])
}
cat("==== Muscle density without the continuous age term ====\n")
cat("  OS : ", hr1(coxph(Surv(os_months,  os_event)  ~ smd_z + sex + ipi, data = df)), "\n")
cat("  PFS: ", hr1(coxph(Surv(pfs_months, pfs_event) ~ smd_z + sex + ipi, data = df)), "\n\n")

# ---- Bootstrap confidence intervals for the primary exposure ----------------
# Nonparametric check of the model-based interval: 5,000 resamples drawn within
# strata defined by the categorical covariates (sex x IPI group), hazard ratio
# for muscle density refitted on each, bias-corrected and accelerated interval.
# The seed reproduces the published intervals exactly; resamples in which the
# model fails to converge are discarded and counted.
set.seed(20260522)
strata <- interaction(df$sex, df$ipi, drop = TRUE)
boot_hr <- function(time, event) {
  f <- as.formula(sprintf("Surv(%s, %s) ~ smd_z + sex + ipi + age_z", time, event))
  function(data, i) {
    tryCatch(as.numeric(exp(coef(coxph(f, data = data[i, ]))["smd_z"])),
             error = function(e) NA_real_)
  }
}
cat("==== Bootstrap (R = 5,000, stratified, BCa) ====\n")
for (ep in list(c("os_months", "os_event", "OS"),
                c("pfs_months", "pfs_event", "PFS"))) {
  b <- boot(df, boot_hr(ep[1], ep[2]), R = 5000, strata = strata)
  failed <- sum(is.na(b$t))
  b$t <- matrix(b$t[!is.na(b$t)], ncol = 1)
  ci <- boot.ci(b, type = "bca")
  cat(sprintf("  %-4s HR %.3f  BCa 95%% CI %.3f-%.3f  (discarded resamples: %d/5000)\n",
              ep[3], b$t0, ci$bca[4], ci$bca[5], failed))
}

# ---- Restricted cubic splines: test for non-linearity -----------------------
# The exposure enters a four-knot spline on its original scale (HU); the other
# covariates are unchanged. The reported p-value tests the joint null that the
# non-linear spline components are zero.
df$smd_hu <- img$smd_mean_hu
dd <- datadist(df[, c("smd_hu", "sex", "ipi", "age_z",
                      "os_months", "os_event", "pfs_months", "pfs_event")])
options(datadist = "dd")
cat("\n==== Restricted cubic splines (4 knots) ====\n")
for (ep in list(c("os_months", "os_event", "OS"),
                c("pfs_months", "pfs_event", "PFS"))) {
  f <- cph(as.formula(sprintf("Surv(%s, %s) ~ rcs(smd_hu, 4) + sex + ipi + age_z",
                              ep[1], ep[2])), data = df, x = TRUE, y = TRUE)
  a <- anova(f)
  p_nl <- a[grep("Nonlinear", rownames(a))[1], "P"]
  cat(sprintf("  %-4s test for non-linearity p = %.2f\n", ep[3], p_nl))
}

# ---- Events per variable ----------------------------------------------------
# Six estimated parameters: exposure, sex, three IPI contrasts, age.
cat("\n==== Events per variable (6 parameters) ====\n")
cat(sprintf("  OS : %d/6 = %.1f\n  PFS: %d/6 = %.1f\n",
            sum(df$os_event), sum(df$os_event) / 6,
            sum(df$pfs_event), sum(df$pfs_event) / 6))
