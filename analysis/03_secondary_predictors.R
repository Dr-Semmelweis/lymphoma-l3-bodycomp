# -----------------------------------------------------------------------------
# 03_secondary_predictors.R
#
# Secondary body-composition predictors.
#
# Computes hazard ratios per 1 within-cohort SD for the nine prespecified
# secondary predictors, on both endpoints, each modelled separately with the
# covariates of the co-primary models. Benjamini-Hochberg false-discovery-rate
# correction is applied within each endpoint, the two endpoints forming
# separate correction families.
#
# The skeletal muscle gauge is exploratory, lies outside the correction family
# and is reported with uncorrected p-values.
#
# The predictors are strongly intercorrelated; they estimate one correlated
# signal, not independent effects.
#
# Input:  clinical and imaging tables (docs/data_dictionary.md).
# Output: standard output.
#
# R >= 4.5; package survival.
# -----------------------------------------------------------------------------

library(survival)

# ---- Input ------------------------------------------------------------------
# Adapt these paths to your own data. The columns each table must provide are
# listed in docs/data_dictionary.md; nothing else in the script depends on how
# the files are named or where they are kept.
clin <- read.csv(file.path("data", "clinical_dataset_n155.csv"), na.strings = c("NA", ""))
img  <- read.csv(file.path("data", "imaging_L3_quantitative_n155.csv"), na.strings = c("NA", ""))
stopifnot(nrow(clin) == nrow(img), identical(clin$patient_id, img$patient_id))

# Category labels this script reads; recode your data or change them here.
sex_levels <- c("M", "F")
ipi_levels <- c("Low", "L-Int", "H-Int", "High")
stopifnot(all(clin$sex %in% sex_levels), all(clin$ipi_group %in% ipi_levels))

z  <- function(x) (x - mean(x)) / sd(x)
df <- data.frame(
  os_months  = clin$os_months,  os_event  = clin$os_event,
  pfs_months = clin$pfs_months, pfs_event = clin$pfs_event,
  age_z = z(clin$age_years),
  sex   = factor(clin$sex, levels = sex_levels),
  ipi   = factor(clin$ipi_group, levels = ipi_levels)
)

# The nine prespecified predictors (the correction family) and, listed last,
# the exploratory skeletal muscle gauge. Each is standardised within the
# cohort so that hazard ratios are per 1 SD.
predictors <- c(
  "Skeletal muscle area"           = "sma_cm2",
  "Normal-attenuation muscle area" = "nama_cm2",
  "Low-attenuation muscle area"    = "lama_cm2",
  "NAMA/TAMA index"                = "nama_tama_index_pct",
  "Mean TAMA attenuation"          = "tama_density_hu",
  "Subcutaneous adipose tissue"    = "sat_cm2",
  "Visceral adipose tissue"        = "vat_cm2",
  "Intermuscular adipose tissue"   = "imat_cm2",
  "VAT/SAT ratio"                  = "vat_sat_ratio"
)
smg <- "smg_au"

# Fits one adjusted model and returns the estimate for the predictor.
fit_one <- function(pred_z, time, event) {
  d <- df; d$x <- pred_z
  m <- coxph(as.formula(sprintf("Surv(%s, %s) ~ x + sex + ipi + age_z", time, event)),
             data = d)
  s <- summary(m)
  c(hr = s$conf.int["x", 1], lo = s$conf.int["x", 3],
    hi = s$conf.int["x", 4], p = s$coefficients["x", 5])
}

# All nine predictors on both endpoints; the false-discovery-rate correction
# is applied within each endpoint (nine tests per endpoint, not eighteen
# jointly), as prespecified.
run_family <- function(time, event) {
  est <- t(sapply(predictors, function(col) fit_one(z(img[[col]]), time, event)))
  cbind(est, fdr = p.adjust(est[, "p"], method = "BH"))
}
os  <- run_family("os_months",  "os_event")
pfs <- run_family("pfs_months", "pfs_event")

cat("==== Secondary predictors: HR per 1 SD, adjusted for sex, IPI, age ====\n")
cat(sprintf("%-32s %-24s %-8s %-24s %-8s\n",
            "Predictor", "OS HR (95% CI)", "OS FDR", "PFS HR (95% CI)", "PFS FDR"))
for (i in seq_along(predictors)) {
  cat(sprintf("%-32s %-24s %-8.3f %-24s %-8.3f\n", names(predictors)[i],
              sprintf("%.2f (%.2f-%.2f)", os[i, "hr"],  os[i, "lo"],  os[i, "hi"]),
              os[i, "fdr"],
              sprintf("%.2f (%.2f-%.2f)", pfs[i, "hr"], pfs[i, "lo"], pfs[i, "hi"]),
              pfs[i, "fdr"]))
}

# The skeletal muscle gauge (SMI x SMD) integrates muscle quantity and quality
# in a single measure; it was prespecified as exploratory and is therefore
# reported with uncorrected p-values, outside the correction family.
cat("\n==== Exploratory: skeletal muscle gauge (uncorrected) ====\n")
for (ep in list(c("os_months", "os_event", "OS"),
                c("pfs_months", "pfs_event", "PFS"))) {
  e <- fit_one(z(img[[smg]]), ep[1], ep[2])
  cat(sprintf("  %-4s HR %.2f (%.2f-%.2f)  uncorrected p = %.4f\n",
              ep[3], e["hr"], e["lo"], e["hi"], e["p"]))
}

# Raw p-values of the family, for completeness.
cat("\n==== Uncorrected p-values of the nine-predictor family ====\n")
for (i in seq_along(predictors)) {
  cat(sprintf("  %-32s OS p = %.4f   PFS p = %.4f\n",
              names(predictors)[i], os[i, "p"], pfs[i, "p"]))
}
