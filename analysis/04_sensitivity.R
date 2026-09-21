# -----------------------------------------------------------------------------
# 04_sensitivity.R
#
# Sensitivity analyses.
#
# Prespecified:
#   1. restriction to DLBCL-NOS and high-grade B-cell lymphoma;
#   2. survival measured from start of treatment rather than from staging;
#   3. the NAMA/TAMA index as the exposure;
#   4. informative censoring, as a Cox model on the censoring indicator;
#   5. restriction to one reconstruction protocol group.
#
# Post hoc:
#   6. exposure by reconstruction protocol, with the between-protocol
#      comparison of muscle density;
#   7. muscle density and end-of-treatment complete response, by Firth
#      penalised logistic regression;
#   8. exposure by sex.
#
# Interactions are tested by likelihood ratio. Z-scores are computed once on
# the full cohort and carried unchanged into every subset, so that "per 1 SD"
# denotes the same quantity throughout.
#
# Input:  clinical and imaging tables (docs/data_dictionary.md).
# Output: standard output.
#
# R >= 4.5; packages survival, logistf.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(survival)
  library(logistf)
})

# ---- Input ------------------------------------------------------------------
# Adapt these paths to your own data. The columns each table must provide are
# listed in docs/data_dictionary.md; nothing else in the script depends on how
# the files are named or where they are kept.
clin <- read.csv(file.path("data", "clinical_dataset_n155.csv"), na.strings = c("NA", ""))
img  <- read.csv(file.path("data", "imaging_L3_quantitative_n155.csv"), na.strings = c("NA", ""))
stopifnot(nrow(clin) == nrow(img), identical(clin$patient_id, img$patient_id))

# Category labels this script reads; recode your data or change them here.
# 'protocol_levels' takes the restricted group first, and 'histology_keep'
# the histologies retained in the first sensitivity analysis.
sex_levels      <- c("M", "F")
ipi_levels      <- c("Low", "L-Int", "H-Int", "High")
protocol_levels <- c("5mm", "3mm")
histology_keep  <- c("DLBCL-NOS", "HGBCL")
stopifnot(all(clin$sex %in% sex_levels), all(clin$ipi_group %in% ipi_levels),
          length(protocol_levels) == 2,
          all(clin$reconstruction_protocol %in% protocol_levels),
          any(clin$histology_group %in% histology_keep))

z  <- function(x) (x - mean(x)) / sd(x)
df <- data.frame(
  os_months  = clin$os_months,  os_event  = clin$os_event,
  pfs_months = clin$pfs_months, pfs_event = clin$pfs_event,
  os_tx  = clin$os_months_from_treatment,
  pfs_tx = clin$pfs_months_from_treatment,
  smd_z  = z(img$smd_mean_hu),
  nt_z   = z(img$nama_tama_index_pct),
  age_z  = z(clin$age_years),
  sex    = factor(clin$sex, levels = sex_levels),
  ipi    = factor(clin$ipi_group, levels = ipi_levels),
  histology = clin$histology_group,
  protocol  = factor(clin$reconstruction_protocol, levels = protocol_levels),
  eot       = clin$eot_complete_response,
  smd_hu    = img$smd_mean_hu
)

# Fits the standard adjusted model on the data provided and returns the
# formatted estimate for the exposure (column 'x').
run <- function(d, x, time, event) {
  d$x <- d[[x]]
  m <- coxph(as.formula(sprintf("Surv(%s, %s) ~ x + sex + ipi + age_z", time, event)),
             data = d)
  s <- summary(m)
  sprintf("HR %.2f (%.2f-%.2f)  p=%.3f  [n=%d, events=%d]",
          s$conf.int["x", 1], s$conf.int["x", 3], s$conf.int["x", 4],
          s$coefficients["x", 5], m$n, m$nevent)
}

both <- function(label, d, x = "smd_z", tos = "os_months", eos = "os_event",
                 tpfs = "pfs_months", epfs = "pfs_event") {
  cat(sprintf("%s\n  OS : %s\n  PFS: %s\n", label,
              run(d, x, tos, eos), run(d, x, tpfs, epfs)))
}

cat("==== Prespecified sensitivity analyses ====\n")

both(sprintf("1. Restriction to %s", paste(histology_keep, collapse = " + ")),
     df[df$histology %in% histology_keep, ])

# Time origin at the start of first-line treatment. The pre-computed columns
# are missing for patients without systemic treatment; non-positive
# times, which cannot enter a survival model, are excluded by the same filter.
both("2. Time origin at treatment start",
     df[!is.na(df$os_tx) & df$os_tx > 0 & !is.na(df$pfs_tx) & df$pfs_tx > 0, ],
     tos = "os_tx", tpfs = "pfs_tx")

both("3. NAMA/TAMA index as the exposure", df, x = "nt_z")

# Cox model on the censoring indicator: the "event" is being censored, and the
# question is whether muscle density predicts it. A null result supports
# censoring that is non-informative with respect to the exposure.
cens <- df; cens$cens <- 1 - cens$os_event
mc <- coxph(Surv(os_months, cens) ~ smd_z + sex + ipi + age_z, data = cens)
cat(sprintf("4. Censoring model: muscle density on the censoring indicator, p=%.2f\n",
            summary(mc)$coefficients["smd_z", 5]))

both(sprintf("5. Restriction to the %s protocol group", protocol_levels[1]),
     df[df$protocol == protocol_levels[1], ])

cat("\n==== Post-hoc analyses ====\n")

# 6. Interaction with the reconstruction protocol, tested on the whole cohort
# by likelihood-ratio comparison of nested models, and the direct comparison
# of muscle density between the two protocol groups.
cat("6. Muscle density x reconstruction protocol\n")
for (ep in list(c("os_months", "os_event", "OS"),
                c("pfs_months", "pfs_event", "PFS"))) {
  f0 <- as.formula(sprintf("Surv(%s, %s) ~ smd_z + sex + ipi + age_z + protocol", ep[1], ep[2]))
  f1 <- as.formula(sprintf("Surv(%s, %s) ~ smd_z * protocol + sex + ipi + age_z", ep[1], ep[2]))
  a  <- anova(coxph(f0, data = df), coxph(f1, data = df))
  cat(sprintf("  %-4s likelihood-ratio p=%.2f\n", ep[3], a[2, "Pr(>|Chi|)"]))
}
tt <- t.test(smd_hu ~ protocol, data = df)
cat(sprintf("  density %s vs %s: difference %.2f HU (%.2f SD), p=%.2f\n",
            protocol_levels[1], protocol_levels[2],
            diff(tt$estimate), diff(tt$estimate) / sd(df$smd_hu), tt$p.value))

# 7. End-of-treatment complete response, evaluable patients only. Firth
# penalisation is used because of the small number of non-responders.
de <- df[!is.na(df$eot), ]
fe <- logistf(eot ~ smd_z + sex + ipi + age_z, data = de)
i  <- which(names(coef(fe)) == "smd_z")
cat(sprintf("7. End-of-treatment response (Firth): OR %.2f (%.2f-%.2f)  p=%.2f  [n=%d, CR=%d]\n",
            exp(coef(fe)[i]), exp(fe$ci.lower[i]), exp(fe$ci.upper[i]),
            fe$prob[i], nrow(de), sum(de$eot)))

# 8. Interaction with sex, likelihood-ratio test.
cat("8. Muscle density x sex\n")
for (ep in list(c("os_months", "os_event", "OS"),
                c("pfs_months", "pfs_event", "PFS"))) {
  f0 <- as.formula(sprintf("Surv(%s, %s) ~ smd_z + sex + ipi + age_z", ep[1], ep[2]))
  f1 <- as.formula(sprintf("Surv(%s, %s) ~ smd_z * sex + ipi + age_z", ep[1], ep[2]))
  a  <- anova(coxph(f0, data = df), coxph(f1, data = df))
  cat(sprintf("  %-4s likelihood-ratio p=%.2f\n", ep[3], a[2, "Pr(>|Chi|)"]))
}
