# -----------------------------------------------------------------------------
# 04_sensitivity.R
#
# Sensitivity analyses, prespecified and post hoc.
#
# Reproduces (Results, "Sensitivity analyses"):
#   Prespecified -
#   1. restriction to DLBCL-NOS and high-grade B-cell lymphoma;
#   2. survival measured from the start of treatment rather than from staging;
#   3. the NAMA/TAMA index as the exposure;
#   4. informative-censoring checks (Cox model on the censoring indicator);
#   5. restriction to the 5 mm reconstruction-protocol group.
#   Post hoc -
#   6. interaction between muscle density and reconstruction protocol, with the
#      between-protocol comparison of muscle density;
#   7. muscle density and end-of-treatment complete response (Firth penalised
#      logistic regression);
#   8. interaction between muscle density and sex.
#
# Standardisation note: z-scores are computed once, on the full cohort of 155,
# and carried unchanged into every subset, so that "per 1 SD" means the same
# quantity in every row of the output.
#
# Input:  data/clinical_dataset_n155.csv
#         data/imaging_L3_quantitative_n155.csv
# Output: printed to standard output; no files are written.
#
# Run from the package root:  Rscript code/04_sensitivity.R
# Requires R >= 4.5 and the packages 'survival', 'logistf'.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(survival)
  library(logistf)
})

clin <- read.csv(file.path("data", "clinical_dataset_n155.csv"), na.strings = c("NA", ""))
img  <- read.csv(file.path("data", "imaging_L3_quantitative_n155.csv"), na.strings = c("NA", ""))
stopifnot(nrow(clin) == 155, identical(clin$patient_id, img$patient_id))

z  <- function(x) (x - mean(x)) / sd(x)
df <- data.frame(
  os_months  = clin$os_months,  os_event  = clin$os_event,
  pfs_months = clin$pfs_months, pfs_event = clin$pfs_event,
  os_tx  = clin$os_months_from_treatment,
  pfs_tx = clin$pfs_months_from_treatment,
  smd_z  = z(img$smd_mean_hu),
  nt_z   = z(img$nama_tama_index_pct),
  age_z  = z(clin$age_years),
  sex    = factor(clin$sex, levels = c("M", "F")),
  ipi    = factor(clin$ipi_group, levels = c("Low", "L-Int", "H-Int", "High")),
  histology = clin$histology_group,
  protocol  = factor(clin$reconstruction_protocol, levels = c("5mm", "3mm")),
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

both("1. Restriction to DLBCL-NOS + HGBCL",
     df[df$histology %in% c("DLBCL-NOS", "HGBCL"), ])

# Time origin at the start of first-line treatment. The pre-computed columns
# are missing for the one patient without systemic treatment; non-positive
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

both("5. Restriction to the 5 mm protocol group", df[df$protocol == "5mm", ])

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
cat(sprintf("  density 5 mm vs 3 mm: difference %.2f HU (%.2f SD), p=%.2f\n",
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
