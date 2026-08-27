# -----------------------------------------------------------------------------
# 01_descriptives.R
#
# Cohort description and body-composition summaries.
#
# Reproduces:
#   - the participant and outcome figures reported in Results ("Participants")
#     and in Table 1: demographics, histology, stage, IPI, cell of origin,
#     treatment, end-of-treatment response, deaths by cause, follow-up;
#   - Table 2: L3 body-composition metrics, overall and by sex;
#   - the Spearman correlations quoted in Results ("Body composition and
#     segmentation quality");
#   - deaths and progression-free survival events by muscle-density tertile.
#
# Input:  data/clinical_dataset_n155.csv
#         data/imaging_L3_quantitative_n155.csv
# Output: printed to standard output; no files are written.
#
# Run from the package root:  Rscript code/01_descriptives.R
# Requires R >= 4.5 and the 'survival' package.
# -----------------------------------------------------------------------------

library(survival)

clin <- read.csv(file.path("data", "clinical_dataset_n155.csv"), na.strings = c("NA", ""))
img  <- read.csv(file.path("data", "imaging_L3_quantitative_n155.csv"), na.strings = c("NA", ""))
stopifnot(nrow(clin) == 155, nrow(img) == 155,
          identical(clin$patient_id, img$patient_id))

# Formatting helpers: median (IQR) and n (%).
med_iqr <- function(x, d = 1) {
  q <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
  sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), q[2], q[1], q[3])
}
n_pct <- function(k, n = 155) sprintf("%d (%d%%)", k, round(100 * k / n))

# ---- Participants -----------------------------------------------------------
cat("==== Participants (n = 155) ====\n")
cat("Age, years, median (IQR):        ", med_iqr(clin$age_years), "\n")
cat("Female sex:                      ", n_pct(sum(clin$sex == "F")), "\n")
cat("Histology:\n")
print(table(clin$histology_group))
cat("Ann Arbor stage III-IV:          ", n_pct(sum(clin$ann_arbor_stage >= 3)), "\n")
cat("Elevated LDH:                    ", n_pct(sum(clin$ldh_elevated)), "\n")
cat("B symptoms:                      ", n_pct(sum(clin$b_symptoms)), "\n")
cat("ECOG >= 2:                       ", n_pct(sum(clin$ecog >= 2)), "\n")
cat("IPI group (Low / L-Int / H-Int / High):\n")
print(table(factor(clin$ipi_group, levels = c("Low", "L-Int", "H-Int", "High"))))
cat("Cell of origin (Hans):\n")
print(table(clin$coo_hans, useNA = "ifany"))
cat("Double expressor:                ", n_pct(sum(clin$double_expressor, na.rm = TRUE)), "\n")
cat("Ki-67, %, median (IQR):          ", med_iqr(clin$ki67_pct, 0), "\n")
cat("BMI, kg/m2, median (IQR):        ", med_iqr(clin$bmi), "\n")
cat("Treatment intensity:\n")
print(table(clin$treatment_intensity, useNA = "ifany"))
cat("End-of-treatment complete response: ",
    sprintf("%d/%d (%d%%)\n", sum(clin$eot_complete_response, na.rm = TRUE),
            sum(!is.na(clin$eot_complete_response)),
            round(100 * mean(clin$eot_complete_response, na.rm = TRUE))))
cat("Reconstruction protocol 5 mm / 3 mm: ",
    sum(clin$reconstruction_protocol == "5mm"), "/",
    sum(clin$reconstruction_protocol == "3mm"), "\n")

# ---- Outcomes and follow-up -------------------------------------------------
# Follow-up is summarised in two ways: the observed distribution of overall-
# survival time, and the reverse Kaplan-Meier estimate, in which the censoring
# indicator is inverted so that censoring becomes the event of interest.
cat("\n==== Outcomes and follow-up ====\n")
cat("Deaths (OS events):              ", sum(clin$os_event), "\n")
cat("  by recorded cause:\n")
print(table(clin$death_cause[clin$vital_status == "deceased"]))
cat("PFS events:                      ", sum(clin$pfs_event), "\n")
cat("  documented progression or relapse:     ",
    sum(clin$pfs_event == 1 & clin$progression_or_relapse == 1), "\n")
cat("  death without documented progression:  ",
    sum(clin$pfs_event == 1 & clin$progression_or_relapse == 0), "\n")
cat("Lost to follow-up:               ",
    n_pct(sum(clin$vital_status == "lost_to_follow_up")), "\n")
cat("Follow-up, months, median (IQR): ", med_iqr(clin$os_months), "\n")
rkm <- survfit(Surv(os_months, 1 - os_event) ~ 1, data = clin)
cat("Reverse Kaplan-Meier median follow-up, months: ",
    sprintf("%.1f", summary(rkm)$table["median"]), "\n")

# ---- Table 2: body composition, overall and by sex --------------------------
cat("\n==== Body composition at L3 (Table 2) ====\n")
metrics <- c(sma_cm2 = 1, smd_mean_hu = 1, smi_cm2_m2 = 1, smg_au = 0,
             nama_cm2 = 1, lama_cm2 = 1, nama_tama_index_pct = 1,
             tama_density_hu = 1, sat_cm2 = 1, vat_cm2 = 1, imat_cm2 = 1)
men   <- clin$sex == "M"
women <- clin$sex == "F"
cat(sprintf("%-22s %-22s %-22s %-22s\n", "Metric",
            sprintf("Overall (n=%d)", nrow(img)),
            sprintf("Men (n=%d)", sum(men)),
            sprintf("Women (n=%d)", sum(women))))
for (m in names(metrics)) {
  d <- metrics[[m]]
  cat(sprintf("%-22s %-22s %-22s %-22s\n", m,
              med_iqr(img[[m]], d), med_iqr(img[[m]][men], d),
              med_iqr(img[[m]][women], d)))
}

# ---- Correlation structure --------------------------------------------------
# The pairs quoted in Results; Spearman rank correlations.
cat("\n==== Spearman correlations quoted in Results ====\n")
rho <- function(a, b) sprintf("%+.2f", cor(a, b, method = "spearman"))
cat("SMD ~ NAMA/TAMA index:  ", rho(img$smd_mean_hu, img$nama_tama_index_pct), "\n")
cat("SMD ~ mean TAMA attn.:  ", rho(img$smd_mean_hu, img$tama_density_hu), "\n")
cat("SMD ~ muscle area:      ", rho(img$smd_mean_hu, img$sma_cm2), "\n")
cat("SMD ~ intermuscular fat:", rho(img$smd_mean_hu, img$imat_cm2), "\n")
cat("SMD ~ age:              ", rho(img$smd_mean_hu, clin$age_years), "\n")

# ---- Events by muscle-density tertile ---------------------------------------
# Descriptive gradient reported alongside the Kaplan-Meier curves.
cat("\n==== Events by SMD tertile ====\n")
# Rank-based tertiles with the larger groups first (52, 52, 51), matching the
# grouping used in the original analysis.
r   <- rank(img$smd_mean_hu, ties.method = "first")
ter <- factor(floor(3 * (r - 1) / length(r)) + 1,
              labels = c("lowest", "middle", "highest"))
print(rbind(n         = table(ter),
            deaths     = tapply(clin$os_event,  ter, sum),
            pfs_events = tapply(clin$pfs_event, ter, sum)))
