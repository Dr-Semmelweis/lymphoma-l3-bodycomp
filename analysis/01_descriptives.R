# -----------------------------------------------------------------------------
# 01_descriptives.R
#
# Cohort description.
#
# Computes:
#   - demographics, disease characteristics, treatment and outcome counts;
#   - L3 body-composition metrics, overall and by sex, as median and IQR;
#   - Spearman correlations among the body-composition metrics;
#   - deaths and progression-free survival events by muscle-density tertile.
#     Tertiles are cut on ranks, so the groups are as equal as the sample allows.
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
stopifnot(nrow(clin) == nrow(img),
          identical(clin$patient_id, img$patient_id))

# Category labels this script reads. Recode your data to these labels, or
# change them here; docs/data_dictionary.md lists them for every table.
sex_female  <- "F"
sex_male    <- "M"
ipi_levels  <- c("Low", "L-Int", "H-Int", "High")
dead_label  <- "deceased"
lost_label  <- "lost_to_follow_up"
stopifnot(all(clin$sex %in% c(sex_male, sex_female)),
          all(clin$ipi_group %in% ipi_levels),
          dead_label %in% clin$vital_status)

# Formatting helpers: median (IQR) and n (%).
med_iqr <- function(x, d = 1) {
  q <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
  sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), q[2], q[1], q[3])
}
n     <- nrow(clin)
n_pct <- function(k, denom = n) sprintf("%d (%d%%)", k, round(100 * k / denom))

# ---- Participants -----------------------------------------------------------
cat(sprintf("==== Participants (n = %d) ====\n", n))
cat("Age, years, median (IQR):        ", med_iqr(clin$age_years), "\n")
cat("Female sex:                      ", n_pct(sum(clin$sex == sex_female)), "\n")
cat("Histology:\n")
print(table(clin$histology_group))
cat("Ann Arbor stage III-IV:          ", n_pct(sum(clin$ann_arbor_stage >= 3)), "\n")
cat("Elevated LDH:                    ", n_pct(sum(clin$ldh_elevated)), "\n")
cat("B symptoms:                      ", n_pct(sum(clin$b_symptoms)), "\n")
cat("ECOG >= 2:                       ", n_pct(sum(clin$ecog >= 2)), "\n")
cat(sprintf("IPI group (%s):\n", paste(ipi_levels, collapse = " / ")))
print(table(factor(clin$ipi_group, levels = ipi_levels)))
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
cat("Reconstruction protocol:\n")
print(table(clin$reconstruction_protocol))

# ---- Outcomes and follow-up -------------------------------------------------
# Follow-up is summarised in two ways: the observed distribution of overall-
# survival time, and the reverse Kaplan-Meier estimate, in which the censoring
# indicator is inverted so that censoring becomes the event of interest.
cat("\n==== Outcomes and follow-up ====\n")
cat("Deaths (OS events):              ", sum(clin$os_event), "\n")
cat("  by recorded cause:\n")
print(table(clin$death_cause[clin$vital_status == dead_label]))
cat("PFS events:                      ", sum(clin$pfs_event), "\n")
cat("  documented progression or relapse:     ",
    sum(clin$pfs_event == 1 & clin$progression_or_relapse == 1), "\n")
cat("  death without documented progression:  ",
    sum(clin$pfs_event == 1 & clin$progression_or_relapse == 0), "\n")
cat("Lost to follow-up:               ",
    n_pct(sum(clin$vital_status == lost_label)), "\n")
cat("Follow-up, months, median (IQR): ", med_iqr(clin$os_months), "\n")
rkm <- survfit(Surv(os_months, 1 - os_event) ~ 1, data = clin)
cat("Reverse Kaplan-Meier median follow-up, months: ",
    sprintf("%.1f", summary(rkm)$table["median"]), "\n")

# ---- Table 2: body composition, overall and by sex --------------------------
cat("\n==== Body composition at L3 ====\n")
metrics <- c(sma_cm2 = 1, smd_mean_hu = 1, smi_cm2_m2 = 1, smg_au = 0,
             nama_cm2 = 1, lama_cm2 = 1, nama_tama_index_pct = 1,
             tama_density_hu = 1, sat_cm2 = 1, vat_cm2 = 1, imat_cm2 = 1)
men   <- clin$sex == sex_male
women <- clin$sex == sex_female
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
# Rank-based tertiles, larger groups first when the cohort does not divide
# evenly.
r   <- rank(img$smd_mean_hu, ties.method = "first")
ter <- factor(floor(3 * (r - 1) / length(r)) + 1,
              labels = c("lowest", "middle", "highest"))
print(rbind(n         = table(ter),
            deaths     = tapply(clin$os_event,  ter, sum),
            pfs_events = tapply(clin$pfs_event, ter, sum)))
