# -----------------------------------------------------------------------------
# 05_validation_agreement.R
#
# Agreement between the automated pipeline and blinded manual segmentation.
#
# Computes:
#   - Dice and Jaccard coefficients from the mask pixel counts;
#   - ICC(A,1), two-way random effects, absolute agreement, single measures,
#     with percentile bootstrap intervals;
#   - Bland-Altman analysis for area and density: mean difference with its
#     paired test, 95% limits of agreement, and proportional bias as the slope
#     of the differences on the pair means;
#   - the width of the limits of agreement in units of the between-patient SD
#     of the cohort;
#   - intra-observer repeatability from the repeated readings.
#
# Input:  the two validation tables, plus the imaging table for the cohort SDs
#         (docs/data_dictionary.md).
# Output: standard output.
#
# R >= 4.5; base R only.
# -----------------------------------------------------------------------------

# ---- Input ------------------------------------------------------------------
# Adapt these paths to your own data. The columns each table must provide are
# listed in docs/data_dictionary.md; nothing else in the script depends on how
# the files are named or where they are kept.
val <- read.csv(file.path("data", "manual_validation_n30.csv"), na.strings = c("NA", ""))
rep <- read.csv(file.path("data", "manual_repeatability_n10.csv"), na.strings = c("NA", ""))
img <- read.csv(file.path("data", "imaging_L3_quantitative_n155.csv"), na.strings = c("NA", ""))
stopifnot(nrow(val) > 0, nrow(rep) > 0, nrow(img) > 0,
          all(rep$patient_id %in% val$patient_id))

# Intraclass correlation coefficient ICC(A,1): two-way random effects,
# absolute agreement, single measures (McGraw & Wong, case 2A).
icc_a1 <- function(a, b) {
  Y <- cbind(a, b); n <- nrow(Y); k <- 2; gm <- mean(Y)
  msr <- k * sum((rowMeans(Y) - gm)^2) / (n - 1)
  msc <- n * sum((colMeans(Y) - gm)^2) / (k - 1)
  mse <- sum((Y - outer(rowMeans(Y), rep(1, k)) -
              outer(rep(1, n), colMeans(Y)) + gm)^2) / ((n - 1) * (k - 1))
  (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
}

# Percentile bootstrap interval for the ICC (patients resampled with
# replacement; the seed makes the interval reproducible).
icc_ci <- function(a, b, R = 5000, seed = 7) {
  set.seed(seed); n <- length(a)
  q <- replicate(R, { i <- sample.int(n, n, replace = TRUE); icc_a1(a[i], b[i]) })
  quantile(q, c(.025, .975))
}

# ---- Spatial overlap --------------------------------------------------------
dice    <- 2 * val$intersection_px / (val$automated_mask_px + val$manual_mask_px)
jaccard <- val$intersection_px /
           (val$automated_mask_px + val$manual_mask_px - val$intersection_px)
cat(sprintf("==== Spatial overlap (n = %d) ====\n", nrow(val)))
cat(sprintf("Dice:    median %.3f (IQR %.3f-%.3f), minimum %.3f\n",
            median(dice), quantile(dice, .25), quantile(dice, .75), min(dice)))
cat(sprintf("Jaccard: median %.3f\n\n", median(jaccard)))

# ---- Agreement for area and density ----------------------------------------
# Bland-Altman analysis on the manual-minus-automated differences: systematic
# bias with its paired test, 95% limits of agreement (bias +/- 1.96 SD), and
# proportional bias as the slope of the differences on the pair means.
agreement <- function(auto, man, unit, cohort_sd, label) {
  d <- man - auto; m <- (man + auto) / 2
  bias <- mean(d); sdd <- sd(d)
  lo <- bias - 1.96 * sdd; hi <- bias + 1.96 * sdd
  tt <- t.test(man, auto, paired = TRUE)
  prop <- summary(lm(d ~ m))$coefficients[2, 4]
  i  <- icc_a1(auto, man); ci <- icc_ci(auto, man)
  cat(sprintf("==== %s ====\n", label))
  cat(sprintf("automated %.2f +/- %.2f | manual %.2f +/- %.2f %s\n",
              mean(auto), sd(auto), mean(man), sd(man), unit))
  cat(sprintf("bias %+.2f %s (%+.2f%%), paired t p=%.2f\n",
              bias, unit, 100 * bias / mean(auto), tt$p.value))
  cat(sprintf("95%% limits of agreement %+.2f to %+.2f %s\n", lo, hi, unit))
  cat(sprintf("proportional bias: slope p=%.2f\n", prop))
  cat(sprintf("ICC(A,1) %.3f (bootstrap 95%% CI %.3f-%.3f)\n", i, ci[1], ci[2]))
  cat(sprintf("LoA width / cohort SD (%.1f %s): %.2f\n\n",
              cohort_sd, unit, (hi - lo) / cohort_sd))
}
agreement(val$automated_sma_cm2, val$manual_sma_cm2, "cm2",
          sd(img$sma_cm2), "Skeletal muscle area")
agreement(val$automated_smd_hu, val$manual_smd_hu, "HU",
          sd(img$smd_mean_hu), "Skeletal muscle density")

# ---- Intra-observer repeatability ------------------------------------------
# Images re-segmented after at least two weeks. The per-pair coefficient
# of variation uses the SD of the two readings (|difference| / sqrt(2)).
d2 <- rep$second_manual_sma_cm2 - rep$first_manual_sma_cm2
cv <- 100 * abs(d2) / sqrt(2) /
      ((rep$first_manual_sma_cm2 + rep$second_manual_sma_cm2) / 2)
cat(sprintf("==== Intra-observer repeatability (n = %d) ====\n", nrow(rep)))
cat(sprintf("coefficient of variation: mean %.2f%%, median %.2f%%, max %.2f%%\n",
            mean(cv), median(cv), max(cv)))
cat(sprintf("ICC(A,1) %.4f | bias second - first %+.2f cm2\n",
            icc_a1(rep$first_manual_sma_cm2, rep$second_manual_sma_cm2), mean(d2)))
