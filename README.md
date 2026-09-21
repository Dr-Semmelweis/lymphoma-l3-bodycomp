# lymphoma-l3-bodycomp

Automated L3 body composition from the low-dose CT of staging [18F]FDG PET/CT, and its association with survival in aggressive B-cell lymphoma.

Code accompanying *Automated CT body composition on staging PET/CT: myosteatosis and survival in aggressive B-cell lymphoma*. Archived on Zenodo: [10.5281/zenodo.22121857](https://doi.org/10.5281/zenodo.22121857).

## Method

**Images.** One non-contrast CT per patient: the low-dose series acquired for attenuation correction of the staging PET/CT, converted to a single compressed NIfTI. Slices must be axial, with the cranio-caudal direction on the third axis; the extraction script stops on volumes that are not.

**Segmentation.** Four TotalSegmentator tasks, pre-trained weights, no retraining: `total` restricted to `vertebrae_L3`, `body` for the trunk mask, `tissue_types` for muscle and fat, `tissue_4_types` for intermuscular fat.

**Analysis slice.** The median z of the L3 vertebral mask. Every tissue mask is intersected with the trunk mask: on whole-body PET/CT the arms lie within the abdominal field of view, where they are labelled as muscle. In our cohort this affected 92 of 155 patients and overestimated skeletal muscle area by a median of 22%, up to 65% in a single patient. The unmasked values are reported alongside the masked ones.

**Tissue quantification.** Hounsfield windows applied within the masks: skeletal muscle −29 to 150; normal-attenuation muscle 30 to 150; low-attenuation muscle −29 to 29; adipose tissue −190 to −30. Total abdominal muscle area is the sum of normal-attenuation muscle, low-attenuation muscle and intermuscular fat. Areas are pixel counts times the in-plane pixel area.

**Analysis.** The exposure is mean muscle attenuation at L3, standardised within the cohort under analysis; hazard ratios are therefore per 1 within-cohort standard deviation and do not depend on our absolute values. Overall and progression-free survival are co-primary and no correction is applied across them; the nine secondary predictors are corrected for false discovery rate within each endpoint. Models are adjusted for sex, IPI group and age as a continuous term.

## Scripts

|Script|Computes|
|-|-|
|`pipeline/01_run_segmentation.py`|the four TotalSegmentator tasks|
|`pipeline/02_extract_l3_metrics.py`|L3 metrics, trunk-restricted, fixed Hounsfield windows|
|`analysis/01_descriptives.R`|summaries, Spearman correlations, events by tertile|
|`analysis/02_primary_models.R`|co-primary Cox models, BCa bootstrap, restricted cubic splines|
|`analysis/03_secondary_predictors.R`|nine secondary predictors, FDR within endpoint|
|`analysis/04_sensitivity.R`|five prespecified and three post hoc analyses|
|`analysis/05_validation_agreement.R`|Dice, ICC, Bland–Altman, intra-observer repeatability|

## Input

Scripts `01`–`04` read a clinical table and a table of L3 metrics; script `05` reads the two validation tables and the imaging table. Required columns are in `docs/data_dictionary.md`. The clinical and imaging tables must have the same number of rows, with `patient_id` in the same order.

Two columns of the imaging table are not produced by the pipeline and must be computed from patient height: `smi_cm2_m2` and `smg_au`.

Paths are set in the `Input` block at the top of each script.

## Requirements

Pipeline: Python 3.12, TotalSegmentator 2.15.0, PyTorch 2.10, CUDA 12.8, nibabel, numpy, scipy. Conversion from DICOM to NIfTI happens before the pipeline; we used dcm2niix 1.0.20250505. The `tissue_types` and `tissue_4_types` tasks require a TotalSegmentator academic licence, free for non-commercial research.

Analysis: R 4.5 with `survival`, `rms`, `Hmisc`, `boot`, `logistf`. Versions used for the manuscript are in `analysis/sessionInfo.txt`.

## Usage

```bash
python pipeline/01_run_segmentation.py --nifti-dir <dir> --seg-dir <dir> [--device gpu]
python pipeline/02_extract_l3_metrics.py --nifti-dir <dir> --seg-dir <dir> --out metrics.csv
Rscript analysis/02_primary_models.R
```

One sub-directory per patient under `--nifti-dir`, each containing a single `*.nii.gz`.

## Citation

Cite via the Zenodo DOI until the article is published; see `CITATION.cff`.

## License

MIT — see `LICENSE`.

