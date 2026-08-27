# lymphoma-l3-bodycomp

Automated L3 body composition from the low-dose CT of staging [18F]FDG PET/CT, and its association with survival in aggressive B-cell lymphoma.

Code accompanying: *Automated CT body composition on staging PET/CT: myosteatosis and survival in aggressive B-cell lymphoma* (manuscript under review). Archived release: Zenodo DOI `10.5281/zenodo.XXXXXXX` (to be minted at first release).

## What is here

| Path | Purpose | Produces |
|---|---|---|
| `pipeline/01_run_segmentation.py` | Runs TotalSegmentator (tasks `total`, `body`, `tissue_types`, `tissue_4_types`) on every patient CT | Segmentation masks |
| `pipeline/02_extract_l3_metrics.py` | Single-slice L3 body-composition metrics with trunk restriction and fixed Hounsfield windows | One CSV row per patient |
| `analysis/01_descriptives.R` | Cohort description | Tables 1–2 |
| `analysis/02_primary_models.R` | Co-primary Cox models (OS, PFS) | Table 3, Figure 3 |
| `analysis/03_secondary_predictors.R` | Nine secondary predictors, FDR correction, exploratory gauge | Table 4 |
| `analysis/04_sensitivity.R` | Five prespecified sensitivity analyses, interactions, post hoc | Results, Sensitivity |
| `analysis/05_validation_agreement.R` | Agreement with blinded manual segmentation (Dice, ICC, Bland–Altman) | Validation sub-study, Suppl. Fig. S3 |

Patient-level data are **not** included and cannot be shared (see the manuscript's Data availability statement). `docs/data_dictionary.md` documents the exact structure the scripts expect.

## Requirements

Pipeline (tested): Python 3.12, TotalSegmentator 2.15.0 (pre-trained weights, no retraining), PyTorch 2.10, CUDA 12.8, `nibabel`, `numpy`, `scipy`. The `tissue_types` and `tissue_4_types` tasks require a TotalSegmentator academic licence (free for non-commercial research; see the TotalSegmentator repository).

Analysis: R 4.5 with `survival`, `rms`, `boot`, `logistf`. Package versions used for the manuscript are frozen in `analysis/sessionInfo.txt`.

## Usage

```bash
python pipeline/01_run_segmentation.py --nifti-dir <dir> --seg-dir <dir> [--device gpu]
python pipeline/02_extract_l3_metrics.py --nifti-dir <dir> --seg-dir <dir> --out metrics.csv
```

Input layout: one sub-directory per patient under `--nifti-dir`, each containing a single `*.nii.gz` (non-contrast low-dose CT). The analysis scripts read the merged patient-level CSV described in `docs/data_dictionary.md` from the working directory.

## Citation

Until the article is published, please cite this repository via its Zenodo DOI. A `CITATION.cff` file is provided.

## License

MIT — see `LICENSE`.
