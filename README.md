# lymphoma-l3-bodycomp

Automated L3 body composition from the low-dose CT of staging \[18F]FDG PET/CT, and its association with survival in aggressive B-cell lymphoma.

Code accompanying: *Automated CT body composition on staging PET/CT: myosteatosis and survival in aggressive B-cell lymphoma* (manuscript under review). Archived release: Zenodo DOI `10.5281/zenodo.XXXXXXX` (to be minted at first release).

## What is here

|Path|Purpose|Produces|
|-|-|-|
|`pipeline/01\_run\_segmentation.py`|Runs TotalSegmentator (tasks `total`, `body`, `tissue\_types`, `tissue\_4\_types`) on every patient CT|Segmentation masks|
|`pipeline/02\_extract\_l3\_metrics.py`|Single-slice L3 body-composition metrics with trunk restriction and fixed Hounsfield windows|One CSV row per patient|
|`analysis/01\_descriptives.R`|Cohort description|Tables 1–2|
|`analysis/02\_primary\_models.R`|Co-primary Cox models (OS, PFS)|Table 3, Figure 3|
|`analysis/03\_secondary\_predictors.R`|Nine secondary predictors, FDR correction, exploratory gauge|Table 4|
|`analysis/04\_sensitivity.R`|Five prespecified sensitivity analyses, interactions, post hoc|Results, Sensitivity|
|`analysis/05\_validation\_agreement.R`|Agreement with blinded manual segmentation (Dice, ICC, Bland–Altman)|Validation sub-study, Suppl. Fig. S3|



## Requirements

Pipeline: Python 3.12, TotalSegmentator 2.15.0 (pre-trained weights), PyTorch 2.10, CUDA 12.8, `nibabel`, `numpy`, `scipy`. The `tissue\_types` and `tissue\_4\_types` tasks require a TotalSegmentator academic licence.

Analysis: R 4.5 with `survival`, `rms`, `boot`, `logistf`, `gtsummary`. Package versions used for the manuscript are frozen in `analysis/sessionInfo.txt`.

## Usage

```bash
python pipeline/01\_run\_segmentation.py --nifti-dir <dir> --seg-dir <dir> \[--device gpu]
python pipeline/02\_extract\_l3\_metrics.py --nifti-dir <dir> --seg-dir <dir> --out metrics.csv
```

Input layout: one sub-directory per patient under `--nifti-dir`, each containing a single `\*.nii.gz` (non-contrast low-dose CT). The analysis scripts read the merged patient-level CSV described in `docs/data\_dictionary.md` from the working directory.

## License

MIT — see `LICENSE`.

