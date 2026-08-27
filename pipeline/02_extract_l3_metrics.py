"""Extract L3 body-composition metrics from the TotalSegmentator masks.

For each patient: the analysis slice is the median z of the L3 vertebral
mask; every tissue mask is intersected with the trunk mask (body_trunc) to
exclude the arms; tissues are then quantified within established Hounsfield
ranges. Areas are pixel counts times the in-plane pixel area.

Hounsfield windows:
  skeletal muscle        -29..150   (Mourtzakis 2008)
  NAMA                    30..150   (Aubrey 2014)
  LAMA                   -29..29
  adipose (SAT/VAT/IMAT) -190..-30
TAMA = NAMA + LAMA + IMAT. VAT is torso_fat within the adipose window; IMAT
comes from the dedicated tissue_4_types mask. The unmasked values quantify
the overestimation that arm inclusion would cause and serve as quality
control only.

Input layout: --nifti-dir and --seg-dir as produced by 01_run_segmentation.py.
Output: one CSV row per patient (--out).

Usage:  python 02_extract_l3_metrics.py --nifti-dir <dir> --seg-dir <dir> --out metrics.csv
"""
import argparse
import csv
from pathlib import Path

import nibabel as nib
import numpy as np
import scipy.stats

FIELDS = [
    "patient_id", "l3_slice_index", "l3_mask_n_slices", "pixel_area_cm2",
    "voxel_z_spacing_mm",
    "sma_cm2", "smd_mean_hu", "smd_median_hu", "smd_sd_hu", "smd_skewness",
    "nama_cm2", "lama_cm2", "tama_cm2", "tama_density_hu", "nama_tama_index_pct",
    "sat_cm2", "sat_density_hu", "vat_cm2", "vat_density_hu", "imat_cm2",
    "tat_cm2", "vat_sat_ratio", "total_body_area_cm2",
    "sma_unmasked_cm2", "smd_unmasked_hu",
    "trunk_mask_effect_sma_pct", "trunk_mask_effect_smd_hu", "arms_in_fov",
]


def mask_slice(path: Path, z: int) -> np.ndarray:
    return np.asarray(nib.load(path).dataobj[:, :, z]) > 0


def extract(pid: str, nifti_dir: Path, seg_dir: Path) -> dict:
    seg = seg_dir / pid
    nii = sorted((nifti_dir / pid).glob("*.nii.gz"))
    if len(nii) != 1:
        return {"patient_id": pid, "error": "expected exactly one CT NIfTI"}
    ct = nib.load(nii[0])
    sp = ct.header.get_zooms()
    px = (sp[0] * sp[1]) / 100                      # cm^2 per pixel

    l3 = nib.load(seg / "total" / "vertebrae_L3.nii.gz").get_fdata()
    if not np.any(l3 > 0):
        return {"patient_id": pid, "error": "empty L3 mask"}
    zs = np.where(l3 > 0)[2]
    z = int(np.median(zs))

    trunk = mask_slice(seg / "body" / "body_trunc.nii.gz", z)
    sm_raw = mask_slice(seg / "tissue_types" / "skeletal_muscle.nii.gz", z)
    sat_raw = mask_slice(seg / "tissue_types" / "subcutaneous_fat.nii.gz", z)
    vat_raw = mask_slice(seg / "tissue_types" / "torso_fat.nii.gz", z)
    imat_raw = mask_slice(seg / "tissue_4_types" / "intermuscular_fat.nii.gz", z)

    sm, sat, vat, imat_m = (m & trunk for m in (sm_raw, sat_raw, vat_raw, imat_raw))
    s = np.asarray(ct.dataobj[:, :, z])
    area = lambda m: float(np.sum(m) * px)

    sm_c = sm & (s >= -29) & (s <= 150)
    sat_c = sat & (s >= -190) & (s <= -30)
    vat_c = vat & (s >= -190) & (s <= -30)
    imat_c = imat_m & (s >= -190) & (s <= -30)
    nama = sm & (s >= 30) & (s <= 150)
    lama = sm & (s >= -29) & (s <= 29)

    sma, sat_a, vat_a = area(sm_c), area(sat_c), area(vat_c)
    nama_a, lama_a, imat_a = area(nama), area(lama), area(imat_c)
    tama_a = nama_a + lama_a + imat_a
    tama_mask = sm_c | imat_c

    # unmasked variant: same windows without the trunk restriction (QC only)
    sm_nm = sm_raw & (s >= -29) & (s <= 150)
    sma_nm = area(sm_nm)
    smd_nm = float(np.mean(s[sm_nm])) if sm_nm.any() else np.nan
    smd = float(np.mean(s[sm_c])) if sm_c.any() else np.nan

    return {
        "patient_id": pid,
        "l3_slice_index": z,
        "l3_mask_n_slices": int(len(np.unique(zs))),
        "pixel_area_cm2": round(px, 6),
        "voxel_z_spacing_mm": round(float(sp[2]), 2),
        "sma_cm2": sma,
        "smd_mean_hu": smd,
        "smd_median_hu": float(np.median(s[sm_c])) if sm_c.any() else np.nan,
        "smd_sd_hu": float(np.std(s[sm_c])) if sm_c.any() else np.nan,
        "smd_skewness": float(scipy.stats.skew(s[sm_c])) if sm_c.sum() > 100 else np.nan,
        "nama_cm2": nama_a,
        "lama_cm2": lama_a,
        "tama_cm2": tama_a,
        "tama_density_hu": float(np.mean(s[tama_mask])) if tama_mask.any() else np.nan,
        "nama_tama_index_pct": (nama_a / tama_a * 100) if tama_a > 0 else np.nan,
        "sat_cm2": sat_a,
        "sat_density_hu": float(np.mean(s[sat_c])) if sat_c.any() else np.nan,
        "vat_cm2": vat_a,
        "vat_density_hu": float(np.mean(s[vat_c])) if vat_c.any() else np.nan,
        "imat_cm2": imat_a,
        "tat_cm2": sat_a + vat_a,
        "vat_sat_ratio": (vat_a / sat_a) if sat_a > 0 else np.nan,
        "total_body_area_cm2": sma + sat_a + vat_a + imat_a,
        "sma_unmasked_cm2": sma_nm,
        "smd_unmasked_hu": smd_nm,
        "trunk_mask_effect_sma_pct": (sma_nm / sma - 1) * 100 if sma > 0 else np.nan,
        "trunk_mask_effect_smd_hu": smd_nm - smd,
        "arms_in_fov": bool(sma_nm > sma * 1.05),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nifti-dir", required=True, type=Path)
    ap.add_argument("--seg-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    ids = sorted(d.name for d in args.seg_dir.iterdir() if d.is_dir())
    rows, failed = [], 0
    for i, pid in enumerate(ids, 1):
        try:
            r = extract(pid, args.nifti_dir, args.seg_dir)
        except Exception as e:
            r = {"patient_id": pid, "error": f"{type(e).__name__}: {e}"}
        rows.append(r)
        if "error" in r:
            failed += 1
            print(f"[{i}/{len(ids)}] {pid}: {r['error']}")

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS + ["error"], extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows) - failed}/{len(rows)} patients extracted -> {args.out}")


if __name__ == "__main__":
    main()
