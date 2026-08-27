"""Run TotalSegmentator on every patient CT.

Input layout: one sub-directory per patient under --nifti-dir, each containing
a single compressed NIfTI (*.nii.gz) with the non-contrast low-dose CT.
Output: masks under --seg-dir/<patient>/<task>/.

Four tasks are run per patient. The 'total' task is restricted to the L3
vertebra, which defines the analysis level; 'body' provides the trunk mask
used to exclude the arms; 'tissue_types' provides the muscle and fat masks;
'tissue_4_types' provides the dedicated intermuscular-fat mask. The
tissue_types and tissue_4_types tasks require a TotalSegmentator academic
licence.

Usage:  python 01_run_segmentation.py --nifti-dir <dir> --seg-dir <dir> [--device gpu]

Tested with TotalSegmentator 2.15.0 (pre-trained weights, no retraining),
PyTorch 2.10, CUDA 12.8, Python 3.12.
"""
import argparse
import subprocess
import sys
from pathlib import Path

TASKS = [
    ("total",          ["--roi_subset", "vertebrae_L3"]),
    ("body",           []),
    ("tissue_types",   []),
    ("tissue_4_types", []),
]


def find_ct(patient_dir: Path) -> Path | None:
    nii = sorted(patient_dir.glob("*.nii.gz"))
    return nii[0] if len(nii) == 1 else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nifti-dir", required=True, type=Path)
    ap.add_argument("--seg-dir", required=True, type=Path)
    ap.add_argument("--device", default="gpu")
    args = ap.parse_args()

    patients = sorted(d for d in args.nifti_dir.iterdir() if d.is_dir())
    failures = 0
    for i, pdir in enumerate(patients, 1):
        ct = find_ct(pdir)
        if ct is None:
            print(f"[{i}/{len(patients)}] {pdir.name}: expected exactly one *.nii.gz, skipped")
            failures += 1
            continue
        for task, extra in TASKS:
            out = args.seg_dir / pdir.name / task
            if out.exists() and any(out.glob("*.nii.gz")):
                continue  # already segmented
            cmd = ["TotalSegmentator", "-i", str(ct), "-o", str(out),
                   "--task", task, "--device", args.device] + extra
            p = subprocess.run(cmd, capture_output=True, text=True)
            if p.returncode != 0:
                print(f"[{i}/{len(patients)}] {pdir.name}/{task}: FAILED "
                      f"(exit {p.returncode})\n{(p.stderr or '')[-300:]}")
                failures += 1
                break
        else:
            print(f"[{i}/{len(patients)}] {pdir.name}: done")

    print(f"\n{len(patients) - failures}/{len(patients)} patients segmented")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
