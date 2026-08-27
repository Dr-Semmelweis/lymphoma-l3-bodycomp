# Data dictionary

Patient-level data are **not** distributed with this repository. The tables below describe the structure the scripts expect, so that the pipeline and the analyses can be run on an equivalently formatted dataset.

## Clinical database

|Column|Description|
|---|---|
|`patient_id`|Pseudonymised patient identifier.|
|`sex`|Sex.|
|`age_years`|Age at first staging.|
|`height_cm`|Height from the clinical record.|
|`weight_kg`|Body weight at staging.|
|`bmi`|Body-mass index.|
|`histology`|Histological diagnosis as recorded.|
|`histology_group`|Histology grouped for analysis.|
|`ann_arbor_stage`|Ann Arbor stage.|
|`b_symptoms`|B symptoms at diagnosis.|
|`ecog`|ECOG performance status.|
|`ldh_elevated`|Serum LDH above the institutional upper limit of normal (214 U/L).|
|`ipi_score`|International Prognostic Index, raw score.|
|`ipi_group`|IPI risk group.|
|`coo_hans`|Cell of origin, Hans algorithm.|
|`double_expressor`|Concurrent MYC and BCL2 overexpression.|
|`ki67_pct`|Ki-67 proliferation index.|
|`treatment_intensity`|First-line regimen intensity category.|
|`eot_complete_response`|Complete metabolic response at end of treatment (Lugano 2014, Deauville 1-3).|
|`progression_or_relapse`|Documented progression or relapse at any time during follow-up.|
|`vital_status`|Vital status at last contact.|
|`death_cause`|Recorded cause of death (deceased patients only).|
|`os_event`|Overall-survival event (death from any cause).|
|`os_months`|Time from first staging to death or last contact.|
|`pfs_event`|Progression-free-survival event (progression, relapse, or death).|
|`pfs_months`|Time from first staging to first PFS event or last contact.|
|`os_months_from_treatment`|Time from start of first-line treatment to death or last contact; sensitivity time origin. NA where no systemic treatment was recorded.|
|`pfs_months_from_treatment`|Time from start of first-line treatment to first PFS event or last contact.|
|`reconstruction_protocol`|CT reconstruction protocol of the baseline PET/CT.|

## Imaging database (pipeline output, L3 metrics)

|Column|Description|
|---|---|
|`patient_id`|Pseudonymised patient identifier.|
|`reconstruction_protocol`|CT reconstruction protocol.|
|`nominal_slice_thickness_mm`|Nominal reconstructed slice thickness.|
|`voxel_z_spacing_mm`|z spacing of the analysed NIfTI volume (reconstruction interval).|
|`pixel_area_cm2`|In-plane pixel area of the analysed grid.|
|`l3_slice_index`|Index of the analysis slice (median z of the L3 vertebral mask).|
|`l3_mask_n_slices`|Cranio-caudal extent of the L3 vertebral mask.|
|`sma_cm2`|Skeletal muscle area within -29..150 HU, trunk-restricted.|
|`smd_mean_hu`|Mean skeletal muscle density; primary exposure.|
|`smd_median_hu`|Median skeletal muscle attenuation.|
|`smd_sd_hu`|Standard deviation of muscle attenuation.|
|`smd_skewness`|Skewness of the muscle attenuation distribution.|
|`nama_cm2`|Normal-attenuation muscle area, 30..150 HU.|
|`lama_cm2`|Low-attenuation muscle area, -29..29 HU.|
|`tama_cm2`|Total abdominal muscle area = NAMA + LAMA + IMAT.|
|`tama_density_hu`|Mean attenuation over the TAMA mask.|
|`nama_tama_index_pct`|NAMA / TAMA.|
|`sat_cm2`|Subcutaneous adipose tissue, -190..-30 HU.|
|`sat_density_hu`|Mean SAT attenuation.|
|`vat_cm2`|Visceral adipose tissue, -190..-30 HU.|
|`vat_density_hu`|Mean VAT attenuation.|
|`imat_cm2`|Intermuscular adipose tissue (dedicated 4-class mask), -190..-30 HU.|
|`tat_cm2`|Total adipose tissue = SAT + VAT + IMAT.|
|`vat_sat_ratio`|VAT / SAT.|
|`total_body_area_cm2`|Whole-body cross-sectional area at L3.|
|`smi_cm2_m2`|Skeletal muscle index = SMA / height^2.|
|`smg_au`|Skeletal muscle gauge = SMI x SMD.|
|`sma_unmasked_cm2`|SMA without trunk restriction (quality-control value).|
|`smd_unmasked_hu`|SMD without trunk restriction (quality-control value).|
|`trunk_mask_effect_sma_pct`|Relative SMA overestimation when the trunk mask is omitted.|
|`arms_in_fov`|Arms detected within the abdominal field of view.|
|`qc_status`|Visual quality-control verdict.|

## Manual-segmentation validation database

|Column|Description|
|---|---|
|`patient_id`|Pseudonymised patient identifier.|
|`reconstruction_protocol`|CT reconstruction protocol.|
|`pixel_area_cm2`|In-plane pixel area.|
|`automated_mask_px`|Automated muscle mask size at the analysis slice.|
|`automated_sma_cm2`|Automated skeletal muscle area.|
|`automated_smd_hu`|Automated mean muscle density.|
|`manual_mask_px`|Manual muscle mask size, blinded radiologist, same slice.|
|`manual_sma_cm2`|Manual skeletal muscle area.|
|`manual_smd_hu`|Manual mean muscle density.|
|`intersection_px`|Pixels common to the two masks.|
|`auto_only_px`|Pixels in the automated mask only.|
|`manual_only_px`|Pixels in the manual mask only.|
|`first_manual_mask_px / _sma_cm2`|First manual reading (repeatability subset).|
|`second_manual_mask_px / _sma_cm2`|Second manual reading, at least two weeks later.|

## Analysis input (`clinical_dataset_n155.csv`)

One row per patient; the merge of the clinical database with the L3 metrics on `patient_id`. Column names match the dictionaries above.
