# goFlux reprocessing of the HF stem CH4 dataset — 25 Sep 2026

Record of the reprocessing done at the source (this repo) after three problems were found in
`ch4-data-filtering` (see its `WORKLOG_2026-09.md`). Numbers below are old → new, where "old"
is the current scripts run on the previous `flux_with_quality_flags.csv` and "new" is the same
scripts run on the goFlux-reprocessed dataset (legacy rows only, n = 1,640).

## What changed in the pipeline

1. **New `scripts/01_import/09_goflux_reprocess.R`.** Every closure with a raw 1-Hz trace is refit
   with goFlux (`best.flux`, LM/HM, Hüppi et al. 2018 criteria) and flagged with fluxqc 0.2.3
   (`flag_detection(precision = "mad", mdf = "wassmann", conf = 0.95)`, `qc_screens`).
   Windows: field-log comp start/end (LGR/UGGA, Jun 2023 – Mar 2025, xlsx logs from Sep 2024 on
   are read, including the Jan–Mar 2025 LGR days); LI-7810 REMARK span + 20 s deadband.
   Output `data/input/HF_2023-2025_tree_flux_goflux.csv` (2,008 closures; legacy and goFlux
   values side by side; `flux_source` says which fills the canonical columns).
2. **`09_quality_flags.R` → `10_quality_flags.R`**, rewritten. MDF = 1.96 σ / t · flux.term with σ =
   record-wide MAD per analyzer and t = closure seconds (the old script used t = nb.obs, per-closure
   Allan σ ×3·t_crit labelled "Christiansen", and datasheet precision). The old columns are kept
   under the same names for the sensitivity ridges, relabelled "Datasheet", "Campaign σ",
   "Per-closure σ ×3t".
3. **SE fix.** Legacy 2023-24 `CH4_SE` = SE_slope × mol in ppm with no /area; corrected ×1000/area
   (`CH4_SE_legacy`). goFlux LM SE / corrected legacy SE = 0.99 (LGR), i.e. the correction is right.
   For HM fits goFlux's SE is ~4× the LM SE; `CH4_SE` carries the SE of the selected model.
4. **Chamber volume.** Collar volume (`tree_volumes.csv`) + analyzer loop: 0.200 L LGR/UGGA
   (`Matthes_Lab/.../diurnal_flux_processing.R`, confirmed by the implied legacy volume on every
   tree), 0.028 L LI-7810. With these, legacy / goFlux LM flux on the same window = 1.01 (median;
   the 1 % is goFlux's dry-air (1 − H2O) term).
5. `06_filtering_snr.R` no longer multiplies `CH4_SE` by 1000; `01_bgs_model.R` and
   `03_ems_model_B.R` no longer hard-code the 282 h / 132 h windows (they moved to 285 h / 135 h).
6. EDI: `01_build_edi_package.R` now publishes the goFlux table first (attribute templates from
   `00_goflux_attributes.R`); the legacy table stays for provenance. `data/edi/` and
   `scripts/04_publish/` are gitignored in this repo, so those edits live only on disk.

## Coverage

| analyzer | year | closures | in legacy | fitted | legacy rows without trace | closures not in legacy (fitted) |
|---|---|---|---|---|---|---|
| LGR/UGGA | 2023 | 305 | 238 | 299 | 6 (2023-09-26 duplicates, no log) | 67 |
| LGR/UGGA | 2024 | 1078 | 904 | 1016 | 29 (2024-05-01, no LGR file) | 141 |
| LGR/UGGA | 2025 | 119 | 0 | 89 | 0 | 89 (Jan–Mar 2025 LGR days) |
| LI-7810 | 2025 | 506 | 498 | 440 | 64 (Oct 2025 raw file not archived; 3 unmatched) | 6 |

1,541 of the 1,640 legacy rows now carry a goFlux flux; 99 keep the legacy value. 368 additional
closures have a raw trace but were dropped by the legacy pipelines (the upstream 2023-24 output
had 179 NA-flux rows; Sep 2024 and Jan–Mar 2025 were never processed). They are in the goFlux CSV
with `in_legacy_dataset = FALSE` and are excluded from the analysis dataset unless
`GOFLUX_INCLUDE_NEW=1`.

## goFlux vs legacy (legacy rows with a trace)

| analyzer | season | n | r | median ratio | same sign | within 20 % | HM chosen | below MDF |
|---|---|---|---|---|---|---|---|---|
| LGR/UGGA | fall | 378 | 0.990 | 1.00 | 100 % | 78 % | 67 % | 62 % |
| LGR/UGGA | spring | 57 | 0.898 | 1.00 | 91 % | 72 % | 74 % | 67 % |
| LGR/UGGA | summer | 548 | 0.988 | 0.99 | 97 % | 66 % | 70 % | 40 % |
| LGR/UGGA | winter | 124 | 0.981 | 1.00 | 100 % | 77 % | 73 % | 81 % |
| LI-7810 | spring | 145 | 0.994 | 1.08 | 99 % | 48 % | 94 % | 10 % |
| LI-7810 | summer | 289 | 0.991 | 1.05 | 97 % | 57 % | 84 % | 2 % |

Means move more than medians because HM fits raise the large summer wetland fluxes
(LGR summer mean 2.52 → 2.97). LI-7810 fluxes are ~5 % higher than the legacy best-60-s-window
values and their SE is 3× smaller (full remark used). Campaign σ: 3.95 ppb (LGR), 0.139 ppb
(LI-7810). Full table: `outputs/tables/goflux_vs_legacy_by_analyzer_season.csv`.

## QC review (nothing clicked by hand)

`outputs/tables/goflux_qc_review.csv`: 528 closures (386 in the legacy set), sorted by number of
reasons. Screens fired: CO2 tracer 129 (slope ≤ 0 or undetected rise), convex trace 292, noisy
closure 113, C0 10, window < 60 s 1; field-log "bad" 64; legacy/goFlux LM ratio outside
0.75–1.33 (different window or volume upstream) 50; sign differs from legacy 13. The
ambient-start screen is off (field-log windows start after closure). A per-day clock-offset check
was tried and rejected for the same reason. Trace plots: `outputs/figures/goflux/`.

## Draft numbers that moved (old → new)

Results, temporal/spatial patterns
- Wetland mean 1.96 ± 0.27 → 2.29 ± 0.32 nmol m⁻² s⁻¹; median 0.15 → 0.15; IQR 0.04–0.59 → 0.04–0.67.
- Upland mean 0.05 ± 0.01 → 0.05 ± 0.01; median 0.02; IQR −0.003–0.08 → −0.004–0.08.
- ~40-fold → ~43-fold; positive fluxes 73.5 % / 90 % → 72.8 % / 89.0 %; upland negative 26.5 % → 27.2 %.
- N. sylvatica tree-level mean 5.67 ± 1.61 (median 4.17) → 6.68 ± 1.87 (median 5.20).
- Contrasts: N. sylvatica − A. rubrum 5.29 (p = 0.001) → 6.24 (p = 0.0008); − T. canadensis 5.48 (p < 0.001) → 6.48 (p = 0.0005); A. rubrum − T. canadensis p 0.99 → 0.99.
- Upland: A. rubrum > Q. rubra p 0.027 → 0.020; > T. canadensis p 0.036 → 0.034.
- A. rubrum wetland/upland 3.4-fold (p = 0.009) → 3.6-fold (p = 0.010); T. canadensis 6.3-fold (p = 0.16) → 5.4-fold (p = 0.26).
- August peak: wetland 7.04 ± 1.80 → 8.59 ± 2.18; upland 0.14 ± 0.03 → 0.16 ± 0.04. Minima unchanged (Jan 0.01; Apr 0.003).
- Summer 2024 wetland 5.30 ± 0.90 → 6.27 ± 1.10 (2023 0.78 ± 0.18 → 0.85 ± 0.19; 2025 0.75 ± 0.13 → 0.79 ± 0.14); N. sylvatica summer 14.8 / 1.7 / 2.08 → 17.6 / 1.85 / 2.19. Summer-2024 round peak 13.9 → 17.0.

Repeatability
- ICC range 0–0.11 (mean 0.05) → 0–0.096 (mean 0.047). N. sylvatica 0.108 → 0.094 (LRT p < 0.001); A. rubrum wetland 0.104 → 0.096 (p < 0.001), upland 0.052 → 0.053 (p = 0.056 → 0.051).
- Spearman ρ range 0.56–0.87 → 0.27–0.87. N. sylvatica 0.87 (p = 0.003) unchanged; A. rubrum upland 0.75 (p = 0.018) → 0.77 (p = 0.014); A. rubrum wetland 0.60 (p = 0.07) → 0.71 (p = 0.028, now significant); T. canadensis wetland 0.62 → 0.27 (p = 0.45).
- BLUP ranges (repeatability figure) shift with the larger wetland fluxes; recheck the figure caption values.

Environmental drivers (rolling correlations)
- Significant predictors 26 vs 5 → 26 vs 9 (upland gains CO2 dry NEON, CO2 dry Ha2, WTD BVS, pressure).
- Wetland strongest: CO2 concentration r = −0.26 → −0.25 (11.25 d); sensible heat −0.23 at 5.9 d → −0.22 at 5.75 d; FCO2 −0.21 at 9 h → −0.21 at 6 h; precipitation +0.21 → +0.20 (8.25 d); shallow SWC +0.21 at 2.6 d → +0.21 at 15 h; WTD +0.18–0.20 at 18 h unchanged.
- Variables significant at both sites 5 → 9 (adds CO2 dry NEON, CO2 dry Ha2, WTD BVS, pressure); wetland-only 21 → 17.
- Optimal windows: TS_Ha2 282 h → 285 h, WTD (BVS) 132 h → 135 h, SWC 129 h → 132 h, LE 123 h unchanged.

Wetland model
- Core model R² 65.4 % → 64.6 %, AIC 831.0 → 932.4, BIC 889.9 → 991.3; ICC 0.47 → 0.415.
- Core standardized effects: N. sylvatica temp 0.552 → 0.580, WTD 0.534 → 0.513, interaction 0.415 → 0.438; T. canadensis 0.078 / 0.074 / 0.038 → 0.080 / 0.071 / 0.036; A. rubrum 0.111 / 0.190 / 0.165 → 0.125 / 0.185 / 0.162.
- Alternatives: SWC-core R² 61.7 % (AIC 914.9, ΔAIC +83.8) → 61.2 % (AIC 1003.2, ΔAIC +70.8); LE-core 56.5 % (AIC 1016.0, ΔAIC +185.0) → 54.9 % (AIC 1112.5, ΔAIC +180.1).
- Full model R² 67.6 % (AIC 784.8, ΔAIC −46.3) → 67.7 % (AIC 883.9, ΔAIC −48.5); it now has 24 parameters (forward selection also kept FC_Ha1 and s10t), BIC no longer improves (+1.9).
- Full-model species effects: N. sylvatica temp 0.779 → 0.965, WTD 0.927 → 0.960, interaction 0.422 → 0.411; T. canadensis 0.025 / 0.154 / 0.084 → 0.012 / 0.157 / 0.085; A. rubrum 0.152 / 0.254 / 0.163 → 0.250 / 0.231 / 0.148. LE full −0.392 → −0.503; SWC full −0.281 → −0.300 (sign reversals hold).
- Predicted extremes (4–21 °C): N. sylvatica dry −0.18 to −0.04 → −0.59 to 0.41; wet 0.07 to 94.9 → −0.15 to 153; T. canadensis max 0.63 → 0.59; A. rubrum max 1.68 → 2.01.
- Species-averaged predictions (Discussion / Figure): dry 0.14–0.43 → 0.20–0.48; wet 0.02–5.69 → 0.00–6.39; upland max 0.079 → 0.09; wetland/upland maximum ratio 72× → 74×.

Upland model
- Model A R² 8.8 % (AIC −398) → 9.8 % (AIC −359.7); Model B R² 8.2 % (AIC −455) → 9.1 % (AIC −414.1); ΔR² 0.6 → 0.7 points; ΔAIC 57 → 54.
- ICC 0.067 / 0.067 → recheck (tree SD 0.042 → 0.046).
- A. rubrum temperature effect (Model B) 0.056 (p = 0.005) → 0.061 (t = 2.93); other species remain non-significant.

Tomography
- Upland ERT CV vs flux r = 0.597 (p < 0.001) → 0.571 (p = 0.001); Q. rubra 0.802 (p = 0.005) → 0.844 (p = 0.002); A. rubrum upland 0.665 (p = 0.036) → 0.611 (p = 0.060, no longer significant).
- Wetland r = −0.404 (p = 0.027) → −0.409 (p = 0.025); N. sylvatica −0.68 (p = 0.031) → −0.69 (p = 0.026).

Variance partitioning and detection
- Site share of variance 21.4 % → 21.9 % (recomputed; the 44.5/18.5/54.8 and 4.9/10.2/86.3 splits in the draft are not produced by any script in the repo and could not be reproduced: a nested species/tree model gives wetland 27.5 / 5.6 / 66.8 and upland 1.2 / 8.6 / 90.2 on the new data).
- Wetland/upland variance ratio 1,200 → 1,700.
- Below the reference MDF: 34.5 % → 40.0 % (LGR 45.4 → 54.1 %, LI-7810 9.4 → 7.6 %); negatives among detected fluxes 7.9 % → 7.2 %.

Methods text that no longer matches: fluxes are no longer "linear rate of change" only; the QC
paragraph (CO2 R² < 0.8 exclusion, SNR < 2, CH4 < −1 exclusion, "817 measurements") describes the
legacy pipeline; the measurement period runs to October 2025 with two analyzers.
