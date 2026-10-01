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
4. **Chamber volume.** Collar volume (`tree_volumes.csv`, which includes cap and 28.9 cm³ of
   tubing) + 0.028 L analyzer internal volume for both analyzers (user decision 30 Sep 2026, shared
   with ch4-data-filtering). Until 1 Oct 2026 this repo used 0.200 L for the LGR, copied from
   `Matthes_Lab/.../diurnal_flux_processing.R` (`lgr_volume <- .2`, no documented source), which made
   legacy / goFlux LM flux on the same window = 1.01. With 0.028 L the LGR fluxes, MDFs and flux
   terms scale by ~0.70 and the legacy comparison for 2023–24 shows that constant volume factor.
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

## Update 30 Sep 2026: window rules, per-day precision, trace QC, decay definitions

Changes after inspecting every CO₂/CH₄ trace (`scripts/01_import/12_trace_qc.R`; plots in
`outputs/figures/trace_qc/`):

- **Deadband.** The UGGA field-log windows started at chamber closure and included the placement
  transient (e.g. a ~100 ppb CH₄ jump in the first 15 s). Both analyzers now discard the first 20 s
  after closure (the LI-7810 already did). `13_deadband_sensitivity.R` refits every closure at 0–45 s:
  medians change ≤ 6 %, the N. sylvatica mean 6.81 → 6.63 (20 s) → 6.46 (30 s); conclusions unchanged.
- **Window end at chamber removal.** The LI-7810 remark often continued after the chamber was lifted
  (CO₂/CH₄ crash to ambient inside the window; e.g. tree 380 on 8 May 2025 fitted −0.15 with CH₄ clearly
  rising). Windows now end at an abrupt, sustained CO₂ fall (7-point running median; > max(10 ppm,
  8σ, 25 % of the rise) within 10 s; never recovers). 52 closures trimmed (41 UGGA, 11 LI-7810).
  A first, naive detector over-trimmed 181 closures and was replaced after visual checks.
- **Precision per analyzer × day** from the full day record (fluxqc::precision_mad_runs), matching the
  filtering paper; UGGA 2.9–6.9 ppb CH₄ across days (median 3.8), LI-7810 0.11–0.26 ppb.
- **Review list.** 40 analysis-set closures (21 above MDF) in
  `outputs/tables/trace_qc_clickpeak_shortlist.csv` / `traces_clickpeak_shortlist.pdf`; one suspect day
  (2024-05-24, 22 of 30 closures with window problems; team's "Timing Updates.xlsx" flags that period)
  that needs a single clock correction rather than clicking.
- **Decay definitions** (`scripts/02_analysis/09_decay_definitions.R`): the within-species relationships
  in the two site specialists hold across ERT definitions and in rank tests (N. sylvatica r ≈ −0.67 to
  −0.69; Q. rubra r = 0.71 index, 0.85 CV); SoT structural loss is unrelated in every group; pooled
  site-level correlations do not survive species adjustment (mixed models on all measurements). Text and
  Figure 5 were reframed to the within-species result.
- Optimal driver windows returned to the preprint values (282 h / 132 h). Final headline numbers:
  wetland 2.23 ± 0.30 vs upland 0.06 ± 0.01 (~40-fold); N. sylvatica 6.50 ± 1.85; core wetland model
  R² 66.1 %, ICC 0.48; upland models 9.2 % / 7.7 %; variance shares wetland 27.8/6.5/65.7, upland
  1.4/8.1/90.4. Full list: `outputs/tables/manuscript/manuscript_numbers.txt`.

## 1 Oct 2026 update: manual review, decay metric, driver models

- **Manual review applied.** 55 shortlisted closures reviewed by hand (`data/input/manual_windows.csv`):
  50 kept, 2 refitted over a clicked window, 3 excluded (initial concentration spike). Analysis n = 1,637.
- **Decay metric.** ERT CV is the primary within-species metric (Figure 5 strips and scatters).
  Tomography-paper classes stay as strip labels. Q. rubra: r = 0.79 with CV (robust; LOO 0.67–0.85);
  r = 0.54 with the species-normalized PC1, which depends on tree 175. Legacy fluxes give the same
  correlations (`10_decay_robustness.R`), so the earlier r = 0.71 under PC1 came from runs that ignored the
  team's corrected start times. Figure 5 fit lines and stars at p < 0.1.
- **Time base.** Flux times were joined to drivers 5 h apart in the screening vs the models; NEON SWC was
  in UTC while every other series is EST; 146 fluxes had a 00:00 placeholder time and 3 an AM/PM slip.
  `10_quality_flags.R` now writes `sample_hour_est` (EST hour of sampling; analyzer clock minus its offset
  where no real time was logged), used by every join; `08_align.R` shifts NEON to EST.
- **Screening.** Within-tree asinh flux (deviation from tree mean); rolling means need ≥ 50 % of hours;
  date-block permutation null of max |r| over all windows (`07_model_checks.R`).
- **Models.** Selection on the asinh scale (was raw flux), deterministic greedy selection (the 100 random
  orders were identical), sampling date as a random effect throughout, final models refitted on every
  row with their own predictors, Nakagawa marginal R², Satterthwaite tests.
- **New headline numbers.** Wetland core: TS_Ha2 69 h × water table 114 h × species, n = 619,
  R²m 63.0 % (species alone 33.6 %), leave-one-date-out R² 0.77 (species 0.44); temperature is 96 %
  seasonal, and season × WTD fits as well on shared rows. Upland A/B: 7.8 % / 7.1 % (species 3.7 / 3.5 %),
  no out-of-sample skill. Window permutation p: wetland TS < 0.002, WTD 0.018; upland TS 0.004;
  upland moisture n.s.
- Draft rebuilt: `DRAFT_ ... - v2 tracked 2026-10-01.docx` (168 tracked edits, 31 replies, 8 new comments).

## 1 Oct 2026, second update: driver data releases and analyzer volume

- **Driver data:** AmeriFlux US-Ha1 28-5 (to Apr 2026) and US-Ha2 17-5 (to May 2026); HF001 and HF070
  15-min tables to Sep 2026 (`data/raw/hf_archive/`); NEON HARV soil moisture, soil temperature, IR canopy
  temperature, throughfall, soil heat flux, wind and the eddy-covariance bundle, RELEASE-2026 to Jun 2025 and
  provisional to Dec 2025 (`data/raw/NEON_2026/`, release per month in `release_log.csv`). The xHA variables
  are rebuilt from NEON products (`01_import/07b_neon_xha.R`; QC-filtered fluxes; six-height profile means);
  agreement with AmeriFlux xHA over 2023–24 is r ≥ 0.98 except the CO2 profile composite (0.67) and wind (unused).
  `06_preprocess_soil_moisture.R` and `07b_neon_xha.R` run outside `run_pipeline.sh` (rerun only when NEON
  downloads change). `stackByTable` can delete the files it unpacks, so both stack a temporary copy.
- **Analyzer volume** 0.028 L for both analyzers (see item 4): MGGA fluxes and MDFs × 0.715.
- **Headline numbers:** wetland 1.67 ± 0.22 vs upland 0.04 ± 0.01 nmol m-2 s-1 (37-fold); N. sylvatica 4.87;
  wetland core model n = 912 to Oct 2025, R2m 57.6 % (species 31.3 %), leave-one-date-out R2 0.71;
  temperature explains flux beyond its seasonal cycle (chi2 = 22.0, p < 0.001); upland models 5.6–5.8 %.
- The LGR is the Microportable GGA (GLA131; raw files `micro_*.txt`), not the UGGA. Data column values still
  read "LGR/UGGA" because ch4-data-filtering imports them.
- Draft: `DRAFT_ ... - v2 tracked 2026-10-01b.docx` (170 tracked edits); EDI package rebuilt.
- **Model accuracy and bias** (`07_model_checks.R`, `outputs/tables/model_checks/accuracy_bias.csv`): on the
  asinh scale the models are unbiased (wetland leave-one-date-out RMSE 0.49, bias −0.03); back-transformed,
  the wetland model underestimates mean flux by about half (mean predicted/observed 0.52; 0.57 with Duan
  smearing) because the largest N. sylvatica emissions are underpredicted. Upland: 0.90 (0.99 smeared).
- **Analyzers:** LGR Microportable (field logs say LGR1; LGR3 files also exist on 8 dates, refits there match
  the original fluxes at r = 0.995) and LI-7810 TG10-01861 (2025). Units, dry correction, windows, clock,
  daily precision/MDF and volume are handled per analyzer; no side-by-side period, so an analyzer offset is
  confounded with 2025.

## 1 Oct 2026, SI and figures
- The whole SI is generated: `scripts/05_manuscript/01_build_si.R` → `outputs/manuscript/Supporting_Information.docx`
  (13 tables, 9 figures) and `si_items.csv`; runs last in `run_pipeline.sh`.
- New figures: `02_analysis/11_site_map.R` (site map from tree GPS, DEM, swamp outline, towers),
  `02_analysis/12_driver_timeseries_si.R` (daily drivers, replacing a static image). `05_combined_driver_timeseries.R`
  now writes `combined_flux_drivers_05.png` (it was overwritten by `04_interaction_plots.R`).
- VPD fix: `00_download/03_download_hf_met_hydro.R` passed pressure in Pa to `plantecophys::RHtoVPD`, which expects kPa;
  VPD was ~4× too high. Screening correlations unchanged (r 0.14 at the wetland).
- `06_manuscript_numbers.R` now gives Satterthwaite tests (refitting with lmerTest when conversion fails).
