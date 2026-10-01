# Repo restructure — status at pause (1 Oct 2026)

Not committed yet; the working tree holds the restructure in progress.

## Done
- New data layout:
  - `data/package/`: our primary data. `0_data/01_assemble_data_package.R` built it: 452 LGR files, 16 LI-7810 files, 12 field logs, chamber volumes, `trees.csv` (60 trees), the earlier dataset, manual windows, tomography and stand data.
  - `data/external/`: public data (AmeriFlux, NEON, HF archive, PhenoCam, GIS, ForestGEO).
  - `data/interim/` and `data/final/`: rebuilt by the pipeline.
  - Old `data/raw`, `data/input` and `data/processed` were moved to `data/_old/`; nothing was deleted.
- Scripts reorganized into stages `0_data` … `7_publish` with `git mv`. Unused scripts were moved to `archive/scripts_2026-10/`. All paths and cross-references were updated.
- Flux stage:
  - The old goFlux script is split into `01_closure_table` and `02_fit_fluxes`, which share `flux_settings.R`.
  - Every cleaning rule now records its count via `log_step()`.
  - `06_flux_dataset` writes `data/final/stem_ch4_flux.csv`, the data dictionary and `flux_processing_log.csv`, and joins the trace-QC flags.
- Rewritten: `run_pipeline.sh` (stages), `README.md`, `data/README.md`, `4_analysis/01_data_summary.R`, `0_data/02_download_data_package.R`, `7_publish/00_attributes.R` and `01_build_edi_package.R`.
- Stage 1 rerun: `environment_hourly.csv` is identical to the baseline.

## Next
1. Check stage 2. The background run of `bash scripts/run_pipeline.sh 2` was started at pause; see `outputs/logs/STATUS.txt`.
   - Compare `data/final/stem_ch4_flux.csv` with the baseline in the scratchpad: 1,637 rows, same fluxes and MDF.
   - Read `flux_processing_log.csv`.
2. Run stages 3–6 and compare the tables, models and `manuscript_numbers.txt` with the baseline.
3. Rebuild the draft (manuscript_work). Reply 39 and the script paths in the replies were updated.
4. Commit. Track `7_publish`? (It was untracked by choice in March.)
5. Delete `data/_old/` after checking.

## Stage 2 result (after the pause)
- 01–05 ran without errors; the fit took 381 s.
- 06_flux_dataset stopped at its tree-metadata check: PLOT disagrees with trees.csv. The earlier dataset codes some upland trees by sub-plot (e.g. "E5" for tree 903), while trees.csv uses "EMS". Decide which PLOT coding the final dataset uses (likely keep the earlier sub-plot codes in the flux table and compare site + species only), then rerun stage 2 with `SKIP_FIT=1`.

## Verified (2 Oct 2026)
- Stage 2: `stem_ch4_flux.csv` is identical to the pre-restructure analysis table on every key column (1,637 rows: fluxes, SE, MDF, detection, sampling hour, tree, species). The processing log has 35 rules.
- Stages 3–6:
  - All model, check and manuscript tables agree to within 2×10⁻⁸ relative (optimizer noise).
  - `manuscript_numbers.txt` is identical.
  - Figures are pixel-identical, except Figures 1–2, whose jitter is now seeded (`set.seed(1)`) so reruns are byte-stable.
- Draft rebuilt and validated.
