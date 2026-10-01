# Data

Nothing in `data/` is tracked by git. Four folders, in the order data flow through them:

| Folder | Contents | Made by |
|---|---|---|
| `package/` | the study's own primary data, as recorded; published as the data package | `scripts/0_data/02_download_data_package.R` (or `01_assemble_data_package.R`, authors) |
| `external/` | public data from other providers | `scripts/0_data/03-05`, `1_environment/01_met_hydro.R`, by hand (below) |
| `interim/` | intermediate files | stages 1–2 |
| `final/` | compiled datasets read by every analysis | stages 1–3 |

`_old/` holds the pre-October-2026 folders (`raw/`, `input/`, `processed/`) until the new layout has been checked; it can then be deleted. `edi/` holds EML templates and the assembled package (`scripts/7_publish/`).

## `package/` — primary data (published)

| Path | Description |
|---|---|
| `trees.csv` | the 60 study trees: ForestGEO tag, site, plot, species, DBH, microtopography (wetland), GPS position (BGS survey 2025-04-17, EMS 2026-01-16) |
| `chamber_volumes.csv` | measured collar depths and collar + cap + tubing volume per tree; old wetland tags (1–32) |
| `field_logs/` | field logs as recorded: closure start/end (analyzer clock), real time, analyzer, notes. `Field_Data_Monthly_Summer2023.csv`, `Field_Data_Monthly_updated2.csv` (2023–Jan 2024, with the team's checks), `summer_2024.csv`, `Timing Updates.xlsx` (team's hand-adjusted windows, May–Jul 2024), `Stem_flux_*.xlsx` (Sep 2024–Mar 2025) |
| `analyzer_raw/lgr/LGR{1,2,3}/` | LGR microportable GLA131 (UGGA) 1-Hz day files, May 2023–Apr 2025 |
| `analyzer_raw/li7810/` | LI-COR LI-7810 (SN TG10-01861) 1-Hz day files, 2025; closures tagged in the REMARK field |
| `previous_processing/HF_2023-2025_tree_flux_v1.csv` | the earlier processed dataset (1,640 measurements). Used for the curated tree tags, the met at each measurement and the fluxes of measurements without an archived 1-Hz record (October 2025, one day in May 2024, six closures in September 2023) |
| `qc_decisions/manual_windows.csv` | windows set or excluded after inspecting the concentration record (`scripts/2_flux/04_review_windows.R`) |
| `tomography/` | `tomography_results_compiled.csv` (ERT and sonic metrics), `ERT_application_results.csv` and `Tree_ID_info.csv` (companion tomography study, Thompson et al. 2026), `images/` (ERT and sonic cross-sections) |
| `stand/` | `BGS_VRP_2025.csv` (June 2025 prism survey of the swamp, BAF 10 ft²/acre, 40 plots), `Black_Gum_Swamp.kmz` (swamp outline) |

## `external/` — public data

| Path | Source | How to get it |
|---|---|---|
| `ameriflux/` | AmeriFlux BASE: US-Ha1 (v28-5), US-Ha2 (v17-5), US-xHA (v11-5) | `0_data/03_download_ameriflux.R`, or download from ameriflux.lbl.gov |
| `neon/` | NEON HARV, RELEASE-2026 + provisional: DP1.00001, 00005, 00040, 00041, 00046, 00094, DP4.00200 | `0_data/04_download_neon.R` (`NEON_TOKEN`) |
| `hf_archive/` | Harvard Forest Data Archive: HF001 table hf001-10 (Fisher met, 15 min), HF070 table hf070-04 (hydrology incl. BGS/BVS wells, 15 min) | `1_environment/01_met_hydro.R` downloads the newest revision if absent |
| `phenocam/` | PhenoCam harvardems2 NDVI | `0_data/05_download_phenocam.R` |
| `gis/` | Harvard Forest GIS (HF000): `elevation_ned/` (30-m DEM, MA State Plane), `elevation_contours.geojson`, `tracts.*`, `stands_1986_93.*`; US Census `cb_2018_us_state_500k.zip` | by hand from the Harvard Forest GIS archive and census.gov |
| `forestgeo/hf253-06-stems-2019.csv` | Harvard Forest ForestGEO census, 2019 (HF253 v6) | by hand from the Harvard Forest Data Archive (HF253) |

## `interim/` — intermediate files

`wtd_met.csv` (hourly water table + Fisher met), `neon_swc_hourly.csv`, `tower_swc_ts_hourly.csv`, `neon_sensor_depths.csv`, `neon_xha_hourly.csv`, `neon_ec_extract.rds` (cache of the NEON eddy-covariance extract), `flux_closures.rds`, `flux_fits.csv` (every closure, including those outside the study design), `flux_traces.rds`, `flux_trace_qc.csv`, `flux_log_*.csv`.

## `final/` — datasets read by the analyses

| File | Description |
|---|---|
| `stem_ch4_flux.csv` | one row per measurement (1,637): flux, SE, MDF, detection class, QC flags, tree, analyzer, window, chamber, corrected sampling time |
| `stem_ch4_flux_dictionary.csv` | column definitions and units |
| `flux_processing_log.csv` | every cleaning rule and the number of records it touched |
| `flux_processing_settings.json` | goFlux/fluxqc settings and constants |
| `environment_hourly.csv` (+ `_variables.csv`) | hourly drivers on the EST clock (UTC−5) |
| `tomography_classes.csv` | decay class and ERT PC1 per tree (reproduces the companion study) |
