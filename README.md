# Tree CH4 Flux at Harvard Forest

## Overview

This project analyzes tree stem methane (CH4) emissions in upland and wetland forest ecosystems at Harvard Forest. It combines tree-level flux chamber measurements with environmental driver data (soil temperature, water table depth, soil moisture, phenology) to understand species-specific and environment-driven controls on CH4 emissions.

## Project Structure

```
tree-flux-2025/
├── scripts/
│   ├── 00_download/         # Programmatic data downloads
│   │   ├── 00_download_edi.R          # Study data from EDI → data/input/
│   │   ├── 01_download_ameriflux.R    # AmeriFlux towers (Ha1, Ha2, xHA)
│   │   ├── 02_download_phenocam.R     # PhenoCam GCC/NDVI
│   │   └── 03_download_hf_met_hydro.R # Fisher Met + water table → wtd_met.csv
│   ├── 01_import/           # Data preprocessing and alignment
│   │   ├── 01_tower_flux.R
│   │   ├── 02_tower_temperature.R
│   │   ├── 03_tower_moisture.R
│   │   ├── 04_neon_download.R
│   │   ├── 05_neon_moisture.R
│   │   ├── 06_preprocess_soil_moisture.R
│   │   ├── 07_phenocam.R
│   │   ├── 08_align.R               # Produces aligned_hourly_dataset.csv
│   │   ├── 09_goflux_reprocess.R    # Raw traces → goFlux/fluxqc → HF_2023-2025_tree_flux_goflux.csv
│   │   └── 10_quality_flags.R       # Analysis dataset with MDF flags → flux_with_quality_flags.csv
│   ├── 02_analysis/         # Core analyses and figures
│   │   ├── 00_data_summary.R         # QC/MDF summary statistics
│   │   ├── 01_timeseries.R           # Temporal flux plots by species
│   │   ├── 02_flux_summaries.R       # Main boxplot figure + mixed models
│   │   ├── 03_repeatability.R        # ICC, Spearman, z-score tracks
│   │   ├── 04_rolling_corrs.R        # Rolling-window correlations
│   │   ├── 05_combined_driver_timeseries.R
│   │   ├── 06_filtering_snr.R        # QC visualization
│   │   ├── 07_filter_sensitivity_ridges.R  # MDF filter sensitivity analysis
│   │   └── 08_tomography.R           # ERT/Sonic imaging + flux
│   ├── 03_modeling/         # Statistical models
│   │   ├── 01_bgs_model.R            # Wetland mixed-effects model
│   │   ├── 02_ems_model_A.R          # Upland instantaneous drivers
│   │   ├── 03_ems_model_B.R          # Upland BGS-style drivers
│   │   ├── 04_interaction_plots.R
│   │   └── 05_compare_models.R
│   └── helpers/             # Shared utilities
│       └── find_ameriflux.R          # Version-agnostic AmeriFlux path lookup
├── data/                    # All data gitignored (see data/README.md)
│   ├── raw/                 # Source data downloads
│   ├── input/               # Study data (from EDI package)
│   └── processed/           # Script-generated intermediates
├── outputs/                 # All gitignored
│   ├── figures/
│   ├── tables/
│   └── models/              # Saved .rds model objects
├── archive/                 # Legacy scripts for reference
├── .gitignore
├── tree-flux-2025.Rproj
└── README.md
```

## Workflow

Scripts are numbered to indicate execution order. Run them sequentially within each phase:

### Phase 0: Data Download (`scripts/00_download/`)
Programmatically downloads all data sources. Run in order:

| Script | Source | Requires |
|--------|--------|----------|
| `00_download_edi.R` | Study data (flux, tomography) from [EDI](https://portal.edirepository.org/) | `EDIutils` |
| `01_download_ameriflux.R` | AmeriFlux towers (Ha1, Ha2, xHA) | `amerifluxr` + free account ([register here](https://ameriflux-data.lbl.gov/Pages/RequestAccount.aspx)) |
| `02_download_phenocam.R` | PhenoCam (harvardems2) | `phenocamr` |
| `03_download_hf_met_hydro.R` | Harvard Forest LTER (Fisher Met + hydro) | `plantecophys` (downloads from [EDI/PASTA](https://pasta.lternet.edu/)) |

NEON data is downloaded directly within the import scripts (`04_neon_download.R`, `05_neon_moisture.R`, `06_preprocess_soil_moisture.R`) via `neonUtilities::loadByProduct()`.

### Phase 1: Data Import (`scripts/01_import/`)
Preprocesses raw downloads and aligns everything into a single hourly dataset. `08_align.R` produces `data/processed/aligned_hourly_dataset.csv`.

`09_goflux_reprocess.R` refits every stem-chamber closure that has a raw 1-Hz trace through [goFlux](https://github.com/Qepanna/goFlux) and [fluxqc](https://github.com/jongewirtzman/fluxqc) (>= 0.2.3), so that the published dataset and the papers built on it share one flux fit:

- windows come from the field-log start/end times (LGR/UGGA, Jun 2023 - Mar 2025; the Sep 2024 onward logs are the xlsx files in `data/raw/upland_wetland/Sept2024_onwards/`) or from the analyzer REMARK span with a 20-s deadband (LI-7810, Apr 2025 onward); nothing is clicked by hand;
- flux = `goFlux::best.flux` (linear or Hutchinson-Mosier, Hüppi et al. 2018 criteria);
- detection: MDF = 1.96 σ / t · flux.term, σ = MAD of first differences over each analyzer's whole record per constant-interval run (`fluxqc::flag_detection(precision = "mad")`), t = closure length in seconds; retain-and-flag;
- chamber volume = measured collar volume (`tree_volumes.csv`) + analyzer loop (0.200 L LGR/UGGA, 0.028 L LI-7810, the same conventions as the legacy pipelines); collar area = π·5.08² cm²;
- physical QC screens (`fluxqc::qc_screens`: C0, CO2 tracer, convexity, window length, noisy closure) only flag; the flagged closures are listed in `outputs/tables/goflux_qc_review.csv` for review.

It writes `data/input/HF_2023-2025_tree_flux_goflux.csv` (one row per closure, legacy and goFlux fluxes and SEs side by side, `flux_source` says which is canonical), trace plots in `outputs/figures/goflux/`, and comparison tables in `outputs/tables/goflux_*.csv`. Raw traces are read from `Matthes_Lab/stem-CH4-flux/Raw LGR Data` (LGR day files) and `data/raw/7810_Processed/Tree_Fluxes/raw-7810-tree_flux_data` (+ the Aug 2025 files in `Matthes_Lab/stem-CH4-flux/raw-7810-data`); the Oct 2025 LI-7810 files are not archived, so those closures keep their legacy flux.

`10_quality_flags.R` builds the analysis dataset `data/processed/flux_with_quality_flags.csv` from it: the canonical flux (`CH4_flux_nmolpm2ps`, goFlux where a trace exists, legacy otherwise), the reference detection flag `CH4_below_MDF` (campaign-σ 95 %), and the comparison MDFs used by the sensitivity analysis (datasheet, campaign σ at 90/95/99 %, per-closure σ ×3·t). By default only closures that were in the legacy dataset are promoted to the analysis dataset (n = 1,640); set `GOFLUX_INCLUDE_NEW=1` to include the closures with traces that the legacy pipeline dropped. All downstream scripts load this flagged dataset.

### Phase 2: Analysis (`scripts/02_analysis/`)
Runs independently once Phase 1 is complete. Produces publication figures and summary statistics.

### Phase 3: Modeling (`scripts/03_modeling/`)
Depends on Phase 2 output (particularly `04_rolling_corrs.R` for optimal window sizes). Fits mixed-effects models for wetland and upland sites.

## Data

All data files are gitignored. See `data/README.md` for sources and instructions on obtaining the raw data.

**Key input files:**
- `data/input/HF_2023-2025_tree_flux_corrected.csv` -- Legacy tree-level CH4 flux measurements (upstream linear fits; the 2023-24 `CH4_SE` in this file is SE_slope × n in ppm without the /area division and is corrected in `09_goflux_reprocess.R`)
- `data/input/HF_2023-2025_tree_flux_goflux.csv` -- goFlux-reprocessed closures with legacy values side by side (generated by `09_goflux_reprocess.R`; this is the file to publish and to import into dependent projects)
- `data/processed/flux_with_quality_flags.csv` -- Analysis dataset with MDF flags, precision, and SNR (generated by `10_quality_flags.R`)
- `data/processed/aligned_hourly_dataset.csv` -- Hourly environmental drivers (generated by Phase 1)

## Requirements

All scripts assume the working directory is the project root (`tree-flux-2025/`).

R packages:

**Download:** `EDIutils`, `amerifluxr`, `phenocamr`, `plantecophys`, `neonUtilities`

**Flux processing:** `goFlux` (>= 0.4.0), `fluxqc` (>= 0.2.3, `remotes::install_github("jongewirtzman/fluxqc")`), `readxl`

**Analysis:** `tidyverse`, `lubridate`, `lme4`, `emmeans`, `performance`, `patchwork`, `zoo`, `RcppRoll`, `ggtext`, `scales`, `magick`, `ggpointdensity`, `ggridges`, `viridis`, `car`, `cowplot`, `pheatmap`, `readxl`

## Sites

- **BGS** (Black Gum Swamp) -- Wetland site
- **EMS** (Environmental Measurement Station) -- Upland site

## Species

| Code | Species | Common Name | Strategy |
|------|---------|-------------|----------|
| bg | *Nyssa sylvatica* | Black gum | Wetland specialist |
| rm | *Acer rubrum* | Red maple | Generalist |
| hem | *Tsuga canadensis* | Eastern hemlock | Generalist |
| ro | *Quercus rubra* | Red oak | Upland specialist |
