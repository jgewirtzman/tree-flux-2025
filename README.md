# Tree CH4 Flux at Harvard Forest

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23421659.svg)](https://doi.org/10.5281/zenodo.23421659)

Stem methane (CH4) fluxes of 60 trees in a forested wetland (Black Gum Swamp) and an upland forest (EMS tower footprint) at Harvard Forest, 2023–2025, with the environmental drivers, tree species and internal wood condition (tomography) that control them. This repository holds the code that turns the raw data into the published dataset, figures, tables, manuscript numbers and Supporting Information.

## How the pipeline is organized

```
data/package/    our primary data, exactly as published (read-only for the pipeline)
data/external/   public data from other providers (AmeriFlux, NEON, Harvard Forest archive, ...)
data/interim/    intermediate files (regenerated)
data/final/      the compiled datasets every analysis reads (regenerated)
outputs/         figures, tables, models, manuscript numbers, SI (regenerated)
```

Scripts run in numbered stages; within a stage, in numbered order. Shared functions are in
`scripts/helpers/` (rolling windows and anomalies, AmeriFlux file lookup) and
`scripts/5_models/model_helpers.R` (predictor screening and forward selection for the three driver
models). `scripts/exploratory/` holds analyses that are not part of the paper and are not run by the pipeline.

| Stage | Folder | What it does | Main outputs |
|---|---|---|---|
| 0 | `scripts/0_data/` | get the data (not run by the pipeline) | `data/package/`, `data/external/` |
| 1 | `scripts/1_environment/` | met and water table, soil moisture, NEON tower variables → hourly drivers | `data/final/environment_hourly.csv` |
| 2 | `scripts/2_flux/` | raw analyzer records + measurement times → fluxes, detection limits, QC flags | `data/final/stem_ch4_flux.csv` (+ dictionary, settings) |
| 3 | `scripts/3_trees/` | tomography decay classes; stand composition | `data/final/tomography_classes.csv`, Table 1 inputs |
| 4 | `scripts/4_analysis/` | descriptive analyses, figures and SI figures | `outputs/figures/`, `outputs/tables/` |
| 5 | `scripts/5_models/` | driver models, model checks, interaction figure, manuscript numbers | `outputs/models/`, `outputs/tables/manuscript/` |
| 6 | `scripts/6_manuscript/` | Supporting Information document | `outputs/manuscript/Supporting_Information.docx` |
| 7 | `scripts/7_publish/` | assemble the data package and its EML for the Harvard Forest Data Archive (`00_attributes.R`, then `01_build_edi_package.R`; metadata text in `metadata/`) | `data/edi/package/` |

Run everything after stage 0:

```bash
bash scripts/run_pipeline.sh
```

`DRIVERS=published bash scripts/run_pipeline.sh` skips rebuilding the environmental drivers and uses the hourly drivers archived in the data package (`environment_hourly.csv`), so stages 2–6 need only the data package plus the ForestGEO census and GIS layers (stand table and site map; see `data/README.md`). This also avoids the public sources that are updated in place (AmeriFlux, NEON, Harvard Forest archive, PhenoCam); the versions behind the published drivers are listed in `data/README.md`; `bash scripts/run_pipeline.sh 4 5` runs selected stages; `SKIP_FIT=1` skips steps 1–2 of stage 2 and reuses the closures and flux fits (the goFlux fit takes ~6 min; a full run takes ~25 min). Logs go to `outputs/logs/`.

### Stage 0: getting the data

- **Our data** (`data/package/`): `0_data/01_download_data_package.R` downloads the two published datasets (stem flux; tomography) from the Harvard Forest Data Archive (knb-lter-hfr) and unpacks them. The datasets are under review at the Harvard Forest Data Archive; until their IDs are issued, they are available from the corresponding author.
- **Public data** (`data/external/`): `02_download_ameriflux.R` (US-Ha1, US-Ha2, US-xHA BASE), `03_download_neon.R` (NEON HARV, RELEASE-2026 + provisional; a NEON API token in `NEON_TOKEN` is optional and speeds up the download), `04_download_phenocam.R` (PhenoCam harvardems2 NDVI summary). The Harvard Forest Fisher met and hydrology tables (HF001, HF070) are downloaded by `1_environment/01_met_hydro.R` when they are not already in `data/external/hf_archive/`. GIS layers and the ForestGEO census: see `data/README.md`.

### Stage 2: from raw records to the flux dataset

Inputs: the raw analyzer records and `data/package/measurements.csv` (one row per measurement: tree, analyzer, chamber closure and opening times, deadband, sampling time, air temperature and pressure).

1. `01_closure_table.R` — each measurement's window (closure + deadband to the chamber opening, ended early at a detected chamber opening), its trace from the raw 1-Hz record, chamber volume, and the analyzer precision per analyzer and day (`goFlux::empirical.prec`).
2. `02_fit_fluxes.R` — goFlux linear and Hutchinson–Mosier fits, `best.flux` selection (Hüppi et al. 2018) and screens (`qc.flags`). MDF = 1.96 σ / t × flux term, σ = analyzer precision on that day.
3. `03_trace_qc.R` — automated checks of every CO2 and CH4 trace (flags carried into the dataset).
4. `04_deadband_sensitivity.R` — refits with 0–45 s deadbands (sensitivity check).
5. `05_flux_dataset.R` — final dataset: detection limits, sampling campaigns, QC flags; checked against `trees.csv`.

## Requirements

Working directory: the project root. The results were produced with R 4.4.3 and these package versions:

- **Flux calculation:** `goFlux` 0.5.0.9002 (fork: `remotes::install_github("jgewirtzman/goFlux@v0.5.0.9002")`, doi:10.5281/zenodo.23256675), `jsonlite` 2.0.0
- **Data handling:** `tidyverse` 2.0.0 (`dplyr` 1.2.0, `tidyr` 1.3.2, `readr` 2.2.0, `stringr` 1.6.0, `tibble` 3.3.1, `ggplot2` 4.0.2), `lubridate` 1.9.5, `zoo` 1.8.15, `RcppRoll` 0.3.1
- **Environmental data:** `neonUtilities` 3.0.2, `rhdf5` 2.50.2 (Bioconductor, `BiocManager::install("rhdf5")`; used by `neonUtilities::stackEddy` for the NEON eddy-covariance bundle), `plantecophys` 1.4.6
- **Models:** `lme4` 1.1-35.5, `lmerTest` 3.1-3, `pbkrtest` 0.5.3 (Kenward–Roger tests of the tree-level wood-condition metrics), `performance` 0.12.4, `emmeans` 1.10.6, `car` 3.1-3
- **Figures:** `patchwork` 1.3.2, `cowplot` 1.1.3, `ggtext` 0.1.2, `scales` 1.4.0, `png` 0.1.8, `magick` 2.9.0, `pheatmap` 1.0.12, `ggnewscale` 0.5.2, `ggrepel` 0.9.6
- **Maps:** `sf` 1.0-24, `terra` 1.8-60, `ggspatial` 1.1.9
- **Supporting Information:** `officer` 0.6.7, `flextable` 0.9.7
- **Stage 0 and 7 only:** `EDIutils` 2.0.0, `EMLassemblyline` 3.6.0, `amerifluxr`

## Sites and species

| Code | Species | Site |
|---|---|---|
| bg | *Nyssa sylvatica* (black gum) | wetland |
| rm | *Acer rubrum* (red maple) | both |
| hem | *Tsuga canadensis* (eastern hemlock) | both |
| ro | *Quercus rubra* (red oak) | upland |

BGS = Black Gum Swamp (wetland); EMS = Environmental Measurement Station (upland). Ten trees per species and site.

## Citation and license

Code: MIT (`LICENSE`); the Figure 7 illustration in `figures_static/` is not covered (see `figures_static/README.md`). Cite the archived release (Gewirtzman 2026, Zenodo, https://doi.org/10.5281/zenodo.23421659 for the latest version; v1.0.0: https://doi.org/10.5281/zenodo.23421660), the paper, and the two datasets in the Harvard Forest Data Archive (stem flux; tomography).
