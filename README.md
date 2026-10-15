# Tree CH4 Flux at Harvard Forest

Stem methane (CH4) fluxes of 60 trees in a forested wetland (Black Gum Swamp) and an upland forest (EMS tower footprint) at Harvard Forest, 2023–2025, with the environmental drivers, tree species and internal wood condition (tomography) that control them. This repository holds the code that turns the raw data into the published dataset, figures, tables, manuscript numbers and Supporting Information.

## How the pipeline is organized

```
data/package/    our primary data, exactly as published (read-only for the pipeline)
data/external/   public data from other providers (AmeriFlux, NEON, Harvard Forest archive, ...)
data/interim/    intermediate files (regenerated)
data/final/      the compiled datasets every analysis reads (regenerated)
outputs/         figures, tables, models, manuscript numbers, SI (regenerated)
```

Scripts run in numbered stages; within a stage, in numbered order.

| Stage | Folder | What it does | Main outputs |
|---|---|---|---|
| 0 | `scripts/0_data/` | get the data (not run by the pipeline) | `data/package/`, `data/external/` |
| 1 | `scripts/1_environment/` | met and water table, soil moisture, NEON tower variables → hourly drivers | `data/final/environment_hourly.csv` |
| 2 | `scripts/2_flux/` | raw analyzer records + field logs → fluxes, detection limits, QC flags | `data/final/stem_ch4_flux.csv` (+ dictionary, processing log, settings) |
| 3 | `scripts/3_trees/` | tomography decay classes; stand composition | `data/final/tomography_classes.csv`, Table 1 inputs |
| 4 | `scripts/4_analysis/` | descriptive analyses, figures and SI figures | `outputs/figures/`, `outputs/tables/` |
| 5 | `scripts/5_models/` | driver models, checks, every number quoted in the manuscript | `outputs/models/`, `outputs/tables/manuscript/` |
| 6 | `scripts/6_manuscript/` | Supporting Information document | `outputs/manuscript/Supporting_Information.docx` |
| 7 | `scripts/7_publish/` | assemble the data package and its EML for the Harvard Forest Data Archive (`00_attributes.R`, then `01_build_edi_package.R`; metadata text in `metadata/`) | `data/edi/package/` |

Run everything after stage 0:

```bash
bash scripts/run_pipeline.sh
```

`bash scripts/run_pipeline.sh 4 5` runs selected stages; `SKIP_FIT=1` skips steps 1–2 of stage 2 and reuses the closures and flux fits (the goFlux fit takes ~20 min). Logs go to `outputs/logs/`.

### Stage 0: getting the data

- **Our data** (`data/package/`): `0_data/02_download_data_package.R` downloads the published package (Harvard Forest Data Archive, knb-lter-hfr) and unpacks it.
- **Public data** (`data/external/`): `03_download_ameriflux.R` (US-Ha1, US-Ha2, US-xHA BASE), `04_download_neon.R` (NEON HARV, RELEASE-2026 + provisional; needs a NEON token in `NEON_TOKEN`), `05_download_phenocam.R`. The Harvard Forest Fisher met and hydrology tables (HF001, HF070) are downloaded by `1_environment/01_met_hydro.R` when they are not already in `data/external/hf_archive/`. GIS layers and the ForestGEO census: see `data/README.md`.

### Stage 2: from raw records to the flux dataset

Inputs: the raw analyzer records and `data/package/measurements.csv` (one row per measurement: tree, analyzer, chamber closure and opening times, deadband, sampling time, air temperature and pressure).

1. `01_closure_table.R` — each measurement's window (closure + deadband to the chamber opening, ended early at a detected chamber opening), its trace from the raw 1-Hz record, chamber volume, and the analyzer precision per analyzer and day (`goFlux::empirical.prec`).
2. `02_fit_fluxes.R` — goFlux linear and Hutchinson–Mosier fits, `best.flux` selection (Hüppi et al. 2018) and screens (`qc.flags`). MDF = 1.96 σ / t × flux term, σ = analyzer precision on that day.
3. `03_trace_qc.R` — automated checks of every CO2 and CH4 trace (flags carried into the dataset).
4. `05_deadband_sensitivity.R` — refits with 0–45 s deadbands (sensitivity check).
5. `06_flux_dataset.R` — final dataset: detection limits, sampling campaigns, QC flags; checked against `trees.csv`.

## Requirements

Working directory: the project root. R ≥ 4.3 with
`tidyverse` (`dplyr`, `tidyr`, `readr`, `stringr`, `tibble`, `ggplot2`), `lubridate`, `readxl`, `jsonlite`, `goFlux` (fork, version 0.5.0.9002: `remotes::install_github("jgewirtzman/goFlux@v0.5.0.9002")`, doi:10.5281/zenodo.23256675), `neonUtilities`, `plantecophys`, `lme4`, `lmerTest`, `performance`, `emmeans`, `RcppRoll`, `zoo`, `patchwork`, `cowplot`, `ggtext`, `scales`, `png`, `magick`, `pheatmap`, `car`, `sf`, `terra`, `ggspatial`, `ggnewscale`, `ggrepel`, `officer`, `flextable`; for stage 0/7 also `EDIutils`, `amerifluxr`, `phenocamr`, `EMLassemblyline`.

## Sites and species

| Code | Species | Site |
|---|---|---|
| bg | *Nyssa sylvatica* (black gum) | wetland |
| rm | *Acer rubrum* (red maple) | both |
| hem | *Tsuga canadensis* (eastern hemlock) | both |
| ro | *Quercus rubra* (red oak) | upland |

BGS = Black Gum Swamp (wetland); EMS = Environmental Measurement Station (upland). Ten trees per species and site.

## Citation and license

Code: MIT (`LICENSE`). Cite the archived release (Zenodo DOI; see `CITATION.cff`), the paper, and the data package in the Harvard Forest Data Archive.
