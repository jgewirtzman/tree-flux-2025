# Switch from fluxqc to the goFlux fork (8 Oct 2026)

The flux scripts (`scripts/2_flux/`) now use only the goFlux fork, release v0.5.0.9001. fluxqc will not be released.

- Install: `remotes::install_github("jgewirtzman/goFlux@v0.5.0.9001")`.
- Zenodo: doi:10.5281/zenodo.23254791.
- Paper citation: "goFlux (Rheault et al. 2024; version 0.5.0.9001 with additions, Gewirtzman 2026a)".

## What changed

| Script | fluxqc | goFlux fork |
|---|---|---|
| 01_closure_table.R | `precision_mad_runs` (whole-day record, first differences / √2) | `empirical.prec(method = "hadamard", by = analyzer × day)` on the closure windows; same `day_sigma` columns |
| 02_fit_fluxes.R | `process_fluxes`, `closure_seconds`, `to_json` | `process.fluxes`, `closure.time`, `jsonlite::toJSON` |
| 03_trace_qc.R | `precision_mad`, `find_rise` | local `precision_mad` helper (same formula), `find.rise` |
| 04_review_windows.R | `click_peak2_stacked` | `click.peak2(gases = ...)` |
| 06_flux_dataset.R | `flux_term` | `flux.term` |

## Kept as local code

- **CO₂ screen.** `co2.tracer` is TRUE when CO₂ rises (pass). The dataset column `qc_co2_tracer` is its negation (TRUE = flagged), as before.
- **Convex screen.** It stays in `02_fit_fluxes.R`, as in fluxqc. A quadratic in time is fitted to the CH₄ window, and the closure is flagged when the quadratic term has the sign of the net trend and p < 0.05.

## Removed

- fluxqc's comparison detection-limit columns (`extra`), which nothing downstream used.

## Effect

- **Fluxes are unchanged.** The maximum difference is 4 × 10⁻¹⁰ nmol m⁻² s⁻¹, with the same model choices and closure lengths. The fit uses the unchanged `*_prec` columns.
- **Noise estimates and detection limits are lower:**

  | Analyzer | σ (ppb), old → new | Change | Typical MDF (nmol m⁻² s⁻¹), old → new |
  |---|---|---|---|
  | LGR | 3.84 → 3.53 | −7% | 0.117 → 0.107 |
  | LI-7810 | 0.135 → 0.098 | 1.43× lower | 0.0037 → 0.0025 |

- **Fluxes smaller in magnitude than the MDF:** 62% → 60% upland, 27% → 25% wetland. The upland analyzer comparison holds; the numbers in the Discussion were updated.
- **Screen counts:**

  | Screen | fluxqc | goFlux fork |
  |---|---|---|
  | CO₂ | 44 | 21 |
  | Noisy | 76 | 102 |
  | Convex (local code) | 203 | 203 |
  | c0 | 13 | 13 |

  The fork's CO₂ and noisy screens are defined slightly differently. Flags only; no data are removed.

## Note for the goFlux fork

`qc.flags()` `qc.convex` tests only `HM.k < 0`. Under the default `k.min` it never fires: it flagged 0 of 1,846 closures, where fluxqc's data-based test flagged 203.

A data-based branch, used when the flagged dataframe is supplied, would restore it, keeping `HM.k < 0` as the fallback. The branch: `lm(conc ~ Etime + I(Etime^2))` inside the window; flag when the sign of the quadratic term equals the sign of the net slope and p < 0.05; at least 6 points.

`process.fluxes(prec = NULL)` also warned on the LI-7810 group about the noise (lag-1 autocorrelation of second differences −0.42; d1c ratio 1.31).
