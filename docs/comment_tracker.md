# Co-author comments on the working draft

Working draft: `DRAFT_ Contrasting controls on tree methane emissions in upland and wetland forests.docx`
(copied from Downloads on 30 Sep 2026; its text matches bioRxiv v1, 10.64898/2026.01.29.702553).
41 comments (Jan–Apr 2026) from Bradford, Marra, Jurado, Thompson, Lutz, Matthes and Gewirtzman.

Status key: **plan** = resolved by the existing workplan; **code** = needs a new analysis or figure change;
**text** = writing only; **done-local** = already coded but not yet committed.

**Update 30 Sep 2026:** every item below has been addressed as tracked changes in
`DRAFT_ Contrasting controls on tree methane emissions - v2 tracked 2026-09-30.docx`, with a threaded
reply (author "Claude (for J. Gewirtzman)") on each co-author comment. Numbers come from
`scripts/5_models/07_manuscript_numbers.R` (outputs/tables/manuscript/). Still open for Jon:
precision estimator vs the guidelines paper; Gelman 2-SD scaling (kept 1 SD); citation details for the
two in-prep companion papers; reference-list entries for Segers 1998 and Christensen et al. 2003 (and other
pre-existing gaps: Terazawa, Dunfield, Jevon, Leung, Marra 2018 etc. should be checked); the SI panel
combining tree trajectories with random effects (no current script); EDI vs HF archive.

**Update 1 Oct 2026:** current draft `DRAFT_ ... - v2 tracked 2026-10-01b.docx`. All 41 comments are covered:
31 thread-starting comments have a reply; the other 10 are follow-ups inside answered threads (8–11 repeat 7 on
other Table 2 rows; 1–2 rephrase 0; 20, 22, 26, 37 are agreements or follow-ups to 19/16, 21, 25, 36). Replies cite
SI items by key, resolved from `outputs/manuscript/si_items.csv`, so they stay correct when the SI is renumbered.
Since 30 Sep, resolved: reference entries (Segers, Christensen, Terazawa, Jeffrey 2021, Jevon 2023, Leung 2026,
Jenkins 2008; soils sentence now Allen 1995 and Davidson et al. 1998); trajectories/random-effects SI figure
(regenerated); Table 2 merged into SI Table S1 with class names in the Methods. Still open: Gelman 2-SD scaling
(kept 1 SD); companion-paper citations ("in prep."/preprint); Harvard Forest archive vs stand-alone EDI package.

## Resolved by the reprocessing / method-alignment workplan

| # | Who | Comment | Resolution |
|---|---|---|---|
| 3 | Gewirtzman | Need to add LI-COR methods | plan: Methods rewrite (two analyzers, goFlux, windows) |
| 4 | Lutz | 7810 fluxes use the "best" minute, not the full closure | plan: goFlux refit now uses the full remark minus a 20-s deadband for the LI-7810; state it |
| 36, 37 | Bradford | Variance partitioning sums to ~120 % | plan: those numbers are not reproducible; recomputed splits sum to 100 % (wetland 27.5/5.6/66.8, upland 1.2/8.6/90.2 on reprocessed data) |
| 7–11 | Jurado | Numeric ranges for the decay categories? | plan: adopt the tomography paper's classes (SoT structural loss > 1 %, species-normalized ERT PC1 > study-set mean); state thresholds, put per-tree values in SI |
| 5 | Jurado | ERT app needs more detail / uncertainty | plan: cite the tomography paper (workflow figure, hemlock core validation r = 0.68) instead of re-describing |
| 6, 34 | Marra | "damage" is an arboriculture term | plan: use "structural loss", as the tomography paper does |
| 33 | Thompson | SoT vs "acoustic tomography" used interchangeably | text: use SoT throughout, matching the tomography paper |
| — | (orphan paragraph, end of SI) | QC paragraph describes R² < 0.8, SNR < 2, CH4 < −1 exclusions | plan: delete; replace with retain-and-flag detection text (guidelines paper, fluxqc) |

## Need a new analysis, table or figure change

| # | Who | Comment | Action |
|---|---|---|---|
| 18 | Matthes | Add y = 0 line to Fig. 2 | done-local: uncommitted edit in `02_flux_summaries.R` |
| 39 | Bradford | Show ERT plots with all fluxes, not just tree means | done-local: uncommitted all-observations plots in `08_tomography.R`; decide main vs SI |
| 16, 19, 20 | Bradford, Marra | Give n behind the means | code: add n (obs, trees) to site/species/month summaries |
| 17 | Bradford | CIs on means? | code: add 95 % CI (tree-level bootstrap) to summaries |
| 27, 40 | Bradford | Supplementary table of all coefficients with SEs; show SEs in species-effects table | code: SI table from `outputs/models/*/coefficients.csv`; species slopes with SE (delta method or emtrends) |
| 30 | Bradford | Do VIFs rise in the full model (sign reversals)? | code/text: yes — reprocessed full model GVIF for TS_Ha2 = 84, s10t = 103 (core ≤ 9.5); report and use to justify core model |
| 15 | Bradford | State VIF thresholds | text + code check: state threshold used in `01_bgs_model.R` refinement step |
| 14 | Bradford | How standardized? Gelman 2008 for binary vs continuous | text: mean-centred / 1 SD; consider 2-SD scaling (Gelman 2008) — decision needed |
| 31 | Bradford | Give n and dates for the wetland model | text: N = 495, 30 trees, dates from model output |
| 35 | Bradford | What about maple in the wetland (ERT)? | text: A. rubrum wetland r = −0.31, p = 0.38 (reprocessed) |

## Writing only

| # | Who | Comment | Action |
|---|---|---|---|
| 0, 1, 2 | Jurado, Thompson, Marra | "typical sizes rather than external indicators of decay" unclear | rephrase: trees chosen by DBH quartiles without regard to external decay signs, so internal decay is not selected for |
| 12 | Gewirtzman | ProCheck soil readings not used | delete sentence |
| 13 | Jurado | Cite the choice of temperature + hydrology as core predictors | add citations |
| 21, 22 | Marra | Italicize species names | global fix |
| 23, 24, 25, 26, 32 | Bradford, Marra | ICC and z-score wording unclear | define z-score; say ICC = between-tree share of variance within a species × site; say what high z-rank consistency means |
| 28 | Bradford | "variation" — in space or time? | clarify |
| 29 | Bradford | "primary drivers" → drivers with most precision/relevance (less attenuation) | soften wording |
| 38 | Thompson | Hook et al. 1971 and Keeley & Franz 1979 not in reference list | add references |

## 1 Oct 2026 — SI trimmed

- SI reduced to 2 texts, 7 tables and 5 figures, all generated by `scripts/05_manuscript/01_build_si.R`.
- Decay class key merged into Table S1, credited to the companion tomography study.
- Per-tree decay metrics are cited to Thompson et al. 2026; ours reproduce them exactly.
- The decay robustness table (S2) now covers each species and site, plus pooled groups, against each metric, giving r (p), ρ and the leave-one-out range.
- Former tables S10 and S11 merged into the diagnostics table (S7). The full-coefficients table (S6) is kept to answer comment 27.
- The overlapping model-check table and figures cut down to the diagnostics table.
- All text comparing original and recalculated fluxes removed.
- Site map (Figure S1) redrawn at publication quality: PDF plus 600-dpi PNG.

## 1 Oct 2026 — Figure 3 kept, reframed around synchrony

- Figure 3 (z-score tracks) stays in the main text.
- Its caption and the repeatability Results now describe how much within-tree variation is shared across trees on the same date. This share comes from adding a sampling-date random effect: *N. sylvatica* 74%, wetland *A. rubrum* 36%, wetland *T. canadensis* 7%, upland species 2–8%. Source: `outputs/tables/synchrony_by_group.csv`, written by `4_analysis/04_repeatability.R`.
- Figure S3 is now the tree random effects only. The trajectory panel was dropped because it duplicated Figure 3.
- JGR-B Key Points and Plain Language Summary added to the draft.
- All repeatability statistics are now on the asinh scale, the same scale as the models, the variance partitioning and the synchrony shares. This covers the ICC, its LRT, the period Spearman correlations, the Figure 3 z-scores and the Figure S3 tree effects.
- Resulting values: ICC 0–0.28 (mean 0.11), with *N. sylvatica* 0.28 and wetland *A. rubrum* 0.20; ρ 0.32–0.90.

## 2 Oct 2026 — Concision and SI moves

- Wordiness and redundancy pass: 30 tracked edits.
- Moved to the SI:
  - Text S3 (new): tomography procedures.
  - Text S1: the flux equation.
  - Text S2: model selection, diagnostics, upland A vs B, variance partitioning.
- Software and data availability moved to a new Open Research section.
- Knowledge-gaps paragraph tightened, with added references: Megonigal et al. 2020, Jeffrey et al. 2021b, Plain & Epron 2021, Iorliam et al. 2026. The bark-methanotroph paper is now Jeffrey et al. 2021a.
- Generalist decay effect sizes added to the Results.
- Old in-document SI removed.
- Length: text about 7,100 words plus 830 words of captions, about 23.8 PU (limit 25).

## 3 Oct 2026 — Wood condition against every measurement (Bradford #39); authorship

- **Test:** wood condition is now tested against every flux measurement, not tree means. Comment 39 asked for this.
  - Model: flux ~ metric (per SD within group) + (1 | tree) + (1 | date), on untransformed flux.
  - Kenward–Roger df. α = 0.10 for these tree-level tests only (10 trees per group), stated in the Methods.
  - Effects are given as nmol m⁻² s⁻¹ per SD of ERT CV.
- **Results (p-values):**

  | Group | p |
  |---|---|
  | *N. sylvatica* | 0.031 |
  | *Q. rubra* | 0.025 |
  | Upland *A. rubrum* | 0.065 |
  | Wetland *A. rubrum* | 0.37 |
  | *T. canadensis* | n.s. |
  | Pooled, species-adjusted | wetland 0.004; upland 0.015 |

  - The conclusions are the same on the asinh scale.
  - The earlier "no metric significant in species-adjusted models" came from an asinh(1000 × flux) analysis. It has been removed.
- **Figures 5 and S5:** grey measurements behind the tree means. Linear axes are zoomed to the tree means, with subtle edge markers for measurements off the axis.
- **Table S2:** same model throughout, with SoT structural loss, ERT CV and the ERT index; ERT mean dropped.
- **Main text:** ERT CV only. The other metrics are referred to Table S2.
- **Authorship:** Mark A. Bradford (Yale School of the Environment; The Forest School) added as a co-author before Marra.
- **Acknowledgements:** now thank Peter Raymond only. Co-author comments in the draft were checked; none are from Raymond.
- **Funding:** NSF DEB-1945921/DEB-2231681 (corrected) and REU DBI-1950364/DBI-2348924 added.
- **Data citations:** AmeriFlux and NEON datasets cited by DOI, with the required acknowledgement statements.

## 3 Oct 2026 (later) — Submission-readiness review; upland framing

### Review fixes

- Figure S4 caption and text: "averaged across species".
- Figure S3 is now cited for the tree effects.
- Abstract and Discussion state the 58% as the model including species, with species alone at 31%.
- Table S6 upland reference species: *T. canadensis*.
- Text S2 sample size: n = 722.
- Window-screening permutation p now uses the +1 correction.
- Repeatability "scale" wording corrected.
- Decay-class denominator: 57 trees.
- Species-comparison methods added: mixed models on untransformed flux, emmeans, KR df, Tukey.
- Barker-Plotkin duplicate merged; Bates `\textbf` artifact removed.
- Corresponding author email: jgewirtz@stanford.edu.
- Variance partitioning recomputed on the model scale: it had used asinh(1000 × nmol).

### Upland framing

- Upland stems are small, persistent net sources, most likely from in-stem production partly offset by methanotrophy. This replaces "competing processes near equilibrium".
- Basis:

  | Analyzer | Detectable uptake | Detectable emission |
  |---|---|---|
  | LI-7810, 2025 (MDF 0.004) | 1% | 83% |
  | LGR, 2023–24 (MDF ~0.1) | 5.5% | 17% |
  | LGR, restricted to 2025 months, soil temperature and moisture | 4.5% | 20% |
  | LI-7810, matched comparison | 1.8% | 86% |

  - LGR uptake shows no soil-moisture dependence and sits just above the MDF.
  - Numbers: `outputs/tables/manuscript/upland_detection_by_analyzer.csv`.
- Detected uptake is treated as mostly analytical noise. Upland stems "rarely, if ever, became net sinks at this height".

### Wording and Key Point 3

- "Smaller in magnitude than the MDF" replaces "below detection" throughout.
- Methods and Text S1 give the rationale for retaining and flagging all fluxes.
- Key Point 3 now reads: "Wetwood/decay associated with higher emissions where stem production dominates and lower emissions where transport dominates."
