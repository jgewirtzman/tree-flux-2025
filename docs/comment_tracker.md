# Co-author comments on the working draft

Working draft: `DRAFT_ Contrasting controls on tree methane emissions in upland and wetland forests.docx`
(copied from Downloads on 30 Sep 2026; its text matches bioRxiv v1, 10.64898/2026.01.29.702553).
41 comments (Jan–Apr 2026) from Bradford, Marra, Jurado, Thompson, Lutz, Matthes and Gewirtzman.

Status key: **plan** = resolved by the existing workplan; **code** = needs a new analysis or figure change;
**text** = writing only; **done-local** = already coded but not yet committed.

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
