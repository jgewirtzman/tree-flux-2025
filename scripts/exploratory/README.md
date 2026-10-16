# Exploratory scripts (not used in the paper)

These scripts are not run by `scripts/run_pipeline.sh`, and nothing in the paper, the
Supporting Information or `07_manuscript_numbers.R` reads their output. They are kept for
reference. Run them from the project root after the pipeline, which makes their inputs.

| Script | Was | What it does | Output |
|---|---|---|---|
| `decay_binned.R` | `4_analysis/10_decay_binned.R` | binned and threshold versions of the wood condition vs CH4 flux analysis (median split, tertiles, decay classes, best split and hinge with permutation p) | `outputs/tables/decay_binned.csv`, `decay_effect_sizes.csv` |
| `compare_upland_models.R` | `5_models/05_compare_upland_models.R` | prints upland Model A vs Model B fit and species slopes (superseded by section 6 of `5_models/06_model_checks.R`, which compares them on shared rows) | console only |
