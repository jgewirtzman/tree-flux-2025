#!/bin/bash
# Run the analysis pipeline in order, from data/package/ and data/external/ to the
# figures, tables, manuscript numbers and Supporting Information.
#
#   bash scripts/run_pipeline.sh                 # stages 1-6
#   bash scripts/run_pipeline.sh 2 4             # only stages 2 and 4
#   SKIP_FIT=1 bash scripts/run_pipeline.sh      # reuse data/interim/flux_fits.csv (the fit takes ~20 min)
#
# Before the first run: get the data package into data/package/ and the public data into
# data/external/ (scripts/0_data/, see README.md). Stage 0 is not run here.
# Not run here: 2_flux/04_review_windows.R (interactive) and 7_publish/ (EDI upload).
# Logs: outputs/logs/<script>.log; one line per script in outputs/logs/STATUS.txt.
cd "$(dirname "$0")/.."
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8   # non-ASCII plot labels (°C, ×) need a UTF-8 locale
LOG=outputs/logs; mkdir -p "$LOG" data/interim data/final; : > "$LOG/STATUS.txt"

stage1=(1_environment/01_met_hydro.R 1_environment/02_soil_moisture.R 1_environment/03_neon_xha.R 1_environment/04_hourly_drivers.R)
stage2=(2_flux/01_closure_table.R 2_flux/02_fit_fluxes.R 2_flux/03_trace_qc.R 2_flux/05_deadband_sensitivity.R 2_flux/06_flux_dataset.R)
stage3=(3_trees/01_tomography_classes.R 3_trees/02_stand_context.R)
stage4=(4_analysis/01_data_summary.R 4_analysis/02_flux_timeseries.R 4_analysis/03_flux_by_species.R 4_analysis/04_repeatability.R
        4_analysis/05_window_screening.R 4_analysis/06_tomography_flux.R 4_analysis/07_decay_definitions.R 4_analysis/08_site_map.R
        4_analysis/09_driver_timeseries.R)
stage5=(5_models/01_wetland_model.R 5_models/02_upland_model_A.R 5_models/03_upland_model_B.R 5_models/04_interaction_plots.R
        5_models/05_compare_upland_models.R 5_models/06_model_checks.R 5_models/07_manuscript_numbers.R)
stage6=(6_manuscript/01_build_si.R)

stages=("$@"); [ ${#stages[@]} -eq 0 ] && stages=(1 2 3 4 5 6)
steps=()
for s in "${stages[@]}"; do
  eval "list=(\"\${stage$s[@]}\")"
  for f in "${list[@]}"; do
    if [ "${SKIP_FIT:-0}" = "1" ] && { [ "$f" = 2_flux/01_closure_table.R ] || [ "$f" = 2_flux/02_fit_fluxes.R ]; }; then continue; fi
    steps+=("scripts/$f")
  done
done

for s in "${steps[@]}"; do
  n=$(basename "$s" .R); t0=$(date +%s)
  Rscript "$s" > "$LOG/$n.log" 2>&1; rc=$?
  echo "$s rc=$rc $(( $(date +%s) - t0 ))s" | tee -a "$LOG/STATUS.txt"
  [ $rc -ne 0 ] && { echo "FAILED: $s (see $LOG/$n.log)"; exit $rc; }
done
echo ALLDONE | tee -a "$LOG/STATUS.txt"
