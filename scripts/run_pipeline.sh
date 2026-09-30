#!/bin/bash
# Rerun the flux-dependent pipeline end to end (environmental imports 01-08
# do not depend on the flux data and are not rerun). Logs: outputs/logs/.
#   bash scripts/run_pipeline.sh            # full run
#   SKIP_GOFLUX=1 bash scripts/run_pipeline.sh   # reuse the existing goFlux table
cd "$(dirname "$0")/.."
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8   # non-ASCII plot labels (°C, ×) need a UTF-8 locale
LOG=outputs/logs; mkdir -p "$LOG"; : > "$LOG/STATUS.txt"
steps=()
[ "${SKIP_GOFLUX:-0}" = "1" ] || steps+=(scripts/01_import/09_goflux_reprocess.R)
steps+=(scripts/01_import/10_quality_flags.R scripts/01_import/11_tomography_classes.R
  scripts/02_analysis/00_data_summary.R scripts/02_analysis/01_timeseries.R scripts/02_analysis/02_flux_summaries.R
  scripts/02_analysis/03_repeatability.R scripts/02_analysis/04_rolling_corrs.R scripts/02_analysis/05_combined_driver_timeseries.R
  scripts/02_analysis/06_filtering_snr.R scripts/02_analysis/07_filter_sensitivity_ridges.R scripts/02_analysis/08_tomography.R
  scripts/03_modeling/01_bgs_model.R scripts/03_modeling/02_ems_model_A.R scripts/03_modeling/03_ems_model_B.R
  scripts/03_modeling/04_interaction_plots.R scripts/03_modeling/05_compare_models.R scripts/03_modeling/06_manuscript_numbers.R)
[ -f scripts/04_publish/00_goflux_attributes.R ] && steps+=(scripts/04_publish/00_goflux_attributes.R)
for s in "${steps[@]}"; do
  n=$(basename "$s" .R); t0=$(date +%s)
  Rscript "$s" > "$LOG/$n.log" 2>&1; rc=$?
  echo "$n rc=$rc $(( $(date +%s) - t0 ))s" | tee -a "$LOG/STATUS.txt"
  [ $rc -ne 0 ] && { echo "FAILED: $s (see $LOG/$n.log)"; exit $rc; }
done
echo ALLDONE | tee -a "$LOG/STATUS.txt"
