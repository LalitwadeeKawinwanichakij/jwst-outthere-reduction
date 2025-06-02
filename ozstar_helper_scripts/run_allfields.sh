#!/bin/bash
# Usage: ./run_allfields.sh FIELDS/fields_leo_to_scl.txt

if [ -z "$1" ]; then
  echo "Usage: $0 <field_list_file>"
  exit 1
fi

FIELD_LIST="$1"

# Create log folder if needed
mkdir -p logs
SUMMARY_LOG="logs2/submission_summary.log"
echo "Submission summary - $(date)" > "$SUMMARY_LOG"

# Read field list into array
mapfile -t FIELDS < "$FIELD_LIST"

# Loop over each field and submit job
for FIELD in "${FIELDS[@]}"; do
  echo "[🚀] Launching pipeline for $FIELD..."

  LOGFILE="logs2/${FIELD}_pipeline.log"

  sbatch --job-name="${FIELD}_chained" \
         --export=FIELD="$FIELD" \
         slurm_submit_chained_jobs_external.sh > "$LOGFILE" 2>&1

  if grep -q "Submitted batch job" "$LOGFILE"; then
    JOB_ID=$(grep "Submitted batch job" "$LOGFILE" | awk '{print $4}')
    echo "[✅] $FIELD submitted successfully with Job ID $JOB_ID" | tee -a "$SUMMARY_LOG"
  else
    echo "[❌] $FIELD submission failed. Check $LOGFILE for details." | tee -a "$SUMMARY_LOG"
  fi

  sleep 0.5  # brief delay to reduce SLURM strain
done

