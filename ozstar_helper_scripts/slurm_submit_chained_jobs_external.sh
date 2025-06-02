#!/bin/bash
#SBATCH -p skylake                   # Partition
#SBATCH --job-name=chain
#SBATCH --output=logs2/%x_%j.out
#SBATCH --error=logs2/%x_%j.err
#SBATCH --time=72:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
echo "[INFO] Running chained jobs for field: $FIELD"
# Example: sbatch submit_chained_jobs_external.sh aqr-00

set -e

skip_downstream=false
# === CONFIGURATION ===
root="/fred/oz408/NISPureParallel"
DL_CPUS=8
DL_MEM=2
DL_TIME="48:00:00"

STAGE1_CPUS=4
STAGE1_MEM=64
STAGE1_TIME="72:00:00"

PRE_CPUS=8
PRE_MEM=1
PRE_TIME="72:00:00"

MOS_CPUS=4
MOS_MEM=8
MOS_TIME="48:00:00"

CONTAM_CPUS=4
CONTAM_MEM=16
CONTAM_TIME="72:00:00"

EXTRACT_CPUS=8
EXTRACT_MEM=8
EXTRACT_TIME="72:00:00"

ZFIT_CPUS=7
ZFIT_MEM=6
ZFIT_TIME="72:00:00"

# BioPage step
BIOPAGE_CPUS=8
BIOPAGE_MEM=2
BIOPAGE_TIME="48:00:00"

# Refit EW Prior step
REFIT_CPUS=8
REFIT_MEM=2
REFIT_TIME="72:00:00"

# === FUNCTIONS ===
is_job_pending_or_running() {
  local job_name="$1"
  squeue -u "$USER" --noheader --format="%.100j" | awk '{$1=$1; print}' | grep "$job_name" >/dev/null
}


check_job_completion() {
  local jobid=$1
  local stepname=$2
  echo "[INFO] Waiting for $stepname job (ID $jobid) to complete..."
  while true; do
    status=$(sacct -j "$jobid" --format=JobID,State --parsable2 --noheader |
             grep -E "^$jobid(\\||$)" | head -n 1 | cut -d'|' -f2 | tr -d '[:space:]')
    if [ -z "$status" ]; then
      echo "[INFO] $stepname: Job status not available yet. Waiting..."
    elif [[ "$status" =~ ^(COMPLETED|FAILED|CANCELLED|TIMEOUT)$ ]]; then
      echo "[INFO] $stepname job finished with status: $status"
      break
    else
      echo "[INFO] $stepname job still running: $status"
    fi
    sleep 30
  done
  if [ "$status" != "COMPLETED" ]; then
    echo "[ERROR] $stepname job did not complete successfully. Aborting."
    exit 1
  fi
}

# === START ===

if [ -z "$FIELD" ]; then
  echo "Usage: $0 FIELDNAME"
  exit 1
fi

echo "[🚀] Starting job submission chain for field: $FIELD"

uncal_dir="FIELDS/$FIELD/UNCAL"
rate_dir="FIELDS/$FIELD/RATE"
proc_log="FIELDS/$FIELD/logs/proc.log"
visits_yaml="FIELDS/$FIELD/Prep/${FIELD}_visits.yaml"
phot_fits="FIELDS/$FIELD/Prep/${FIELD}_phot.fits"
ircat_fits="FIELDS/$FIELD/Prep/${FIELD}-ir.cat.fits"
extractions_dir="FIELDS/$FIELD/Extractions"
extracted_fits="$extractions_dir/${FIELD}-extracted.fits"
redshiftfit_dir="FIELDS/$FIELD/RedshiftFitting"

fits_table="FIELDS/fields.fits"

#expected_count=$(pixi run python -c "from astropy.table import Table; prods = Table.read('$fits_table', '$FIELD'); print(len(prods['productFilename']))")
expected_count=$(pixi run python -c 'from astropy.table import Table; prods = Table.read("'"$fits_table"'", hdu="'"$FIELD"'"); print(len(prods["productFilename"]))')
actual_uncal=$(ls "$uncal_dir"/*uncal.fits 2>/dev/null | wc -l || echo 0)
actual_rate=$(ls "$rate_dir"/*rate.fits 2>/dev/null | wc -l || echo 0)



# === Step 1: Download
download_jobname="${FIELD}_download"
if [ "$actual_uncal" -eq "$expected_count" ]; then
  echo "[INFO] UNCAL complete — skipping download."
elif is_job_pending_or_running "$download_jobname"; then
  echo "[INFO] Download job already in queue — skipping resubmission."
else
  pixi run python generate_slurm_script_download.py --field "$FIELD" --cpus $DL_CPUS --mem $DL_MEM --time $DL_TIME
  jid_download=$(sbatch ${FIELD}_download.slurm | awk '{print $4}')
  check_job_completion "$jid_download" "Download"
fi

# === Step 2: Stage1
stage1_jobname="${FIELD}_stage1"
if [ "$actual_rate" -eq "$expected_count" ]; then
  echo "[INFO] RATE complete — skipping stage1."
elif is_job_pending_or_running "$stage1_jobname"; then
  echo "[INFO] Stage 1 job already in queue — skipping resubmission."
else
  pixi run python generate_slurm_script_stage1.py --field "$FIELD" --cpu $STAGE1_CPUS --mem $STAGE1_MEM --time $STAGE1_TIME
  jid_stage1=$(sbatch ${FIELD}_stage1.slurm | awk '{print $4}')
  check_job_completion "$jid_stage1" "Stage1"
fi

# === Step 3: Preprocess
preproc_jobname="${FIELD}_preprocess"
preprocess_done=false
if [ -f "$proc_log" ]; then
  last_two=$(grep -v '^\s*$' "$proc_log" | tail -n 2 | tr -d '\r')
  if echo "$last_two" | grep -Fxq "kill='preprocess'" && echo "$last_two" | grep -q "Update exposure footprints"; then
    echo "[INFO] Preprocess complete — skipping."
    preprocess_done=true
  fi
fi

if [ "$preprocess_done" = false ]; then
  if is_job_pending_or_running "$preproc_jobname"; then
    echo "[INFO] Preprocess job already in queue — skipping submission."
  else
    echo "[INFO] Submitting preprocess job..."
    pixi run python generate_slurm_script_preprocess.py --field "$FIELD" --cpu $PRE_CPUS --mem $PRE_MEM --time $PRE_TIME
    jid_preprocess=$(sbatch ${FIELD}_preprocess.slurm | awk '{print $4}')
    check_job_completion "$jid_preprocess" "Preprocess"
  fi
fi

# === Step 4: Mosaic
mosaic_jobname="${FIELD}_mosaic"

if [ -f "$phot_fits" ] && [ -f "$ircat_fits" ]; then
  echo "[INFO] Mosaic outputs present — skipping mosaic."
elif [ -f "$visits_yaml" ]; then
  if is_job_pending_or_running "$mosaic_jobname"; then
    echo "[INFO] Mosaic job already in queue — skipping submission."
  else
    echo "[INFO] Submitting mosaic job..."
    pixi run python generate_slurm_script_mosaic.py --field "$FIELD" --cpu $MOS_CPUS --mem $MOS_MEM --time $MOS_TIME
    jid_mosaic=$(sbatch ${FIELD}_mosaic.slurm | awk '{print $4}')
    check_job_completion "$jid_mosaic" "Mosaic"
  fi
fi

# === Step 5: Contam
contam_jobname="${FIELD}_contam"
contam_log="FIELDS/$FIELD/logs/contam.log"
n_clean=$(ls "$extractions_dir"/*grism_clean.fits 2>/dev/null | wc -l || echo 0)

if [ "$n_clean" -gt 0 ]; then
  echo "[INFO] Found $n_clean grism_clean.fits files — skipping contam."
elif [ -f "$contam_log" ] && tail -n 1 "$contam_log" | grep -q "No grism files, skipping"; then
  echo "[⚠️] Skipping contam — no grism exposures available (see contam.log)."
  skip_downstream=true
elif is_job_pending_or_running "$contam_jobname"; then
  echo "[INFO] Contam job already in queue — waiting for completion..."
  jid_contam=$(get_job_id_by_name "$contam_jobname")
  check_job_completion "$jid_contam" "Contam"
elif [ -f "$phot_fits" ] && [ -f "$ircat_fits" ]; then
  echo "[INFO] Submitting contam job..."
  pixi run python generate_slurm_script_contam.py --field "$FIELD" --cpu $CONTAM_CPUS --mem $CONTAM_MEM --time $CONTAM_TIME
  jid_contam=$(sbatch ${FIELD}_contam.slurm | awk '{print $4}')
  check_job_completion "$jid_contam" "Contam"
else
  echo "[INFO] Skipping contam — missing required Prep files."
fi


if [ "$skip_downstream" = true ]; then
  echo "[⚠️] Skipping extract, zfit, biopage, and refit — no grism data available."
else
    # === Extract Completion Check (Step 6 Pre-check) ===
    extract_skip=false
    extracted_exists=false
    beam_count=0
    naxis2=0
    percent=0

    echo "[INFO] Checking extract step completion..."

    if [ -f "$extracted_fits" ]; then
	echo "[DEBUG] Found extracted file: $extracted_fits"
	extracted_exists=true
	beam_count=$(ls "$extractions_dir"/*beams.fits 2>/dev/null | wc -l || echo 0)
	#naxis2=$(pixi run python -c "from astropy.io import fits; h = fits.getheader('$extracted_fits', 1); print(h.get('NAXIS2', 0))")
	naxis2=$(pixi run python -c 'from astropy.io import fits; h = fits.getheader("'$extracted_fits'", 1); print(h.get("NAXIS2", 0))')

	echo "[DEBUG] beam_count=$beam_count"
	echo "[DEBUG] naxis2=$naxis2"

	if [ "$naxis2" -gt 0 ]; then
	    percent=$(awk "BEGIN { print 100 * $beam_count / $naxis2 }")
	    echo "[DEBUG] percent=$percent"
	    if (( $(echo "$percent > 95" | bc -l) )); then
		extract_skip=true
		echo "[INFO] Extract step already completed — skipping extract."
	    fi
	else
	    echo "[WARNING] NAXIS2 is 0 — skipping extract check."
	fi
    else
	echo "[INFO] No extracted.fits found — will run extract."
    fi


    # === Step 6: Extract ===
    extract_jobname="${FIELD}_extract"

    if [ "$extract_skip" = "false" ]; then
	if is_job_pending_or_running "$extract_jobname"; then
	    echo "[INFO] Extract job already in queue — skipping submission."
	else
	    echo "[INFO] Submitting extract job..."
	    pixi run python generate_slurm_script_extract.py --field "$FIELD" --cpu $EXTRACT_CPUS --mem $EXTRACT_MEM --time $EXTRACT_TIME
	    jid_extract=$(sbatch ${FIELD}_extract.slurm | awk '{print $4}')
	    check_job_completion "$jid_extract" "Extract"
	fi
    else
	echo "[INFO] Skipping extract — already complete."
    fi
    # === Re-check percent after extract if it was just run ===
    if [ "$extract_skip" = "false" ]; then
	if [ -f "$extracted_fits" ]; then
	    echo "[DEBUG] Rechecking extract completeness after job..."
	    extracted_exists=true  
	    beam_count=$(ls "$extractions_dir"/*beams.fits 2>/dev/null | wc -l || echo 0)
	    #naxis2=$(pixi run python -c "from astropy.io import fits; h = fits.getheader('$extracted_fits', 1); print(h.get('NAXIS2', 0))")
	    naxis2=$(pixi run python -c 'from astropy.io import fits; h = fits.getheader("'$extracted_fits'", 1); print(h.get("NAXIS2", 0))')
	    echo "[DEBUG] (post-extract) beam_count=$beam_count, naxis2=$naxis2"

	    if [ "$naxis2" -gt 0 ]; then
		percent=$(awk "BEGIN { print 100 * $beam_count / $naxis2 }")
		echo "[DEBUG] (post-extract) percent=$percent"
	    fi
	fi
    fi


    # === Step 7: ZFit ===
    fitresults_file="$extractions_dir/${FIELD}_fitresults.fits"
    zfit_jobname="${FIELD}_zfit"

    if [ -f "$fitresults_file" ]; then
	echo "[INFO] Fitresults already exist — skipping zfit."
    elif is_job_pending_or_running "$zfit_jobname"; then
	echo "[INFO] ZFit job already in queue — skipping resubmission."
    elif [ "$extracted_exists" = true ] && (( $(echo "$percent > 95" | bc -l) )); then
	echo "[INFO] Submitting zfit job..."
	pixi run python generate_slurm_script_zfit.py --field "$FIELD" --cpu $ZFIT_CPUS --mem $ZFIT_MEM --time $ZFIT_TIME
	jid_zfit=$(sbatch ${FIELD}_zfit.slurm | awk '{print $4}')
	check_job_completion "$jid_zfit" "ZFit"
    else
	echo "[INFO] Skipping zfit — either already completed or extract output incomplete."
    fi


    # === Step 8: BioPage Summary ===
    biopage_jobname="${FIELD}_biopage"
    n_full=$(ls "$extractions_dir"/*full.fits 2>/dev/null | wc -l || echo 0)
    fitresults_file="$extractions_dir/${FIELD}_fitresults.fits"
    summaryplot_dir="FIELDS/$FIELD/Summary_Plots"
    summary_tar="summaryplots/${FIELD}_Summary_Plots.tar.gz"
    n_pdf=$(ls "$summaryplot_dir"/*.pdf 2>/dev/null | wc -l || echo 0)

    if { [ "$n_full" -eq "$n_pdf" ] && [ "$n_pdf" -gt 0 ]; } || [ -f "$summary_tar" ]; then
	echo "[INFO] BioPage summary already complete — skipping."
    elif [ "$n_full" -ge 1 ] && [ -f "$fitresults_file" ]; then
	if is_job_pending_or_running "$biopage_jobname"; then
	    echo "[INFO] BioPage job already in queue — skipping submission."
	else
	    echo "[INFO] Submitting BioPage summary job..."
	    pixi run python generate_slurm_script_biopage.py --field "$FIELD" --cpu $BIOPAGE_CPUS --mem $BIOPAGE_MEM --time $BIOPAGE_TIME --root "$root" --useprior 0
	    jid_biopage=$(sbatch ${FIELD}_biopage.slurm | awk '{print $4}')
	    check_job_completion "$jid_biopage" "BioPage"
	fi
    else
	echo "[INFO] Skipping BioPage — full.fits or fitresults.fits missing."
    fi

    # === Step 9: Refit EW Prior ===
    refit_jobname="${FIELD}_zfitew"
    fitresults_ew_file="$redshiftfit_dir/${FIELD}_fitresults_useEWprior.fits"

    if [ -f "$fitresults_ew_file" ]; then
	echo "[INFO] RefitEWprior already completed — skipping."
    elif [ -f "$fitresults_file" ]; then
	if is_job_pending_or_running "$refit_jobname"; then
	    echo "[INFO] RefitEWprior job already in queue — waiting for completion..."
	    jid_zfitew=$(get_job_id_by_name "$zfitew_jobname")
	    check_job_completion "$jid_zfitew" "BioPage \(EW prior\)|RefitEWprior"
	else
	    echo "[INFO] Submitting refitEWprior job..."
	    pixi run python generate_slurm_script_refitEWprior.py --fieldname "$FIELD" --cpu $REFIT_CPUS --mem $REFIT_MEM --time $REFIT_TIME
	    jid_refit=$(sbatch ${FIELD}_zfitew.slurm | awk '{print $4}')
	    check_job_completion "$jid_refit" "RefitEWprior"
	fi
    else
	echo "[INFO] Skipping refitEWprior — fitresults file not found."
    fi


    # === Step 10: BioPage Summary (EW Prior) ===
    biopage_ew_jobname="${FIELD}_biopage"
    summaryplot_ew_dir="FIELDS/$FIELD/Summary_Plots_EWprior"
    n_ew_full=$(find "$redshiftfit_dir" -maxdepth 1 -type f ! -xtype l -name "*full.fits" | wc -l)
    n_ew_pdf=$(ls "$summaryplot_ew_dir"/*.pdf 2>/dev/null | wc -l || echo 0)

    if [ "$n_ew_full" -eq "$n_ew_pdf" ] && [ "$n_ew_pdf" -gt 0 ]; then
	echo "[INFO] BioPage summary EW prior already complete - skipping."
    elif [ "$n_ew_full" -ge 1 ] && [ -f "$fitresults_ew_file" ]; then
	if is_job_pending_or_running "$biopage_ew_jobname"; then
	    echo "[INFO] BioPage EW prior job already in queue - skipping submission."
	else
	    echo "[INFO] Submitting BioPage summary job EW prior..."
	    pixi run python generate_slurm_script_biopage.py --field "$FIELD" --cpu $BIOPAGE_CPUS --mem $BIOPAGE_MEM --time $BIOPAGE_TIME --root "$root" --useprior 1
	    jid_biopage_ew=$(sbatch ${FIELD}_biopage.slurm | awk '{print $4}')
	    check_job_completion "$jid_biopage_ew" "BioPage EW prior"
	fi
    else
	echo "[INFO] Skipping BioPage EW prior — full.fits or fitresults_useEWprior.fits missing."
    fi

fi 
contam_log="FIELDS/$FIELD/logs/contam.log"
extract_log="FIELDS/$FIELD/logs/extr.log"

if [ -f "$contam_log" ] && tail -n 1 "$contam_log" | grep -q "No grism files, skipping"; then
  echo "[⚠️] Pipeline incomplete for field: $FIELD — no grism exposures found."
elif [ -f "$extract_log" ] && tail -n 1 "$extract_log" | grep -q "Stellar field, skipping"; then
  echo "[⚠️] Pipeline incomplete for field: $FIELD — identified as stellar field skipped."
else
  echo "[✅] Pipeline complete for field: $FIELD"
fi
