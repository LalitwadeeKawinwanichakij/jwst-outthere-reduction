import argparse
import os
from astropy.table import Table

def generate_slurm_script(field, cpus_per_task, mem_per_cpu, time, expected_files):
    script = f"""#!/bin/bash
#SBATCH -p trevor
#SBATCH --job-name={field}_download
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem-per-cpu={mem_per_cpu}G
#SBATCH --time={time}
#SBATCH --output=FIELDS/{field}/logs/download_%j.out
#SBATCH --error=FIELDS/{field}/logs/download_%j.err

mkdir -p FIELDS/{field}/logs

echo "Starting download for {field}..."
pixi run --environment jwst python ./workflow/scripts/download.py {field} \\
  --ncpu {cpus_per_task} > FIELDS/{field}/logs/download.log 2>&1

echo "Starting CRDS sync for {field}..."
pixi run --no-lockfile-update --environment jwst crds sync \\
  --contexts "$CRDS_CONTEXT" \\
  --fetch-references \\
  --dataset-files FIELDS/{field}/UNCAL/* \\
  >> FIELDS/{field}/logs/download.log 2>&1

echo "Job completed!"

# --- Check if all expected UNCAL files are present ---
echo "Verifying UNCAL files for {field}..."
uncal_dir="$(pwd)/FIELDS/{field}/UNCAL"
uncal_count=$(ls "$uncal_dir"/*uncal.fits 2>/dev/null | wc -l)
expected_count={expected_files}

echo "Found $uncal_count / $expected_count UNCAL files."

if [ "$uncal_count" -eq "$expected_count" ]; then
  echo "[INFO] File count matches. Generating and submitting stage1 job..."
  python generate_slurm_script_stage1.py --field {field} --cpu 4 --mem 32 --time 72:00:00
  if [ $? -eq 0 ]; then
    sbatch {field}_stage1.slurm
    echo "[INFO] stage1 job submitted."
  else
    echo "[ERROR] Failed to generate stage1 script."
  fi
else
  echo "[INFO] UNCAL file count mismatch — skipping stage1 submission."
fi
"""
    return script

def main():
    parser = argparse.ArgumentParser(description="Generate SLURM script for data download.")
    parser.add_argument("--field", type=str, required=True, help="Field name (e.g., sex-02)")
    parser.add_argument("--cpus", type=int, default=8, help="Number of CPUs per task")
    parser.add_argument("--mem", type=int, default=2, help="Memory per CPU (in GB)")
    parser.add_argument("--time", type=str, default="48:00:00", help="Time limit (e.g., 48:00:00)")

    args = parser.parse_args()

    # Read expected number of products from the FITS table
    main_dir = os.getcwd()
    fields_dir = os.path.join(main_dir, "FIELDS")
    fname = args.field
    fits_path = os.path.join(fields_dir, "fields.fits")

    try:
        prods = Table.read(fits_path, fname)
        expected_count = len(prods["productFilename"])
    except Exception as e:
        print(f"Error reading product count from FITS file: {e}")
        return

    slurm_script = generate_slurm_script(
        field=args.field,
        cpus_per_task=args.cpus,
        mem_per_cpu=args.mem,
        time=args.time,
        expected_files=expected_count
    )

    script_filename = f"{args.field}_download.slurm"
    with open(script_filename, "w") as file:
        file.write(slurm_script)

    print(f"SLURM script for {args.field} created: {script_filename}")

if __name__ == "__main__":
    main()

