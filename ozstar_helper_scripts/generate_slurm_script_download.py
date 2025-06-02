import os
import argparse
from astropy.table import Table
def generate_slurm_script(field, cpus_per_task, mem_per_cpu, time):
    script = f"""#!/bin/bash
#SBATCH -p trevor                   # Partition
#SBATCH --job-name={field}_download   # Job name
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --cpus-per-task={cpus_per_task}           # Number of CPU cores
#SBATCH --mem-per-cpu={mem_per_cpu}G            # Memory per CPU
#SBATCH --time={time}             # Max runtime
#SBATCH --output=FIELDS/{field}/logs/download_%j.out  # Standard output log
#SBATCH --error=FIELDS/{field}/logs/download_%j.err   # Error log

# Create necessary directories
mkdir -p FIELDS/{field}/logs

# Run download.py script
echo "Starting download for {field}..."
pixi run --environment jwst python ./workflow/scripts/download.py {field} \
  --ncpu {cpus_per_task} > FIELDS/{field}/logs/download.log 2>&1

# Run CRDS sync to fetch reference files
echo "Starting CRDS sync for {field}..."
pixi run --no-lockfile-update --environment jwst crds sync \
  --contexts "$CRDS_CONTEXT" \
  --fetch-references \
  --dataset-files FIELDS/{field}/UNCAL/* \
  >> FIELDS/{field}/logs/download.log 2>&1

echo "Job completed!"
"""
    return script
def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Generate SLURM script for data download.")
    parser.add_argument("--field", type=str, required=True, help="Field name (e.g., sex-02)")
    parser.add_argument("--cpus", type=int, default=8, help="Number of CPUs per task")
    parser.add_argument("--mem", type=int, default=2, help="Memory per CPU (in GB)")
    parser.add_argument("--time", type=str, default="48:00:00", help="Time limit (e.g., 48:00:00)")

    # Parse arguments
    args = parser.parse_args()

    # Example usage:
    #field = "sex-02"
    #cpus_per_task = 8
    #mem_per_cpu = 2  # in GB
    #time = "48:00:00"

    slurm_script = generate_slurm_script(args.field, args.cpus, 
            args.mem, args.time)

    # Save the script to a file
    script_filename = f"{args.field}_download.slurm"
    with open(script_filename, "w") as file:
        file.write(slurm_script)

    print(f"SLURM script for {args.field} created: {args.field}_download.slurm")
if __name__ == "__main__":
    main()
