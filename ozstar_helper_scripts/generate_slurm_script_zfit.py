import argparse
def generate_zfit_slurm(field, partition='skylake', cpus_per_task=8, mem_per_cpu=4, time='48:00:00'):
    script = f"""#!/bin/bash
#SBATCH -p {partition}                   # Partition
#SBATCH --job-name={field}_zfit         # Job name
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task={cpus_per_task}           # Number of CPUs
#SBATCH --mem-per-cpu={mem_per_cpu}G             # Memory per CPU
#SBATCH --time={time}                   # Max runtime
#SBATCH --output=FIELDS/{field}/logs/zfit_%j.out  # Standard output log
#SBATCH --error=FIELDS/{field}/logs/zfit_%j.err   # Standard error log

echo "Starting redshift fitting for field: {field}"

# Create necessary directories
mkdir -p FIELDS/{field}/logs

# Activate Pixi environment
eval "$(pixi activate -e grizli)"

# Run redshift fitting script with Pixi
echo "Running redshiftFit.py for {field} with {cpus_per_task} CPUs"
pixi run --no-lockfile-update --environment grizli \
  ./workflow/scripts/redshiftFit.py {field} \
  --ncpu {cpus_per_task} \
  > FIELDS/{field}/logs/zfit.log 2>&1

echo "Redshift fitting completed for {field}."
"""
    return script

# Example usage:
#field = "sex-01"
#partition = "skylake"
#cpus_per_task =16
#mem_per_cpu = 1 # in GB
#time = "72:00:00"
def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Generate SLURM script for data redshift fitting.")

    # `--field` is required
    parser.add_argument("--field", type=str, required=True, help="Field name (e.g., sex-02)")
    # Other arguments have default values
    parser.add_argument("--cpus", type=int, default=16, help="Number of CPUs per task (default: 16)")
    parser.add_argument("--mem", type=int, default=4, help="Memory per CPU in GB (default: 4GB)")
    parser.add_argument("--time", type=str, default="72:00:00", help="Time limit (default: 72:00:00)")
    # Parse arguments
    args = parser.parse_args()
    slurm_script = generate_zfit_slurm(args.field, cpus_per_task=args.cpus, mem_per_cpu=args.mem, time=args.time)

    # Save the script to a file
    with open(f"{args.field}_zfit.slurm", "w") as file:
        file.write(slurm_script)

    print(f"SLURM script for {args.field} created: {args.field}_zfit.slurm")
if __name__ == "__main__":
    main()
