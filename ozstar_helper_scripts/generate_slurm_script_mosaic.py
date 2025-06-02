import argparse
def generate_mosaic_slurm(field, partition='skylake', cpus_per_task=8, mem_per_cpu=4, time='48:00:00'):
    script = f"""#!/bin/bash
#SBATCH -p {partition}                   # Partition
#SBATCH --job-name={field}_mosaic       # Job name
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task={cpus_per_task}           # Number of CPUs
#SBATCH --mem-per-cpu={mem_per_cpu}G             # Memory per CPU
#SBATCH --time={time}                   # Max runtime
#SBATCH --output=FIELDS/{field}/logs/mosaic_%j.out  # Standard output log
#SBATCH --error=FIELDS/{field}/logs/mosaic_%j.err   # Standard error log

echo "Starting mosaic creation for field: {field}"

# Create necessary directories
mkdir -p FIELDS/{field}/logs

# Activate Pixi environment
eval "$(pixi activate -e grizli)"

# Run mosaic script with Pixi
echo "Running mosaic.py for {field} with {cpus_per_task} CPUs"
pixi run --no-lockfile-update --environment grizli \
  ./workflow/scripts/mosaic.py {field} \
  > FIELDS/{field}/logs/mos.log 2>&1

echo "Mosaic creation completed for {field}."
"""
    return script

# Example usage:
#field = "sex-01"
#partition = "skylake"
#cpus_per_task =4
#mem_per_cpu = 8  # in GB
#time = "72:00:00"
def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Generate SLURM script for data download.")
    # `--field` is required
    parser.add_argument("--field", type=str, required=True, help="Field name (e.g., sex-02)")
    # Other arguments have default values
    parser.add_argument("--cpus", type=int, default=4, help="Number of CPUs per task (default: 4)")
    parser.add_argument("--mem", type=int, default=8, help="Memory per CPU in GB (default: 8GB)")
    parser.add_argument("--time", type=str, default="48:00:00", help="Time limit (default: 48:00:00)")
    args = parser.parse_args()
    slurm_script = generate_mosaic_slurm(args.field,cpus_per_task=args.cpus, mem_per_cpu=args.mem,time=args.time)
    # Save the script to a file
    script_filename = f"{args.field}_mosaic.slurm"
    with open(script_filename, "w") as file:
        file.write(slurm_script)

    print(f"SLURM script for {args.field} created: {args.field}_mosaic.slurm")
if __name__ == "__main__":
    main()
