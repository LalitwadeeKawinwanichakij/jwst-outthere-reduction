import argparse
def generate_preprocess_slurm(field, partition='trevor', cpus_per_task=8, mem_per_cpu=4, time='48:00:00'):
    script = f"""#!/bin/bash
#SBATCH -p {partition}                   # Partition
#SBATCH --job-name={field}_preprocess   # Job name
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task={cpus_per_task}           # Number of CPUs
#SBATCH --mem-per-cpu={mem_per_cpu}G             # Memory per CPU
#SBATCH --time={time}                   # Max runtime
#SBATCH --output=FIELDS/{field}/logs/preprocess_%j.out  # Standard output log
#SBATCH --error=FIELDS/{field}/logs/preprocess_%j.err   # Standard error log

echo "Starting preprocessing for field: {field}"

# Create necessary directories
mkdir -p FIELDS/{field}/logs

# Run preprocess script with Pixi
echo "Running preprocess.py for {field} with {cpus_per_task} CPUs"
pixi run --no-lockfile-update --environment grizli \
  ./workflow/scripts/preprocess.py {field} \
  --ncpu {cpus_per_task} \
  > FIELDS/{field}/logs/proc.log 2>&1

echo "Preprocessing completed for {field}."
"""
    return script
def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Generate SLURM script for data preprocessing.")
    # `--field` is required
    parser.add_argument("--field", type=str, required=True, help="Field name (e.g., sex-02)")

    # Other arguments have default values
    parser.add_argument("--cpus", type=int, default=8, help="Number of CPUs per task (default: 8)")
    parser.add_argument("--mem", type=int, default=1, help="Memory per CPU in GB (default: 1GB)")
    parser.add_argument("--time", type=str, default="72:00:00", help="Time limit (default: 72:00:00)")
    # Parse arguments
    args = parser.parse_args()
    # Example usage:
    #field = "sex-01"
    #partition = "trevor"
    #cpus_per_task = 8
    #mem_per_cpu = 1  # in GB
    #time = "72:00:00"

    slurm_script = generate_preprocess_slurm(args.field, cpus_per_task=args.cpus, mem_per_cpu=args.mem, time=args.time)

    # Save the script to a file
    script_filename = f"{args.field}_preprocess.slurm"
    with open(script_filename, "w") as file:
        file.write(slurm_script)

    print(f"SLURM script for {args.field} created: {args.field}_preprocess.slurm")
if __name__ == "__main__":
    main()
