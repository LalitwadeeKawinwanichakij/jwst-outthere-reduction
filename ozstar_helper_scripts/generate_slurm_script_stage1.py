import argparse
def generate_stage1_slurm(field, cpus_per_task, mem_per_cpu, time):
    script = f"""#!/bin/bash
#SBATCH -p skylake                   # Partition
#SBATCH --job-name={field}_stage1   # Job name
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --cpus-per-task={cpus_per_task}           # Number of CPUs
#SBATCH --mem-per-cpu={mem_per_cpu}G             # Memory per CPU
#SBATCH --time={time}               # Max runtime
#SBATCH --output=FIELDS/{field}/logs/stage1_%j.out  # Standard output log
#SBATCH --error=FIELDS/{field}/logs/stage1_%j.err   # Standard error log

echo "Starting Stage 1 processing for field: {field}"

# Create necessary directories
mkdir -p FIELDS/{field}/RATE
mkdir -p FIELDS/{field}/logs

# Activate Pixi environment
eval "$(pixi activate -e jwst)"

# Run Stage 1 pipeline with GNU Parallel
echo "Running Stage 1 for {field} with {cpus_per_task} CPUs"
pixi run --no-lockfile-update --environment jwst parallel \
  --link -j {cpus_per_task} \
  ./workflow/scripts/runStage1.py --scratch ::: \
  FIELDS/{field}/UNCAL/* ::: \
  $(printf 'FIELDS/{field}/RATE/%s\\n' $(ls FIELDS/{field}/UNCAL/ | sed 's/uncal/rate/')) \
  > FIELDS/{field}/logs/stage1.log 2>&1

echo "Stage 1 processing completed for {field}."
"""
    return script
def main():
    # Example usage:
    #field = "sex-02"
    #cpus_per_task = 8
    #mem_per_cpu = 4  # in GB
    #time = "48:00:00"
    parser = argparse.ArgumentParser(description="Generate SLURM script for data Stage 1.")

    parser.add_argument("--field", type=str, required=True, help="Field name (e.g., sex-02)")
    parser.add_argument("--cpus", type=int, default=4, help="Number of CPUs per task (default: 4)")
    parser.add_argument("--mem", type=int, default=32, help="Memory per CPU in GB (default: 32GB)")
    parser.add_argument("--time", type=str, default="48:00:00", help="Time limit (default: 48:00:00)")
    # Parse arguments
    args = parser.parse_args()
    slurm_script = generate_stage1_slurm(args.field, args.cpus, args.mem, 
            args.time)

    # Save the script to a file
    script_filename = f"{args.field}_stage1.slurm"
    with open(script_filename, "w") as file:
        file.write(slurm_script)

    print(f"SLURM script for {args.field} created: {args.field}_stage1.slurm")
if __name__ == "__main__":
    main()
