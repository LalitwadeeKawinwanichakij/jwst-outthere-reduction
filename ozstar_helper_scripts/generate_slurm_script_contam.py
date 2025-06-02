import argparse
def generate_contam_slurm(field, partition='skylake', cpus_per_task=8, mem_per_cpu=4, time='48:00:00'):
    script = f"""#!/bin/bash
#SBATCH -p {partition}                   # Partition
#SBATCH --job-name={field}_contam       # Job name
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task={cpus_per_task}           # Number of CPUs
#SBATCH --mem-per-cpu={mem_per_cpu}G             # Memory per CPU
#SBATCH --time={time}                   # Max runtime
#SBATCH --output=FIELDS/{field}/logs/contam_%j.out  # Standard output log
#SBATCH --error=FIELDS/{field}/logs/contam_%j.err   # Standard error log

echo "Starting contamination detection for field: {field}"

# Create necessary directories
mkdir -p FIELDS/{field}/logs

# Activate Pixi environment
eval "$(pixi activate -e grizli)"

# Run contamination detection script with Pixi
echo "Running contamination.py for {field} with {cpus_per_task} CPUs"
pixi run --no-lockfile-update --environment grizli /fred/oz041/lkawin/NISPureParallel/workflow/scripts/contamination.py {field} --ncpu {cpus_per_task} > FIELDS/{field}/logs/contam.log 2>&1

echo "Contamination detection completed for {field}."
"""
    return script
def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Generate SLURM script for contamination modeling.")
    parser.add_argument("--field", type=str, required=True, help="Field name (e.g., sex-02)")
    parser.add_argument("--cpus", type=int, default=4, help="Number of CPUs per task")
    parser.add_argument("--mem", type=int, default=16, help="Memory per CPU (in GB)")
    parser.add_argument("--time", type=str, default="72:00:00", help="Time limit (e.g., 72:00:00)")

    # Parse arguments
    args = parser.parse_args()
    # Example usage:
    #field = "sex-01"
    #partition = "skylake"
    #cpus_per_task = 8
    #mem_per_cpu = 4  # in GB
    #time = "72:00:00"

    slurm_script = generate_contam_slurm(args.field, 
            cpus_per_task=args.cpus, mem_per_cpu=args.mem, time=args.time)

    # Save the script to a file
    script_filename = f"{args.field}_contam.slurm"
    with open(script_filename, "w") as file:
        file.write(slurm_script)

    print(f"SLURM script for {args.field} created: {args.field}_contam.slurm")
if __name__ == "__main__":
    main()
