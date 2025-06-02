import argparse

def generate_slurm_script(fieldname, cpu, mem, time):
    script_content = f"""#!/bin/bash
#SBATCH -p skylake                   # Partition
#SBATCH --job-name={fieldname}_zfitew     # Job name
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --cpus-per-task={cpu}        # Number of CPUs
#SBATCH --mem-per-cpu={mem}G         # Memory per CPU
#SBATCH --time={time}          # Max runtime (hours)
#SBATCH --tmp=20GB
#SBATCH --output=FIELDS/{fieldname}/logs/zfitEW_%j.out  # Standard output log
#SBATCH --error=FIELDS/{fieldname}/logs/zfitEW_%j.err   # Standard error log

echo "Starting redshift fitting for field: {fieldname}"

# Create necessary directories
mkdir -p FIELDS/{fieldname}/logs
# Create a directory on the local disk
#mkdir -p $JOBFS/data/FIELDS/{fieldname}
# Activate Pixi environment
#eval "$(pixi activate -e grizli)"
export MPLCONFIGDIR=$JOBFS/data/matplotlib_cache
mkdir -p $MPLCONFIGDIR
# Run redshift fitting script with Pixi
echo "Running redshiftFit.py for {fieldname}"
if [ ! -d "/fred/oz408/NISPureParallel/FIELDS/{fieldname}/RedshiftFitting" ]; then
    ./create_symlink_extractions.sh {fieldname}
fi

pixi run --no-lockfile-update --environment grizli ./resources/scripts/Postprocessing_logLikeEW_4OzStar_parallel.py --field {fieldname} --output_dir /fred/oz408/NISPureParallel/FIELDS/{fieldname}/RedshiftFitting > FIELDS/{fieldname}/logs/zfitEW.log 2>&1

echo "Redshift fitting with EW prior completed for {fieldname}."
# tar up the files and copy it back to /fred
#cd $JOBFS/data/lkawin/FIELDS/ 
#tar -czf {fieldname}_RedshiftFitting.tar.gz {fieldname}/RedshiftFitting
#mv $JOBFS/data/lkawin/FIELDS/{fieldname}_RedshiftFitting.tar.gz /fred/oz041/lkawin/NISPureParallel/FIELDS/{fieldname}
"""

    # Write the script to a file
    output_filename = f"{fieldname}_zfitew.slurm"
    with open(output_filename, "w") as f:
        f.write(script_content)

    print(f"SLURM job submission script generated: {output_filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a SLURM job submission script for redshift fitting.")
    parser.add_argument("--fieldname", required=True, help="Field name (e.g., sex-00)")
    parser.add_argument("--cpu", type=int, default=8, help="Number of CPUs per task (default: 1)")
    parser.add_argument("--mem", type=int, default=2, help="Memory per CPU in GB (default: 16)")
    parser.add_argument("--time", type=str, default="72:00:00", help="Maximum runtime in hours (default: 72)")
    
    args = parser.parse_args()

    # Generate the SLURM script
    generate_slurm_script(args.fieldname, args.cpu, args.mem, args.time)

