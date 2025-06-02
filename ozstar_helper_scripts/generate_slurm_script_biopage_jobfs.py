import argparse

SLURM_TEMPLATE = """#!/bin/bash
#SBATCH -p skylake
#SBATCH --job-name={field}_biopage
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem-per-cpu={mem_per_cpu}G
#SBATCH --time={time}
#SBATCH --tmp=20G
#SBATCH --output={root}/summaryplots/{field}_biopage_%j.out
#SBATCH --error={root}/summaryplots/{field}_biopage_%j.err

export MPLCONFIGDIR=$JOBFS/data/lkawin/matplotlib_cache
mkdir -p $MPLCONFIGDIR
mkdir -p $JOBFS/data/lkawin/{field}/Summary_Plots
mkdir -p {root}/FIELDS/{field}

unshare --user --pid --map-root-user --fork --mount-proc <<EOF
    cd {root}
    squashfuse {field}.sqfs {root}/FIELDS/{field}
    ls {root}/FIELDS/{field}

    pixi run --no-lockfile-update --environment grizli python PlotGrizli_summaryplot_4OzStar.py \\
        --field {field} --output_dir $JOBFS/data/lkawin/{field}/Summary_Plots --cpus {cpus}
EOF

cd $JOBFS/data/lkawin/{field}
tar -czf {field}_Summary_Plots.tar.gz Summary_Plots
mv {field}_Summary_Plots.tar.gz {root}/summaryplots/
"""

def generate_slurm_script(field, cpus=8, mem_per_cpu=16, time="48:00:00", root="/fred/oz408/NISPureParallel"):
    script_content = SLURM_TEMPLATE.format(field=field, cpus=cpus, mem_per_cpu=mem_per_cpu, time=time, root=root)
    script_filename = f"{field}_biopage.slurm"
    with open(script_filename, "w") as f:
        f.write(script_content)
    print(f"SLURM job script generated: {script_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate SLURM job for BioPage summary.")
    parser.add_argument("--field", type=str, required=True, help="Field name (e.g. aqr-00)")
    parser.add_argument("--cpu", type=int, default=8, help="CPUs per task (default: 8)")
    parser.add_argument("--mem", type=int, default=16, help="Memory per CPU in GB (default: 16)")
    parser.add_argument("--time", type=str, default="48:00:00", help="Time limit (e.g. 48:00:00)")
    parser.add_argument("--root", type=str, default="/fred/oz408/NISPureParallel", help="Root path to FIELDS directory")

    args = parser.parse_args()
    generate_slurm_script(args.field, args.cpu, args.mem, args.time, args.root)

