#PBS -q gpuvolta
#PBS -P e14
#PBS -l ncpus=12
#PBS -l ngpus=1
#PBS -l mem=96GB
#PBS -l jobfs=1GB
#PBS -l walltime=24:00:00
#PBS -l storage=gdata/e14+scratch/e14
#PBS -l wd
#PBS -M z5131023@unsw.edu.au

cd /g/data/e14/jb2381/StaircaseExperiments/

# Julia
export JULIA_DEPOT_PATH="/g/data/e14/jb2381/.julia"
export JULIA_NUM_THREADS=auto
module load julia

cd stability_experiments
# Run the experiment
julia --project stablestaircase_test_lineareos.jl > $PBS_JOBID.log