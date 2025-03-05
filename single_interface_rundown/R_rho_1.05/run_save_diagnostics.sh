#!/bin/bash
#PBS -q normalbw
#PBS -P e14
#PBS -l ncpus=1
#PBS -l mem=256GB
#PBS -l walltime=15:00:00
#PBS -l storage=gdata/e14+scratch/e14
#PBS -l wd
#PBS -M z5131023@unsw.edu.au

cd /g/data/e14/jb2381/StaircaseExperiments/single_interface_rundown/R_rho_1.05

# Julia
export JULIA_DEPOT_PATH="/g/data/e14/jb2381/.julia"
export JULIA_NUM_THREADS=auto
module load julia

# Run the experiment
julia --project save_diagnostics.jl > $PBS_JOBID.log