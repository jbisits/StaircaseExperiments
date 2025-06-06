#!/bin/bash
#PBS -q normalbw
#PBS -P e14
#PBS -l ncpus=1
#PBS -l mem=256GB
#PBS -l walltime=4:00:00
#PBS -l storage=gdata/e14+scratch/e14
#PBS -l wd
#PBS -M z5131023@unsw.edu.au

cd /g/data/e14/jb2381/StaircaseExperiments/single_interface_rundown/

# Julia
export JULIA_DEPOT_PATH="/g/data/e14/jb2381/.julia"
export JULIA_NUM_THREADS=auto

# Run the experiment
julia --project re_run_diagnostics.jl > $PBS_JOBID.log