#!/bin/bash
#PBS -q gpuvolta
#PBS -P e14
#PBS -l ncpus=12
#PBS -l ngpus=1
#PBS -l mem=96GB
#PBS -l jobfs=1GB
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/e14+scratch/e14
#PBS -l wd
#PBS -M z5131023@unsw.edu.au

cd /g/data/e14/jb2381/StaircaseExperiments/single_interface_rundown

# Julia
export JULIA_DEPOT_PATH="/g/data/e14/jb2381/.julia"
export JULIA_NUM_THREADS=auto

# Run the experiment
julia +lts --project dns_nl_rundown_equaldiffs.jl > $PBS_JOBID.log