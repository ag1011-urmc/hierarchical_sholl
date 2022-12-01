#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --job-name="Sims"
#SBATCH --mem=8GB
#SBATCH --output=scenario_%a.out
#SBATCH -c 4
#SBATCH -a 1-7
export OMP_NUM_THREADS=4

# This may vary per HPC system. At USC's hpc system
# we use: source /usr/usc/R/default/setup.sh
module load r/4.2.1

Rscript --vanilla MDND_Sim_Fit.R
