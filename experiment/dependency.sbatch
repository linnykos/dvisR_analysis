#! /bin/bash
#SBATCH --job-name=dependency
#SBATCH --output=dep_sbatch.out
#SBATCH --error=dep_error.err
#SBATCH --account=pi-jingshuw
#SBATCH --partition=caslake
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=2000

module load gsl
module load R
Rscript dependency.r
