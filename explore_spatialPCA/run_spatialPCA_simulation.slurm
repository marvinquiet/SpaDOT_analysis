#!/bin/bash
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --time=72:00:00
##SBATCH --partition=main
##SBATCH --mem=100G
##SBATCH --job-name=spatialPCA

### ====
## This script runs spatialPCA on simulated data
### ===

source ~/.bashrc
conda activate /media/sda/wenjinma/envs/spatialATAC
cd /media/sda/wenjinma/projects/spatialATAC-explore/explore_spatialPCA

declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
for seed in "${arr[@]}"
do
    echo "=== Run spatialPCA on simulation data for $seed ==="
    Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 20 -y 50 -s $seed
    #Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 10 -y 50 -s $seed
    Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 20 -y 30 -s $seed
    #Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 10 -y 30 -s $seed
    Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 30 -y 50 -s $seed
    Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 30 -y 30 -s $seed
done
