#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00
#SBATCH --partition=mulan
#SBATCH --mem=100G
#SBATCH --job-name=spatialPCA

### ====
## This script runs spatialPCA on simulated data
### ===

source ~/.bashrc
# conda activate /media/sda/wenjinma/envs/spatialATAC
conda activate ~/envs/spatialATAC
# cd /media/sda/wenjinma/projects/spatialATAC-explore/explore_spatialPCA
cd /net/mulan/home/wenjinma/projects/SpaDOT_analysis/benchmark_study/explore_spatialPCA

# declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
# for seed in "${arr[@]}"
# do
#     echo "=== Run spatialPCA on simulation data for $seed ==="
#     # Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 20 -y 50 -s $seed
#     # #Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 10 -y 50 -s $seed
#     # Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 20 -y 30 -s $seed
#     # #Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 10 -y 30 -s $seed
#     # Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 30 -y 50 -s $seed
#     # Rscript run_spatialPCA_large-scale_simulation_integrated.R -x 30 -y 30 -s $seed
# done

declare -a arr=(1993 2024 1105)
declare -a prop_list=(0.9 0.8 0.7 0.6 0.5)
declare -a DEG_list=(100 200 500)
for seed in "${arr[@]}"; do
    for prop in "${prop_list[@]}"; do
        for DEG in "${DEG_list[@]}"; do
            # Rscript run_spatialPCA_simulation_prop_integrated.R --seed $seed --prop $prop --DEG $DEG
            Rscript run_spatialPCA_simulation_prop_integrated.R --seed $seed --prop $prop --DEG $DEG --mixedNPCs
        done
    done
done
