#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --partition=mulan
#SBATCH --mem=100G
#SBATCH --job-name=PRECAST

### ====
## This script runs PRECAST on simulated data
### ===

source ~/.bashrc
conda activate ~/envs/spatialATAC
# cd /net/mulan/home/wenjinma/projects/spatialATAC-explore/explore_PRECAST
cd /net/mulan/home/wenjinma/projects/SpaDOT_analysis/benchmark_study/explore_PRECAST

# declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
declare -a arr=(1993 2024 1105)
declare -a prop_list=(0.9 0.8 0.7 0.6 0.5)
declare -a DEG_list=(100 200 500)
for seed in "${arr[@]}"; do
    for prop in "${prop_list[@]}"; do
        for DEG in "${DEG_list[@]}"; do
            # Rscript run_PRECAST_simulation_prop.R -s $seed --prop $prop --DEG $DEG
            Rscript run_PRECAST_simulation_prop.R -s $seed --prop $prop --DEG $DEG --mixedNPCs
        done
    done
done
