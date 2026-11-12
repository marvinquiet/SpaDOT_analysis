#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00
#SBATCH --partition=mulan
#SBATCH --mem=100G
#SBATCH --job-name=spaGCN

source ~/.bashrc
# conda activate /media/sda/wenjinma/envs/SpaGCN_py39
# cd /media/sda/wenjinma/projects/spatialATAC-explore/explore_spaGCN
conda activate /net/zootopia/disk1/wenjinma/envs/SpaGCN_py39
cd /net/mulan/home/wenjinma/projects/SpaDOT_analysis/benchmark_study/explore_spaGCN

# declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
declare -a arr=(1993 2024 1105)
declare -a prop_list=(0.9 0.8 0.7 0.6 0.5)
declare -a DEG_list=(100 200 500)
for seed in "${arr[@]}"; do
    for prop in "${prop_list[@]}"; do
        for DEG in "${DEG_list[@]}"; do
            # python run_spaGCN_simulation_prop.py --seed $seed --prop $prop --DEG $DEG
            python run_spaGCN_simulation_prop.py --seed $seed --prop $prop --DEG $DEG --mixedNPCs
        done
    done
done

