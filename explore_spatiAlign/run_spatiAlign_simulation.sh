#!/bin/bash
# source ~/.bashrc
# conda activate /media/sda/wenjinma/envs/STAligner_py38
# cd /media/sda/wenjinma/projects/spatialATAC-explore/explore_spatiAlign

# conda activate /net/zootopia/disk1/wenjinma/envs/spatialign_py38
# cd /net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_spatiAlign

# declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
declare -a arr=(1993 2024 1105)
declare -a prop_list=(0.9 0.8 0.7 0.6 0.5)
declare -a DEG_list=(100 200 500)
for seed in "${arr[@]}"; do
    for prop in "${prop_list[@]}"; do
        for DEG in "${DEG_list[@]}"; do
            # python run_spatiAlign_simulation_prop.py --seed $seed --prop $prop --DEG $DEG
            python run_spatiAlign_simulation_prop.py --seed $seed --prop $prop --DEG $DEG --mixedNPCs
        done
    done
done
