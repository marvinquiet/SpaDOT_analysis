#!/bin/bash
source ~/.bashrc
conda activate /media/sda/wenjinma/envs/GraphST
cd /media/sda/wenjinma/projects/spatialATAC-explore/explore_SPIRAL

declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
for seed in "${arr[@]}"
do
    python run_SPIRAL_simulation.py --seed $seed --DEG_x 20 --DEG_y 50
    python run_SPIRAL_simulation.py --seed $seed --DEG_x 30 --DEG_y 50
    python run_SPIRAL_simulation.py --seed $seed --DEG_x 20 --DEG_y 30
    python run_SPIRAL_simulation.py --seed $seed --DEG_x 30 --DEG_y 30
done
