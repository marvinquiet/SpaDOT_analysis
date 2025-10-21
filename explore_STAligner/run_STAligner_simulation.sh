#!/bin/bash
source ~/.bashrc
conda activate /media/sda/wenjinma/envs/STAligner_py38
cd /media/sda/wenjinma/projects/spatialATAC-explore/explore_STAligner

declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
for seed in "${arr[@]}"
do
    python run_STAligner_simulation.py --seed $seed --DEG_x 20 --DEG_y 50
    python run_STAligner_simulation.py --seed $seed --DEG_x 30 --DEG_y 50
    python run_STAligner_simulation.py --seed $seed --DEG_x 20 --DEG_y 30
    python run_STAligner_simulation.py --seed $seed --DEG_x 30 --DEG_y 30
done
