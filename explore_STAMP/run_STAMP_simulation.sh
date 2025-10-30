#!/bin/bash
#source ~/.bashrc
# conda activate /net/zootopia/disk1/wenjinma/envs/stamp
#cd /net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_STAMP

#declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
declare -a arr=(1993)
for seed in "${arr[@]}"
do
    #python run_STAMP_simulation.py --seed $seed --DEG_x 20 --DEG_y 50
    #python run_STAMP_simulation.py --seed $seed --DEG_x 30 --DEG_y 50
    #python run_STAMP_simulation.py --seed $seed --DEG_x 20 --DEG_y 30
    #python run_STAMP_simulation.py --seed $seed --DEG_x 30 --DEG_y 30
    python run_STAMP_simulation_prop.py --seed $seed --prop 0.9
    python run_STAMP_simulation_prop.py --seed $seed --prop 0.8
    python run_STAMP_simulation_prop.py --seed $seed --prop 0.7
    python run_STAMP_simulation_prop.py --seed $seed --prop 0.6
    python run_STAMP_simulation_prop.py --seed $seed --prop 0.5
done
