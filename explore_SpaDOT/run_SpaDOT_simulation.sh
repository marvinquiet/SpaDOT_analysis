#!/bin/bash
# conda activate ~/envs/SpaDOT
# cd /media/sda/wenjinma/projects/spatialATAC-explore/explore_SpaDOT

#declare -a arr=(1993)
declare -a arr=(2024 1105 327 113 905 789 1314 888 666)
for seed in "${arr[@]}"
do
    #python run_SpaDOT_simulation.py --seed $seed --DEG_x 20 --DEG_y 50
    #python run_SpaDOT_simulation.py --seed $seed --DEG_x 30 --DEG_y 50
    #python run_SpaDOT_simulation.py --seed $seed --DEG_x 20 --DEG_y 30
    #python run_SpaDOT_simulation.py --seed $seed --DEG_x 30 --DEG_y 30
    python run_SpaDOT_simulation_prop.py --seed $seed --prop 0.9
    python run_SpaDOT_simulation_prop.py --seed $seed --prop 0.8
    python run_SpaDOT_simulation_prop.py --seed $seed --prop 0.7
    python run_SpaDOT_simulation_prop.py --seed $seed --prop 0.6
    python run_SpaDOT_simulation_prop.py --seed $seed --prop 0.5
done
