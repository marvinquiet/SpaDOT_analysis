#!/bin/bash
#source ~/.bashrc
#conda activate /media/sda/wenjinma/envs/spatialATAC
#cd /media/sda/wenjinma/projects/spatialATAC-explore/explore_spaVAE

# declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
# for seed in "${arr[@]}"
# do
#     # prepare simulation data
#     # #python prepare_simulation.py --seed $seed --DEG_x 10 --DEG_y 30
#     # python prepare_simulation.py --seed $seed --DEG_x 20 --DEG_y 30
#     # python prepare_simulation.py --seed $seed --DEG_x 30 --DEG_y 30
#     # #python prepare_simulation.py --seed 666 --DEG_x 10 --DEG_y 50
#     # python prepare_simulation.py --seed $seed --DEG_x 20 --DEG_y 50
#     # python prepare_simulation.py --seed $seed --DEG_x 30 --DEG_y 50
    
#     # run spaVAE
#     ##python run_spaVAE_simulation.py --seed 666 --DEG_x 10 --DEG_y 30
#     #python run_spaVAE_simulation.py --seed $seed --DEG_x 20 --DEG_y 30
#     #python run_spaVAE_simulation.py --seed $seed --DEG_x 30 --DEG_y 30
#     ##python run_spaVAE_simulation.py --seed 666 --DEG_x 10 --DEG_y 50
#     #python run_spaVAE_simulation.py --seed $seed --DEG_x 20 --DEG_y 50
#     #python run_spaVAE_simulation.py --seed $seed --DEG_x 30 --DEG_y 50
# done

declare -a arr=(1993 2024 1105)
declare -a prop_list=(0.9 0.8 0.7 0.6 0.5)
declare -a DEG_list=(100 200 500)
for seed in "${arr[@]}"; do
    for prop in "${prop_list[@]}"; do
        for DEG in "${DEG_list[@]}"; do
            # python prepare_simulation_prop.py --seed $seed --prop $prop --DEG $DEG
            # python run_spaVAE_simulation_prop.py --seed $seed --prop $prop --DEG $DEG
            python prepare_simulation_prop.py --seed $seed --prop $prop --DEG $DEG --mixedNPCs
            python run_spaVAE_simulation_prop.py --seed $seed --prop $prop --DEG $DEG --mixedNPCs
        done
    done
done
