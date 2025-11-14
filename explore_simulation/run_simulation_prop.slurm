#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --partition=main
#SBATCH --mem=100G
#SBATCH --job-name=simulated_data

### ====
## This script generate simulated brain data and select SVGs
### ===

source ~/.bashrc
conda activate ~/envs/spatialATAC
cd /net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_simulation

#DEG=200
#DEG=500

#declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
#declare -a arr=(1993 2024 1105)
#for seed in "${arr[@]}"
#do
#    echo "=== Run simulation for $seed ==="
#    # --- simulation for cell type proportions ---
#    Rscript simulate_brain_prop.R --prop 0.9 --seed $seed --DEG $DEG
#    python simulate_brain_prop_generateh5ad.py --prop 0.9 --seed $seed --DEG $DEG
#    Rscript simulate_brain_prop.R --prop 0.8 --seed $seed --DEG $DEG   
#    python simulate_brain_prop_generateh5ad.py --prop 0.8 --seed $seed --DEG $DEG
#    Rscript simulate_brain_prop.R --prop 0.7 --seed $seed --DEG $DEG   
#    python simulate_brain_prop_generateh5ad.py --prop  0.7 --seed $seed --DEG $DEG
#    Rscript simulate_brain_prop.R --prop 0.6 --seed $seed --DEG $DEG   
#    python simulate_brain_prop_generateh5ad.py --prop 0.6 --seed $seed --DEG $DEG
#    Rscript simulate_brain_prop.R --prop 0.5 --seed $seed --DEG $DEG   
#    python simulate_brain_prop_generateh5ad.py --prop 0.5 --seed $seed --DEG $DEG
#done

declare -a arr=(1993 2024 1105)
declare -a prop_list=(0.9 0.8 0.7 0.6 0.5)
declare -a DEG_list=(100 200 500)
for seed in "${arr[@]}"; do
    for prop in "${prop_list[@]}"; do
        for DEG in "${DEG_list[@]}"; do
            Rscript simulate_brain_prop_mixedNPCs.R --prop $prop --seed $seed --DEG $DEG
            python simulate_brain_prop_mixedNPCs_generateh5ad.py --prop $prop --seed $seed --DEG $DEG
        done
    done
done
