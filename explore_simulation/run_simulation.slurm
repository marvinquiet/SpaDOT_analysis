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

declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
for seed in "${arr[@]}"
do
    echo "=== Run simulation for $seed ==="
    # --- simulation for DEGs ---
    Rscript simulate_brain.R -x 20 -y 50 -s $seed
    python simulate_brain_generateh5ad.py --DEG_x 20 --DEG_y 50 --seed $seed 
    Rscript simulate_brain.R -x 20 -y 30 -s $seed
    python simulate_brain_generateh5ad.py --DEG_x 20 --DEG_y 30 --seed $seed 
    Rscript simulate_brain.R -x 30 -y 50 -s $seed
    python simulate_brain_generateh5ad.py --DEG_x 30 --DEG_y 50 --seed $seed 
    Rscript simulate_brain.R -x 30 -y 30 -s $seed
    python simulate_brain_generateh5ad.py --DEG_x 30 --DEG_y 30 --seed $seed 
    #echo "=== Select SVGs"
    #Rscript select_SVGs.R --scenario_dir /net/mulan/home/wenjinma/projects/spatialATAC-explore/results/spatialPCA_simulation_brain/${seed}/10_30
    #Rscript select_SVGs.R --scenario_dir /net/mulan/home/wenjinma/projects/spatialATAC-explore/results/spatialPCA_simulation_brain/${seed}/10_50
    #Rscript select_SVGs.R --scenario_dir /net/mulan/home/wenjinma/projects/spatialATAC-explore/results/spatialPCA_simulation_brain/${seed}/20_30
    #Rscript select_SVGs.R --scenario_dir /net/mulan/home/wenjinma/projects/spatialATAC-explore/results/spatialPCA_simulation_brain/${seed}/20_50
    #Rscript select_SVGs.R --scenario_dir /net/mulan/home/wenjinma/projects/spatialATAC-explore/results/spatialPCA_simulation_brain/${seed}/30_30
    #Rscript select_SVGs.R --scenario_dir /net/mulan/home/wenjinma/projects/spatialATAC-explore/results/spatialPCA_simulation_brain/${seed}/30_50

    # --- simulation for cell type proportions ---
    #Rscript simulate_brain_prop.R --prop 0.9 --seed $seed
    #python simulate_brain_prop_generateh5ad.py --prop 0.9 --seed $seed 
    #Rscript simulate_brain_prop.R --prop 0.8 --seed $seed
    #python simulate_brain_prop_generateh5ad.py --prop 0.8 --seed $seed 
    #Rscript simulate_brain_prop.R --prop 0.7 --seed $seed
    #python simulate_brain_prop_generateh5ad.py --prop  0.7 --seed $seed 
    #Rscript simulate_brain_prop.R --prop 0.6 --seed $seed
    #python simulate_brain_prop_generateh5ad.py --prop 0.6 --seed $seed 
    #Rscript simulate_brain_prop.R --prop 0.5 --seed $seed
    #python simulate_brain_prop_generateh5ad.py --prop 0.5 --seed $seed 
done
