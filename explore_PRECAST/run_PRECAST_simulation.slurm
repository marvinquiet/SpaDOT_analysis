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
# cd /net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_PRECAST

declare -a arr=(1993 2024 1105 327 113 905 789 1314 888 666)
for seed in "${arr[@]}"
do
    echo "=== Run PRECAST on simulation data for $seed ==="
    # Rscript run_PRECAST_simulation.R -x 20 -y 50 -s $seed
    # #Rscript run_PRECAST_simulation.R -x 10 -y 50 -s $seed
    # Rscript run_PRECAST_simulation.R -x 20 -y 30 -s $seed
    # #Rscript run_PRECAST_simulation.R -x 10 -y 30 -s $seed
    # Rscript run_PRECAST_simulation.R -x 30 -y 50 -s $seed
    # Rscript run_PRECAST_simulation.R -x 30 -y 30 -s $seed

    Rscript run_PRECAST_simulation_prop.R -s $seed --prop 0.9
    Rscript run_PRECAST_simulation_prop.R -s $seed --prop 0.8
    Rscript run_PRECAST_simulation_prop.R -s $seed --prop 0.7
    Rscript run_PRECAST_simulation_prop.R -s $seed --prop 0.6
    Rscript run_PRECAST_simulation_prop.R -s $seed --prop 0.5
done
