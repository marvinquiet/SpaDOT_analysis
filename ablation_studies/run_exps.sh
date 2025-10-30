#!/bin/bash

# --- combined
#project_dir="SpaDOT" # original SpaDOT
#project_dir="SpaDOT_VAE_combined" # VAE
#project_dir="SpaDOT_SVGPVAE_combined" # SVGP-VAE
#project_dir="SpaDOT_GCNVAE_combined" # GCN-VAE
#project_dir="SpaDOT_GATVAE_combined" # GAT-VAE
#project_dir="SpaDOT_SVGPGATVAE_combined" # SVGP+GAT-VAE
#project_dir="SpaDOT_SVGPGATVAE_alignment_combined" # SVGP+GAT-VAE
#project_dir="SpaDOT_SVGPGATVAE_alignment_KMeans_combined" # SVGP+GAT-VAE
#project_dir="SpaDOT_full" # SVGP+GAT-VAE
#project_dir="SpaDOT_full_noalignment" # experiment

# --- separate
#project_dir="SpaDOT_VAE_separate"
#project_dir="SpaDOT_SVGPVAE_separate"
#project_dir="SpaDOT_GCNVAE_separate"
#project_dir="SpaDOT_GATVAE_separate"
#project_dir="SpaDOT_SVGPGATVAE_separate"
#project_dir="SpaDOT_SVGPGATVAE_alignment_separate"
#project_dir="SpaDOT_SVGPGATVAE_alignment_KMeans_separate"

#python $project_dir/train.py --data preprocessed_MouseOrganogenesis.h5ad --output_dir $project_dir/output/ --prefix MouseOrganogenesis
#python $project_dir/train.py --data preprocessed_ChickenHeart.h5ad --output_dir $project_dir/output/ --prefix ChickenHeart
#python $project_dir/train.py --data preprocessed_Axolotl.h5ad --output_dir $project_dir/output/ --prefix Axolotl



#python analyze.py --output_dir $project_dir/output/ --prefix MouseOrganogenesis1993
#
#
#
#python $project_dir/train.py --data preprocessed_ChickenHeart.h5ad --output_dir $project_dir/output/ --prefix ChickenHeart
#python analyze.py --output_dir $project_dir/output/ --prefix ChickenHeart1993
#
#python $project_dir/train.py --data preprocessed_Axolotl.h5ad --output_dir $project_dir/output/ --prefix Axolotl
#python analyze.py --output_dir $project_dir/output/ --prefix Axolotl1993

