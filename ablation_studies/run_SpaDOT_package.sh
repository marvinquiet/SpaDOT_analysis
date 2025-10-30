#!/bin/bash

for i in {1..50}; do
    python run_SpaDOT.py --data preprocessed_MouseOrganogenesis.h5ad --prefix MouseOrganogenesis
    #python run_SpaDOT.py --data preprocessed_ChickenHeart.h5ad --prefix ChickenHeart
    #python run_SpaDOT.py --data preprocessed_Axolotl.h5ad --prefix Axolotl
done
