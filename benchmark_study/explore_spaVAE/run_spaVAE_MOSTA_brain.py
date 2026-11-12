# conda activate /media/sda/wenjinma/envs/spatialATAC
import os,sys
import math
import h5py
import anndata
import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib.pyplot as plt
from matplotlib import transforms
import seaborn as sns
import colorcet as cc
sns.set_palette(cc.glasbey)
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_MOSTA"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spaVAE/MOSTA_brain_results"

sys.path.append('/media/sda/wenjinma/projects/spatialATAC-explore/explore_spaVAE/spaVAE/src/spaVAE_Batch/')
# run spaVAE Batch version
os.system('python /media/sda/wenjinma/projects/spatialATAC-explore/explore_spaVAE/spaVAE/src/spaVAE_Batch/run_spaVAE_Batch.py --data_file {0}/combined.h5 --model_file {0}/model.pt --final_latent_file {0}/final_latent.txt --denoised_counts_file {0}/denoised_counts.txt > {0}/log.out 2>&1'.format(result_dir))  