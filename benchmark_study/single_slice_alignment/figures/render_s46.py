"""Render Fig. S46 from aggregate_s46.py output.

Usage: python render_s46.py --input aggregated_curves.json --output-dir DIR
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='1'
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--output-dir', type=Path, required=True)
args = parser.parse_args()
SOURCE, OUT = args.input, args.output_dir
curves=json.loads(SOURCE.read_text()); assert len(curves)==84
lookup={(c['deg'],c['method'],c['transition'],c['mode']):c for c in curves}
methods=['PRECAST','SEDR','spaGCN','GraphST','spatialPCA','spatiAlign','spaVAE']
colors=dict(zip(methods,['#68a756','#f09232','#4c2623','#165451','#d3b300','#ded6c1','#0082c0']))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':18,'axes.titlesize':23,'pdf.fonttype':42,'ps.fonttype':42})
fig,axes=plt.subplots(7,6,figsize=(26,20))
fig.subplots_adjust(left=.078,right=.99,top=.928,bottom=.073,hspace=.57,wspace=.24)
for row,m in enumerate(methods):
 for di,deg in enumerate([100,200,500]):
  for ti,tr in enumerate(['t0_t1','t1_t2']):
   ax=axes[row,2*di+ti]
   for mode,style in [('multi','-'),('single','--')]:
    c=lookup[(deg,m,tr,mode)];x=np.asarray(c['distance'])
    ax.fill_between(x,c['lower_95'],c['upper_95'],facecolor=colors[m],alpha=.18,linewidth=0,zorder=1)
    ax.plot(x,c['mean_cdf'],color=colors[m],lw=2.8,ls=style,zorder=3)
   ax.set_xscale('log');ax.set_xlim(x[0],x[-1]);ax.set_ylim(0,1.015);ax.set_yticks([0,.5,1]);ax.grid(alpha=.15,zorder=0)
   ax.tick_params(axis='both',labelsize=16,pad=4)
   for side in ['top','right']:ax.spines[side].set_visible(False)
   if di==0 and ti==0:ax.set_ylabel(m,fontweight='bold',fontsize=23,labelpad=14)
   else:ax.tick_params(labelleft=False)
   if row==0:ax.set_title(tr.replace('_',' → '),fontsize=23,pad=16)
for di,deg in enumerate([100,200,500]):
 a=axes[0,di*2].get_position();b=axes[0,di*2+1].get_position()
 fig.text((a.x0+b.x1)/2,.982,f'DEG {deg}',ha='center',va='top',fontsize=32,fontweight='bold')
fig.supylabel('Cumulative fraction',fontsize=23,x=.010)
legend=[Line2D([],[],color='black',lw=3.1,ls='-',label='Multi-slice'),Line2D([],[],color='black',lw=3.1,ls='--',label='Single-slice'),Patch(facecolor='#777777',alpha=.2,edgecolor='none',label='95% pointwise CI (3 seeds)')]
fig.legend(handles=legend,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.5,.016),fontsize=21,handlelength=2.6,columnspacing=2.3)
OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT/'S46_Simulation_CDF_aggregated_DEG_CI.pdf')
plt.close(fig)
