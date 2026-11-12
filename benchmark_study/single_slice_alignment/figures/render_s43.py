"""Render Fig. S43 from matched-spot real CDF points produced by prepare_figures.py.

Usage: python render_s43.py --input CSV --output-dir DIR
"""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--output-dir', type=Path, required=True)
args = parser.parse_args()
OUT = args.output_dir
OUT.mkdir(parents=True, exist_ok=True)
curves = pd.read_csv(args.input)
METHODS=['PRECAST','SEDR','spaGCN','GraphST','spatialPCA','spatiAlign','spaVAE']
COLORS={'PRECAST':'#68a756','SEDR':'#f09232','spaGCN':'#4c2623','GraphST':'#165451','spatialPCA':'#d3b300','spatiAlign':'#ded6c1','spaVAE':'#0082c0'}
groups=[('MouseEmbryo','Mouse Organogenesis',['E12.5','E13.5','E15.5']),('ChickenHeart','Chicken Cardiogenesis',['D4','D7','D10','D14']),('AxolotlTelencephalon','Axolotl Telencephalon\nNeurogenesis',['St.44','St.54','St.57'])]
columns=[(ds,stages[i]+' to '+stages[i+1],i) for ds,_,stages in groups for i in range(len(stages)-1)]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14,'pdf.fonttype':42,'axes.linewidth':.8})
fig,axes=plt.subplots(7,7,figsize=(29,20),squeeze=False,sharey=True)
fig.subplots_adjust(left=.062,right=.99,bottom=.08,top=.90,hspace=.47,wspace=.17)
boxes=[[ax.get_position().frozen() for ax in row] for row in axes]
spatial_segments=[(1,200),(1e5,1e12)]
axis_records=[]

def draw(ax,sub,method,limits,ticks,show_y):
    for mode,style in [('Single-slice','--'),('Multi-slice','-')]:
        g=sub[sub['mode']==mode]
        x=g.distance.to_numpy();y=g.cdf.to_numpy()
        assert np.all(np.diff(x)>=0) and np.all(np.diff(y)>=0)
        # Flat ECDF tails are displayed; original distance and CDF values are untouched.
        xx=np.r_[min(limits[0],x[0]/10),x,max(limits[1],x[-1]*10)]
        yy=np.r_[0,y,1]
        line,=ax.plot(xx,yy,color=COLORS[method],ls=style,lw=2.35)
        if method=='spatiAlign':line.set_path_effects([pe.Stroke(linewidth=3.4,foreground='#888888'),pe.Normal()])
    ax.set_xscale('log');ax.set_xlim(*limits);ax.set_ylim(0,1.015)
    ax.xaxis.set_major_locator(ticker.FixedLocator(ticks))
    ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_yticks([0,.25,.5,.75,1]);ax.tick_params(labelsize=12)
    ax.tick_params(axis='y',labelleft=show_y)
    ax.grid(True,which='major',alpha=.18)
    ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)

for row,method in enumerate(METHODS):
    for col,(ds,title,idx) in enumerate(columns):
        sub=curves[(curves.method==method)&(curves.dataset==ds)&(curves.transition_index==idx)]
        assert set(sub['mode'])=={'Single-slice','Multi-slice'}
        ax=axes[row,col]
        if method=='spatialPCA':
            b=boxes[row][col];ax.remove()
            gap=b.width*.085;lw=b.width*.43;rw=b.width-lw-gap
            left=fig.add_axes([b.x0,b.y0,lw,b.height]);right=fig.add_axes([b.x0+lw+gap,b.y0,rw,b.height])
            draw(left,sub,method,spatial_segments[0],[1,10,100],col==0)
            draw(right,sub,method,spatial_segments[1],[1e5,1e8,1e11],False)
            right.spines['left'].set_visible(False);right.tick_params(axis='y',left=False)
            if col==0:left.set_ylabel(method,fontsize=19,fontweight='bold',labelpad=16)
            # A pair of slashes at the gap denotes the omitted distance interval.
            center=b.x0+lw+gap/2
            for offset in [-.0018,.0018]:
                fig.add_artist(Line2D([center+offset-.0012,center+offset+.0012],[b.y0-.004,b.y0+.004],transform=fig.transFigure,color='black',lw=1.2,clip_on=False))
            axis_records.append({'method':method,'dataset':ds,'transition_index':idx,'segments':spatial_segments})
        else:
            draw(ax,sub,method,(1,100),[1,10,100],col==0)
            if col==0:ax.set_ylabel(method,fontsize=19,fontweight='bold',labelpad=16)
            axis_records.append({'method':method,'dataset':ds,'transition_index':idx,'xlim':list(ax.get_xlim())})
        if row==0:ax.set_title(title,fontsize=17,pad=14)
start=0
for ds,title,stages in groups:
    n=len(stages)-1
    left=boxes[0][start].x0;right=boxes[0][start+n-1].x1
    fig.text((left+right)/2,.954,title,ha='center',va='center',fontsize=23,fontweight='bold',linespacing=1.1)
    start+=n
for last in [1,4]:
    x=(boxes[0][last].x1+boxes[0][last+1].x0)/2
    fig.add_artist(Line2D([x,x],[.075,.918],transform=fig.transFigure,color='#b6b6b6',lw=1,ls=':'))
fig.legend(handles=[Line2D([0],[0],color='#333333',lw=2.5,ls='--',label='Single-slice'),Line2D([0],[0],color='#333333',lw=2.5,ls='-',label='Multi-slice')],loc='lower center',bbox_to_anchor=(.525,.014),ncol=2,frameon=False,fontsize=20)
fig.text(.009,.49,'Cumulative fraction',rotation=90,ha='center',va='center',fontsize=19)
widths=[b.width for row in boxes for b in row]
assert max(widths)-min(widths)<1e-10
assert all(r.get('xlim')==[1,100] for r in axis_records if r['method']!='spatialPCA')
sp=curves[curves.method=='spatialPCA']
assert not ((sp.distance>200)&(sp.distance<1e5)).any()
fig.savefig(OUT/'S43_raw_distance_CDF_combined.pdf')
fig.savefig(OUT/'Raw_distance_CDF_combined.png',dpi=110)
plt.close(fig)
