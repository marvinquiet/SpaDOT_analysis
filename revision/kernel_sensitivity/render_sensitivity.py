"""Regenerate summaries with all available prespecified seeds at each point.

Original snapshots and the prior matched-seed presentation remain unchanged.
No training, OT recalculation, or association-correlation analysis is performed.
"""
import os
os.environ['MPLCONFIGDIR']='/tmp/spadot_available_seeds_mpl'
os.environ['XDG_CACHE_HOME']='/tmp/spadot_available_seeds_cache'
from pathlib import Path
import json
import hashlib
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input-dir', type=Path, required=True,
                    help='Audited kernel-sensitivity snapshot root with snapshot/, ot_comparison_snapshot/ and the two baseline summary CSVs')
parser.add_argument('--paper-simulation-ot', type=Path, required=True,
                    help='CSV of paper-parameter per-pair OT scores across the 45 simulation settings')
parser.add_argument('--output-dir', type=Path, required=True)
args = parser.parse_args()
ROOT, OUT = args.input_dir, args.output_dir
OUT.mkdir(parents=True, exist_ok=False)
S=ROOT/'ot_comparison_snapshot'
for folder in ['snapshot','ot_comparison_snapshot']:
    source=ROOT/folder
    for name,expected in json.loads((source/'TRANSFER_SHA256.json').read_text()).items():
        assert hashlib.sha256((source/name).read_bytes()).hexdigest()==expected,name
D=pd.read_csv(S/'domain_metrics_three_real.csv')
T=pd.read_csv(S/'transition_metrics_both_OT.csv')
scales=[.01,.03,.1,.3,1.]
plot_scales=np.array(scales)
seeds={666,1993,2026}
fields=[('ARI','Domain ARI ↑'),('PAS','PAS ↓'),('CHAOS','CHAOS ↓'),('naive_accuracy','Naive transition\naccuracy ↑')]
datasets=[('MouseEmbryo','MouseOrganogenesis','#496b48'),
          ('ChickenHeart','Chicken Cardiogenesis','#236c80'),
          ('AxolotlTelencephalon','Axolotl Telencephalon\nNeurogenesis','#8b4c76')]
assert set(D.seed)==seeds and set(T.seed)==seeds
assert not D.duplicated(['dataset','kernel_scale','seed','timepoint']).any()
assert not T.duplicated(['family','dataset','kernel_scale','seed','i']).any()
simulation=pd.read_csv(ROOT/'final_metric_summary.csv')
simulation=simulation.loc[(simulation.dataset=='all_simulations')&simulation.metric.isin([m for m,_ in fields])].copy()
assert len(simulation)==20
alternative_simulation_path=args.paper_simulation_ot
assert alternative_simulation_path.is_file(), alternative_simulation_path
alternative_simulation_summary=None
if alternative_simulation_path.exists():
    ast=pd.read_csv(alternative_simulation_path)
    arun=ast.groupby(['dataset','kernel_scale','seed'],as_index=False).naive_accuracy.mean()
    aset=arun.groupby(['dataset','kernel_scale'],as_index=False).naive_accuracy.mean()
    assert aset.dataset.nunique()==45
    alternative_simulation_summary=aset.groupby('kernel_scale').naive_accuracy.agg(mean='mean',std='std',count='count').reindex(scales)
    alternative_simulation_summary.to_csv(OUT/'simulation_alternative_naive_summary.csv')
curves={}
records=[]
run_records=[]
for family in ['current','alternative']:
    for ds,_,_ in datasets:
        for metric,_ in fields:
            raw=D.loc[D.dataset==ds] if metric!='naive_accuracy' else T.loc[(T.dataset==ds)&(T.family==family)]
            per=raw.groupby(['kernel_scale','seed'],as_index=False)[metric].mean()
            if metric=='naive_accuracy':
                expected_pairs=3 if ds=='ChickenHeart' else 2
                assert raw.groupby(['kernel_scale','seed']).size().eq(expected_pairs).all()
            for row in per.itertuples(index=False):
                run_records.append(dict(family=family,dataset=ds,kernel_scale=row.kernel_scale,seed=row.seed,metric=metric,value=getattr(row,metric)))
            q=per.groupby('kernel_scale')[metric].agg(['mean','std','count']).reindex(scales)
            curves[(family,ds,metric)]=q
            for scale,row in q.iterrows():
                actual=sorted(per.loc[per.kernel_scale==scale,'seed'].tolist())
                records.append(dict(family=family,dataset=ds,kernel_scale=scale,metric=metric,
                    mean=row['mean'],sd=row['std'],n_seeds=len(actual),included_seeds=','.join(map(str,actual)),
                    policy='All available prespecified seeds at this coefficient; exclude failed cases only'))
summary=pd.DataFrame(records)
summary.to_csv(OUT/'real_metric_summary_both_OT.csv',index=False)
pd.DataFrame(run_records).to_csv(OUT/'real_per_seed_metrics_both_OT.csv',index=False)
primary=summary.loc[summary.family=='current'].copy()
all_summary=pd.concat([primary.rename(columns={'n_seeds':'n_summary_units'}).assign(summary_unit='training seed'),
    simulation.assign(family='current',summary_unit='simulation setting')],ignore_index=True)
all_summary[['dataset','kernel_scale','metric','mean','sd','n_summary_units','summary_unit']].to_csv(OUT/'final_metric_summary.csv',index=False)
ari=primary.loc[primary.metric=='ARI'].pivot(index='kernel_scale',columns='dataset',values='mean').reindex(scales)
ari.to_csv(OUT/'mean_ARI_with_mouse.csv')
counts=summary[['family','dataset','kernel_scale','metric','n_seeds','included_seeds']]
counts.to_csv(OUT/'point_sample_counts.csv',index=False)
old=pd.read_csv(ROOT/'real_metric_summary_both_OT.csv')
old=old.loc[old.metric.isin([m for m,_ in fields])]
changes=summary.merge(old[['family','dataset','kernel_scale','metric','mean','sd','n_seeds']],on=['family','dataset','kernel_scale','metric'],suffixes=('_new','_previous'),validate='one_to_one')
changes['mean_change']=changes.mean_new-changes.mean_previous
changes.to_csv(OUT/'change_from_previous_matched_summary.csv',index=False)

# Full-grid assertions: only the documented failures may reduce n.
for row in summary.itertuples():
    expected=3
    if row.dataset=='MouseEmbryo' and row.kernel_scale==.03:
        expected=0 if row.family=='current' and row.metric=='naive_accuracy' else 2
    if row.dataset=='ChickenHeart' and row.kernel_scale==.1 and row.family=='current' and row.metric=='naive_accuracy':
        expected=2
    assert row.n_seeds==expected,(row,expected)
assert len(D)==147 and len(D.groupby(['dataset','kernel_scale','seed']))==44
assert len(T.loc[T.family=='current'].groupby(['dataset','kernel_scale','seed']))==41
assert len(T.loc[T.family=='alternative'].groupby(['dataset','kernel_scale','seed']))==44

plt.rcParams.update({'font.size':18,'axes.labelsize':20,'xtick.labelsize':18,'ytick.labelsize':18,'pdf.fonttype':42})
def axis(ax,first,last,title,label):
    ax.set_xscale('log')
    ax.set_xticks(plot_scales,['.01','.03','.1','.3','1'])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(plot_scales[0]*.88,plot_scales[-1]*1.12)
    ax.axvline(.1,color='#999999',linestyle='--',linewidth=1)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',alpha=.12)
    if first:ax.set_title(title,fontsize=19,pad=12)
    if last:ax.set_xlabel(r'Kernel $l$',fontsize=20)
    if label:ax.set_ylabel(label,fontsize=20)

def draw(family,with_sim,stem):
    items=datasets+([('all_simulations','Simulation Study','#a67523')] if with_sim else [])
    fig,axs=plt.subplots(len(items),4,figsize=(18,3.5*len(items)+.2),layout='constrained')
    for i,(ds,label,color) in enumerate(items):
        for j,(metric,title) in enumerate(fields):
            ax=axs[i,j]
            if ds=='all_simulations':
                q=simulation.loc[simulation.metric==metric].set_index('kernel_scale').reindex(scales).rename(columns={'sd':'std'})
                if family=='alternative' and metric=='naive_accuracy':
                    q=alternative_simulation_summary if alternative_simulation_summary is not None else pd.DataFrame({'mean':[np.nan]*5,'std':[np.nan]*5},index=scales)
            else:q=curves[(family,ds,metric)]
            ax.errorbar(plot_scales,q['mean'],yerr=q['std'].fillna(0),marker='o',markersize=3.2,
                        markerfacecolor='white',markeredgewidth=1.1,barsabove=True,
                        capsize=5,capthick=1.4,elinewidth=1.5,lw=1.7,color=color)
            axis(ax,i==0,i==len(items)-1,title,label if j==0 else '')
            if ds=='all_simulations' and q['mean'].isna().all():
                ax.set_ylim(0,1)
                ax.set_yticks([0,.5,1])
                ax.text(.5,.5,'Alternative OT\nin progress',ha='center',va='center',transform=ax.transAxes,fontsize=16,color='#9b3636')
            if ds!='all_simulations':
                for scale,row in q.iterrows():
                    if pd.notna(row['mean']) and row['count']<3:
                        ax.annotate(f'n={int(row["count"])}',(scale,row['mean']+row['std']),xytext=(0,9),textcoords='offset points',ha='center',fontsize=12)
                lo=(q['mean']-q['std'].fillna(0)).min()
                hi=(q['mean']+q['std'].fillna(0)).max()
                gap=max(hi-lo,.00001)
                ax.set_ylim(lo-.12*gap,hi+.28*gap)
                if q['mean'].isna().any():
                    ax.text(.04,.05,'OT unavailable\nat 0.03',transform=ax.transAxes,va='bottom',fontsize=16,color='#9b3636')
    for ext in ['png','pdf']:fig.savefig(OUT/f'{stem}.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)
draw('current',True,'kernel_sensitivity_final')
draw('alternative',True,'kernel_sensitivity_alternative_OT_real')

# Comparison uses all available seeds separately under each OT setting.
# No complete-case restriction across OT families is applied to these curves.
fig,axs=plt.subplots(3,1,figsize=(8.5,12),layout='constrained')
for i,(ds,label,_) in enumerate(datasets):
    for family,color,legend in [('current','#236c80','λ₁=0.1, λ₂=5, ε=0.02'),('alternative','#d77e24','λ₁=1, λ₂=50, ε=0.01')]:
        q=curves[(family,ds,'naive_accuracy')]
        axs[i].errorbar(plot_scales,q['mean'],yerr=q['std'].fillna(0),marker='o',markersize=3.2,
                       markerfacecolor='white',markeredgewidth=1.1,barsabove=True,
                       capsize=5,capthick=1.4,elinewidth=1.5,color=color,label=legend,lw=1.7)
        for scale,row in q.iterrows():
            if pd.notna(row['mean']) and row['count']<3:
                axs[i].annotate(f'n={int(row["count"])}',(scale,row['mean']),xytext=(0,12 if family=='alternative' else -19),textcoords='offset points',ha='center',fontsize=12,color=color)
    axis(axs[i],False,i==2,'','Naive accuracy')
    axs[i].set_title(label.replace('\n',' '),fontsize=17)
    if i==0:axs[i].legend(frameon=False,fontsize=12)
for ext in ['png','pdf']:fig.savefig(OUT/f'kernel_sensitivity_OT_comparison_real.{ext}',dpi=180,bbox_inches='tight')
plt.close(fig)

# Pair-specific naive accuracy with all available seeds; no correlation metric.
pairs=T.groupby(['family','dataset','kernel_scale','i']).naive_accuracy.agg(mean='mean',sd='std',n='count').reset_index()
pairs.to_csv(OUT/'naive_accuracy_by_timepoint_pair.csv',index=False)
T.to_csv(OUT/'naive_accuracy_per_seed_pair.csv',index=False)
fig,axs=plt.subplots(1,3,figsize=(18,5),layout='constrained')
for ax,(ds,label,_) in zip(axs,datasets):
    for pair,g in pairs.loc[(pairs.family=='current')&(pairs.dataset==ds)].groupby('i'):
        q=g.set_index('kernel_scale').reindex(scales)
        ax.errorbar(plot_scales,q['mean'],yerr=q.sd,marker='o',markersize=3.2,
                    markerfacecolor='white',barsabove=True,capsize=5,capthick=1.4,
                    label=f't{int(pair)} → t{int(pair)+1}')
    axis(ax,False,True,'','Naive accuracy')
    ax.set_title(label,fontsize=17)
    ax.legend(frameon=False,fontsize=12)
fig.savefig(OUT/'naive_accuracy_by_timepoint_pair.png',dpi=180,bbox_inches='tight')
plt.close(fig)

body='''<!doctype html><meta charset="utf-8"><title>Kernel sensitivity: all available seeds</title>
<style>body{font:16px system-ui;line-height:1.5;max-width:1600px;margin:35px auto;padding:20px;color:#25333d}img{max-width:100%}th,td{padding:7px 12px;border-bottom:1px solid #ddd}</style>
<h1>Kernel sensitivity: all available prespecified seeds</h1>
<p>Updated inclusion policy: at each kernel coefficient and for each metric, include every available result from training seeds666,1993,2026. Exclude only failed/unavailable cases; do not discard a successful result because the same seed failed elsewhere. Association correlation is excluded.</p>
<p>All real domain points have n=3 except Mouse coefficient0.03 (n=2, missing training seed1993). Current-OT points have n=3 except Chicken coefficient0.1 (n=2) and Mouse coefficient0.03 (n=0). Alternative-OT points have n=3 except Mouse coefficient0.03 (n=2). Reduced nonzero counts are annotated on the figures; all counts appear in the tables.</p>
<p>Real points show mean ± sample SD over available training seeds, after averaging timepoints or adjacent pairs within each seed. Simulation points average the three training seeds within each of45 settings, then show mean ± sample SD across settings. Error bars are not confidence intervals; no hypothesis tests were used.</p>
<p>The x-axis is Kernel l, where l denotes the kernel_scale denominator coefficient used in the implementation: K_l(i,j)=exp(-distance(s_i,s_j)^2/l). The plotted values are0.01,0.03,0.1,0.3,1; the dashed line marks the unchanged implementation default l=0.1. No square-root reparameterization is applied to the axis. Hollow smaller markers and wider caps expose small SD bars without altering the data. Fixed inputs,99 executed training epochs and the original saved embeddings are preserved.</p>
<h2>Primary OT: λ₁=0.1, λ₂=5, ε=0.02</h2><img src="kernel_sensitivity_final.png">
<h2>Alternative OT: λ₁=1, λ₂=50, ε=0.01</h2><img src="kernel_sensitivity_alternative_OT_real.png">
<p>The simulation domain metrics are unchanged by post-training OT. The simulation naive-accuracy panel uses the alternative-OT results when available; an unevaluated panel is explicitly marked rather than replaced with current-OT values.</p>
<h2>OT comparison with all available results under each setting</h2><img src="kernel_sensitivity_OT_comparison_real.png">
<p>These are descriptive available-case means. The two OT curves can have different seed counts at a coefficient, as disclosed above; they are not a paired statistical test.</p>
<h2>Pair-specific naive accuracy</h2><img src="naive_accuracy_by_timepoint_pair.png">
<p>The naive score uses annotation-aggregated spot transport, entrywise minimum row/column normalization, and the original benchmark correspondence graph. It applies no visualization or curated-score threshold.</p>
<h2>Primary mean ARI</h2>'''+ari.to_html(float_format=lambda x:f'{x:.6f}')+ '<h2>Per-point counts</h2>'+counts.to_html(index=False)+'''
<p><a href="kernel_sensitivity_final.pdf">Primary figure PDF</a> · <a href="kernel_sensitivity_alternative_OT_real.pdf">Alternative figure PDF</a> · <a href="kernel_sensitivity_OT_comparison_real.pdf">OT comparison PDF</a></p>'''
(OUT/'report.html').write_text(body)

verification=dict(inclusion_policy='All available seeds666,1993,2026 at each coefficient and metric; exclude failed cases only',
    source_transfer_hashes_verified=True,original_snapshots_and_prior_figures_preserved=True,
    valid_training_runs=719,primary_real_domain_runs=44,current_real_OT_runs=41,alternative_real_OT_runs=44,
    new_training_runs=0,alternative_simulation_OT_scores_available=alternative_simulation_summary is not None,association_r_excluded=True,
    axis_label='Kernel l',axis_units='implementation kernel_scale denominator coefficient; no conversion',
    markers='hollow3.2pt; errorbar cap5pt; unchanged SD',
    no_global_complete_case_filter=True,simulation_results_unchanged=True)
(OUT/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
print('PRIMARY DOMAIN ARI\n'+ari.round(6).to_string())
print('PRIMARY NAIVE ACCURACY\n'+primary.loc[primary.metric=='naive_accuracy',['dataset','kernel_scale','mean','sd','n_seeds']].round(6).to_string(index=False))
print('OUTPUT',OUT)
