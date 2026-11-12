"""Prepare sensitivity statistics for the original epsilon=.02 OT protocol.

Use render_sensitivity.py and its alternative-family output for paper settings.
The snapshot root must contain the audited intermediate metric tables. The
results root must contain the all-available seed real_metric_summary_both_OT.csv.
"""
import argparse
from pathlib import Path
import json
import hashlib
import math
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--snapshot-root', required=True, type=Path)
parser.add_argument('--results-root', required=True, type=Path)
parser.add_argument('--output-dir', required=True, type=Path)
args = parser.parse_args()
NEW, ROOT, OUT = args.results_root, args.snapshot_root, args.output_dir
OUT.mkdir(parents=True, exist_ok=False)
for folder in ['snapshot', 'ot_comparison_snapshot']:
    p = ROOT / folder
    for name, expected in json.loads((p / 'TRANSFER_SHA256.json').read_text()).items():
        assert hashlib.sha256((p / name).read_bytes()).hexdigest() == expected, name

real = pd.read_csv(ROOT / 'ot_comparison_snapshot/domain_metrics_three_real.csv')
sim = pd.read_csv(ROOT / 'ot_comparison_snapshot/domain_metrics_all_simulations.csv')
rt = pd.read_csv(ROOT / 'ot_comparison_snapshot/transition_metrics_both_OT.csv')
rt = rt.loc[rt.family == 'current'].copy()
st = pd.read_csv(ROOT / 'snapshot/tables/summary_transition.csv')
st = st.loc[(st.kind == 'simulation') & (st.table == 'WOT_annotation')].copy()
domain = pd.concat([real.assign(kind='real'), sim], ignore_index=True)
trans = pd.concat([rt.assign(kind='real'), st], ignore_index=True)
metrics = ['ARI', 'PAS', 'CHAOS', 'naive_accuracy']
scales = [.01, .03, .1, .3, 1.]
labels = {'MouseEmbryo':'MouseOrganogenesis', 'ChickenHeart':'Chicken Cardiogenesis',
          'AxolotlTelencephalon':'Axolotl Telencephalon Neurogenesis', 'all_simulations':'Simulation Study'}
domain_seeds = {'MouseEmbryo':[666,1993,2026], 'ChickenHeart':[666,1993,2026], 'AxolotlTelencephalon':[666,1993,2026]}
ot_seeds = {'MouseEmbryo':[666,1993,2026], 'ChickenHeart':[666,1993,2026], 'AxolotlTelencephalon':[666,1993,2026]}
keys = ['kind', 'dataset', 'kernel_scale', 'seed']
dr = domain.groupby(keys)[['ARI','PAS','CHAOS']].mean().reset_index()
dr['n_timepoints'] = domain.groupby(keys).size().values
tr = trans.groupby(keys).naive_accuracy.mean().reset_index()
tr['n_pairs'] = trans.groupby(keys).size().values
runs = dr.merge(tr, on=keys, how='outer', validate='one_to_one')
runs['domain_in_figure'] = runs.apply(lambda r: r.kind == 'simulation' or r.seed in domain_seeds[r.dataset], axis=1)
runs['OT_in_figure'] = runs.apply(lambda r: pd.notna(r.naive_accuracy) and (r.kind == 'simulation' or r.seed in ot_seeds[r.dataset]), axis=1)
runs['training_status'] = 'completed_99_epochs'
assert len(runs) == 719 and runs.loc[runs.kind=='simulation'].shape[0] == 675

units = []
for ds in labels:
    for metric in metrics:
        q = runs.loc[runs.kind=='simulation'] if ds=='all_simulations' else runs.loc[runs.dataset==ds]
        q = q.loc[q.OT_in_figure if metric=='naive_accuracy' else q.domain_in_figure].copy()
        for scale in scales:
            s = q.loc[q.kernel_scale==scale]
            if ds=='all_simulations':
                for case, group in s.groupby('dataset'):
                    assert len(group)==3 and set(group.seed)=={666,1993,2026}
                    units.append([labels[ds], metric, scale, case, float(group[metric].mean()), '666,1993,2026'])
            else:
                for row in s.sort_values('seed').itertuples():
                    units.append([labels[ds], metric, scale, f'seed{row.seed}', float(getattr(row,metric)), str(row.seed)])
rep = pd.DataFrame(units, columns=['dataset','metric','kernel_coefficient_2l2','replicate','value','training_seeds'])
summary, groups = [], []
for ds in labels:
    for metric in metrics:
        for scale in scales:
            idx = rep.index[(rep.dataset==labels[ds]) & (rep.metric==metric) & (rep.kernel_coefficient_2l2==scale)].tolist()
            assert not idx or idx == list(range(min(idx),max(idx)+1))
            vals = rep.loc[idx,'value']
            n = len(vals)
            note = 'OT failed for both completed seeds; not imputed' if n==0 else ''
            seedtext = ','.join(sorted(set(','.join(rep.loc[idx,'training_seeds']).split(',')) - {''}, key=int)) if idx else ''
            summary.append([labels[ds],scale,math.sqrt(scale/2),metric,None if n==0 else float(vals.mean()),None if n<2 else float(vals.std(ddof=1)),n,'simulation setting' if ds=='all_simulations' else 'training seed',seedtext,note])
            groups.append(None if not idx else [min(idx)+2,max(idx)+2])
headers=['dataset','kernel_coefficient_2l2','conventional_length_l','metric','mean','sample_SD','n','summary_unit','training_seeds','missing_note']
outsum=pd.DataFrame(summary,columns=headers)

# Reproduce every plotted primary mean and SD, including the Mouse OT gap.
oldreal=pd.read_csv(NEW/'real_metric_summary_both_OT.csv')
oldsim=pd.read_csv(ROOT/'final_metric_summary.csv')
for ds,label in labels.items():
    old=oldsim.loc[(oldsim.dataset==ds)&oldsim.metric.isin(metrics)] if ds=='all_simulations' else oldreal.loc[(oldreal.dataset==ds)&(oldreal.family=='current')&oldreal.metric.isin(metrics)]
    new=outsum.loc[outsum.dataset==label].rename(columns={'kernel_coefficient_2l2':'kernel_scale'})
    c=new.merge(old,on=['kernel_scale','metric'],suffixes=('_new','_old'),validate='one_to_one')
    assert len(c)==20
    assert np.allclose(c.mean_new,c.mean_old,atol=1e-12,equal_nan=True)
    assert np.allclose(c.sample_SD,c.sd,atol=1e-12,equal_nan=True)

notes=[
 ['Scope','Exact descriptive statistics for the primary kernel-sensitivity figure. Association correlation is excluded.'],
 ['Metrics','ARI (higher), PAS (lower), CHAOS (lower), naive transition accuracy (higher).'],
 ['Kernel coefficient','The software kernel_scale equals 2l^2. Conventional RBF length l=sqrt(kernel_scale/2).'],
 ['Kernel covariance','K(i,j)=exp[-distance(s_i,s_j)^2/(2l^2)]. Coordinates are standardized independently within each timepoint.'],
 ['Fixed coefficient','The kernel coefficient is fixed during training; it is not learned.'],
 ['Real domain aggregation','Average timepoint scores within each training run, then calculate mean and sample SD across the included training seeds.'],
 ['Real transition aggregation','Average successive-timepoint-pair scores within each training run, then calculate mean and sample SD across included training seeds.'],
 ['Simulation aggregation','Average timepoints or pairs within each run, then average the three training seeds within each simulation setting. Calculate mean and sample SD over 45 settings.'],
 ['Weighting','Timepoints/pairs, training seeds, and simulation settings receive equal weight at their respective aggregation stages. Scores are not pooled over spots.'],
 ['Error bars','Sample SD (denominator n-1). Error bars are not SEM or confidence intervals.'],
 ['Real n','n is the number of training seeds included in a plotted point. These are computational repeats, not independent biological replicates.'],
 ['Simulation n','n=45 distinct simulation settings per kernel coefficient. Each averages three training seeds. Do not treat 135 training runs as 135 independent settings.'],
 ['Simulation design','3 DEG levels (100,200,500) x 3 data-generation seeds x 5 mixture proportions =45 settings. Training seeds are separate from data-generation seeds.'],
 ['Mouse domain seeds','All available seeds666,1993,2026 at each coefficient. n=3 except coefficient0.03, where seed1993 training failed (n=2).'],
 ['Chicken domain seeds','666,1993,2026 at all five coefficients.'],
 ['Axolotl domain seeds','666,1993,2026 at all five coefficients.'],
 ['Mouse OT seeds','All available seeds666,1993,2026. n=3 except coefficient0.03: seed1993 training failed and the two completed seeds failed post-training OT, so n=0.'],
 ['Chicken OT seeds','All available seeds666,1993,2026. n=3 except coefficient0.1 (n=2) because seed1993 failed post-training OT at that coefficient only.'],
 ['Axolotl OT seeds','666,1993,2026 at all five coefficients.'],
 ['Missing training','Mouse coefficient0.03/seed1993 failed at training epoch70 before producing a final embedding. Original protocol and failure were preserved.'],
 ['No imputation','Use all available prespecified seeds at each coefficient, excluding only the failed case. Missing values are not zero; no seed42 substitution.'],
 ['Naive accuracy definition','Sum of Q entries on reference correspondence edges divided by the sum of all Q entries. No edge threshold is applied.'],
 ['Q definition','Aggregate spot-level transport by original annotation, then take the entrywise minimum of row-normalized and column-normalized weights.'],
 ['Naive interpretation','Agreement with benchmark reference region correspondences. It is not the fraction of individual spots with validated lineage assignments.'],
 ['No curated-score cutoff','The naive score uses all entries. The separate thresholded curated accuracy is not used in this workbook.'],
 ['Analysis OT','Spot-level WOT: lambda1=0.1; lambda2=5; epsilon=0.02; epsilon0=1; growth_iters=3; tau=1000. The third transport map is used.'],
 ['Training OT','Training NumPy OT epsilon=0.05, lambda1=0.1, lambda2=5. Post-training OT is a separate calculation.'],
 ['Training protocol','Official SpaDOT1.0.5; 99 executed epochs (1-99); KL cycle30; strict Torch; one numerical thread. Inputs fixed across kernel scales.'],
 ['Readout','Final KMeans seed2026/n_init10, with dataset-specific fixed cluster counts.'],
 ['Tests','No hypothesis tests, p-values, multiplicity corrections, or significance labels were used to generate the sensitivity figure.'],
 ['Limits','SD across simulation settings describes heterogeneity between settings. Small real-data seed counts limit conclusions about seed variability.'],
 ['Source 1','Hash-verified sensitivity snapshot: domain_metrics_three_real.csv and domain_metrics_all_simulations.csv.'],
 ['Source 2','Hash-verified transition_metrics_both_OT.csv, current family only, and simulation WOT_annotation rows from summary_transition.csv.'],
 ['Verification','All 80 primary mean/SD points independently reproduced from included replication units to tolerance1e-12.'],
]

# Share flat tables without any correlation columns or curated accuracy values.
outsum.to_csv(OUT/'kernel_sensitivity_summary.csv',index=False)
rep.to_csv(OUT/'replication_units.csv',index=False)
runs.sort_values(keys).to_csv(OUT/'per_run_metrics.csv',index=False)
domain[['kind','dataset','kernel_scale','seed','timepoint','ARI','PAS','CHAOS']].to_csv(OUT/'per_timepoint_domain_metrics.csv',index=False)
trans[['kind','dataset','kernel_scale','seed','i','naive_accuracy']].to_csv(OUT/'per_pair_naive_accuracy.csv',index=False)
pd.DataFrame(notes,columns=['item','description']).to_csv(OUT/'statistical_methods.csv',index=False)

run_columns=['kind','dataset','kernel_scale','seed','ARI','PAS','CHAOS','naive_accuracy','n_timepoints','n_pairs','domain_in_figure','OT_in_figure','training_status']
def records(frame):
    return json.loads(frame.to_json(orient='values', double_precision=15))
payload=dict(summary_headers=headers,summary=records(outsum),groups=groups,
    replicate_headers=list(rep.columns),replicates=records(rep),
    run_headers=run_columns,runs=records(runs[run_columns].sort_values(keys)),
    methods_headers=['Item','Description'],methods=notes)
(OUT/'workbook_data.json').write_text(json.dumps(payload))
(OUT/'verification.json').write_text(json.dumps(dict(source_transfer_hashes_verified=True,primary_summary_points_reproduced=80,
    numerical_tolerance=1e-12,summary_rows=len(summary),replication_unit_rows=len(rep),per_run_rows=len(runs),
    per_timepoint_rows=len(domain),per_pair_rows=len(trans),association_r_excluded=True),indent=2)+'\n')
print(outsum.loc[outsum.kernel_coefficient_2l2==.1,['dataset','metric','mean','sample_SD','n']].to_string(index=False))
print('OUTPUT',OUT)
