"""Aggregate raw simulation ECDFs for Fig. S46.

Usage: python aggregate_s46.py --input distance_CDFs.json --output aggregated_curves.json
Proportions and seeds receive equal weights. Intervals are pointwise Student-t.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='1'
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.stats import t
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
source = args.input
records=json.loads(source.read_text())
assert len(records)==630
METHODS=['PRECAST','SEDR','spaGCN','GraphST','spatialPCA','spatiAlign','spaVAE']
COLORS=dict(zip(METHODS,['#68a756','#f09232','#4c2623','#165451','#d3b300','#ded6c1','#0082c0']))
DEGS=[100,200,500];SEEDS=[1105,1993,2024];PROPS=[.5,.6,.7,.8,.9];TRS=['t0_t1','t1_t2']
index={(r['deg'],r['method'],r['transition'],r['seed'],r['proportion']):r for r in records}
assert len(index)==630
# Fixed display range retained for all methods except spatialPCA. Its full positive range
# is shared across all three DEG levels and both transitions, without omitting data.
pca_curves=[r[mode] for r in records if r['method']=='spatialPCA' for mode in ['multi','single']]
pca_pos=[x for c in pca_curves for x in c['distance'] if x>0]
pca_limits=(10.**np.floor(np.log10(min(pca_pos))),10.**np.ceil(np.log10(max(pca_pos))))
limits={m:(pca_limits if m=='spatialPCA' else (1.,100.)) for m in METHODS}
critical=float(t.ppf(.975,2));assert abs(critical-4.302652729696142)<1e-9
summaries=[];lookup={};max_difference=0.
for deg in DEGS:
 for method in METHODS:
  # Include quantile knots for fine rendering on a shared x grid per method.
  lo,hi=limits[method];grid=np.geomspace(lo,hi,601)
  for tr in TRS:
   for mode in ['multi','single']:
    seed_curves=[];all_curves=[];pair_counts=[]
    for seed in SEEDS:
     curves=[]
     for prop in PROPS:
      rec=index[(deg,method,tr,seed,prop)];c=rec[mode]
      assert c['n_pairs']==rec['matched_source']*rec['matched_target']
      xx=np.asarray(c['distance']);yy=np.asarray(c['cdf'])
      assert np.all(np.diff(xx)>=0) and np.all(np.diff(yy)>0)
      curve=np.interp(grid,xx,yy,left=0.,right=1.)
      assert np.isfinite(curve).all() and np.all(np.diff(curve)>=-1e-12)
      curves.append(curve);all_curves.append(curve);pair_counts.append(c['n_pairs'])
     assert len(curves)==5
     seed_curves.append(np.mean(curves,axis=0))
    y=np.stack(seed_curves);mean=y.mean(axis=0);sd=y.std(axis=0,ddof=1)
    half=critical*sd/np.sqrt(3);lower=np.maximum(0,mean-half);upper=np.minimum(1,mean+half)
    direct=np.sum(np.stack(all_curves),axis=0)/15
    diff=float(np.max(abs(direct-mean)));max_difference=max(max_difference,diff);assert diff<1e-12
    # Independent variance calculation and confidence-band coverage of the estimated mean.
    manual_sd=np.sqrt(((y[0]-mean)**2+(y[1]-mean)**2+(y[2]-mean)**2)/2)
    assert np.allclose(sd,manual_sd,rtol=1e-12,atol=1e-14)
    assert np.all(lower<=mean) and np.all(mean<=upper) and np.all(np.diff(mean)>=-1e-12)
    result={'deg':deg,'method':method,'transition':tr,'mode':mode,'distance':grid.tolist(),'mean_cdf':mean.tolist(),'lower_95':lower.tolist(),'upper_95':upper.tolist(),'seed_curves':y.tolist(),'seeds':SEEDS,'proportions':PROPS,'n_independent_seeds':3,'n_conditions_averaged':15,'pair_counts_per_condition':pair_counts}
    summaries.append(result);lookup[(deg,method,tr,mode)]=result
assert len(summaries)==84
args.output.parent.mkdir(parents=True, exist_ok=True)
with args.output.open('x') as stream:
    stream.write(json.dumps(summaries, separators=(',', ':')) + '\n')
