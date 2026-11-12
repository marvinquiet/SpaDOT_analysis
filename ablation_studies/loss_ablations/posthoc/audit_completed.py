"""Independent CPU metric/objective audit of completed cases; never trains."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
import sys
import json
import argparse
from pathlib import Path
ROOT=Path(os.environ['SPADOT_LOSS_STUDY_ROOT']).resolve()
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'training'));sys.path.insert(0,str(ROOT/'source'))
from common import sha,check_frozen,write_json
import torch
import numpy as np
import pandas as pd
import anndata as ad
import yaml
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits
from SpaDOT.utils._train_utils import _beta_cycle_linear_period


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--destination',required=True)
    parser.add_argument('--trial',action='append');args=parser.parse_args()
    destination=Path(args.destination);destination.mkdir(exist_ok=False)
    manifest=json.loads((ROOT/'manifest.json').read_text());check_frozen(manifest)
    status=json.loads((ROOT/'status.json').read_text())
    ids=args.trial or status.get('completed_ids',[])
    assert ids and len(ids)==len(set(ids))
    graphs=np.load(ROOT/'evaluation/fixed_spatial_neighbors.npz');audits=[]
    for trial in ids:
        job=next(j for j in manifest['jobs'] if j['id']==trial)
        directory=ROOT/'runs'/trial;result=json.loads((directory/'result.json').read_text())
        assert result['all_passed'] and result['executed_epoch']==99 and result['retries']==0
        assert all(sha(directory/p)==h for p,h in result['output_hashes'].items())
        predictions=pd.read_csv(directory/'predictions.csv',index_col=0,dtype=str)
        config=yaml.safe_load(Path(job['config']).read_text());info=manifest['datasets'][job['dataset']]
        metrics={};fit_count=0
        for child in sorted(directory.iterdir()):
            if not child.is_dir() or not (child/'result.json').exists():continue
            fit=json.loads((child/'result.json').read_text());fit_count+=1
            assert all(sha(child/p)==h for p,h in fit['output_hashes'].items())
            loss=pd.read_csv(child/'loss.csv',index_col=0)
            assert loss.index.tolist()==list(range(1,100)) and np.isfinite(loss.values).all()
            beta=_beta_cycle_linear_period(99,30,stop=config['beta1'])
            replay=(config['lambda1']*loss.Recon-beta*loss.SVGP_KL+config['beta2']*loss.GAT_KL
                +config['omiga1']*loss.alignment+config['omiga2']*loss.KMeans+config['omiga3']*loss.OT)
            assert np.allclose(replay,loss.elbo,rtol=1e-10,atol=1e-8)
            refresh=list(range(10,100,10))
            if config['omiga3']==0:
                assert fit['inactive_ot_guard'] and fit['ot_refresh_epochs']==[]
                assert fit['inactive_ot_skipped_epochs']==refresh and (loss.OT==0).all()
            else:
                assert not fit['inactive_ot_guard'] and fit['inactive_ot_skipped_epochs']==[]
                assert fit['ot_refresh_epochs']==refresh
            state=torch.load(child/'SpaDOT_model.pth',map_location='cpu',weights_only=True)
            assert all(torch.isfinite(x).all() for x in state.values() if x.is_floating_point())
            latent=ad.read_h5ad(child/'latent.h5ad')
            assert latent.n_vars==20 and np.isfinite(latent.X).all() and latent.obs_names.is_unique
            assert fit['original_class_checkpoint_exact'] and fit['checkpoint_max_abs_error']==0
            for tp in fit['timepoints']:
                part=latent[latent.obs.timepoint.astype(str)==tp]
                k=info['final_clusters'][info['timepoints'].index(tp)]
                with threadpool_limits(limits=1):
                    labels=KMeans(n_clusters=k,n_init=10,random_state=2026,init='k-means++',algorithm='lloyd').fit_predict(part.X)
                assert adjusted_rand_score(predictions.loc[part.obs_names,'kmeans'],labels)==1
                neighbors=graphs[job['dataset']+'__'+tp]
                score=dict(ARI=float(adjusted_rand_score(part.obs.annotation.astype(str),labels)),
                           PAS=float(np.mean(np.sum(labels[neighbors]!=labels[:,None],axis=1)>=6)))
                assert all(score[m]==result['scores'][tp][m] for m in ['ARI','PAS'])
                metrics[tp]=score
        assert fit_count==(1 if job['mode']=='combined' else len(info['timepoints']))
        assert set(metrics)==set(info['timepoints'])
        means={m:float(np.mean([x[m] for x in metrics.values()])) for m in ['ARI','PAS']}
        assert all(means[m]==result['mean_'+m] for m in ['ARI','PAS'])
        audits.append(dict(trial=trial,all_passed=True,dataset=job['dataset'],mode=job['mode'],
            method=job['method']['id'],mean_ARI=means['ARI'],mean_PAS=means['PAS'],scores=metrics,
            independently_verified_fits=fit_count,executed_epoch=99))
    check_frozen(manifest)
    write_json(destination/'audit.json',dict(all_passed=True,cases=audits,training_steps=0,
        cuda_access=False,original_sources_inputs_preserved=True))
    print(json.dumps(audits),flush=True)


if __name__=='__main__':main()
