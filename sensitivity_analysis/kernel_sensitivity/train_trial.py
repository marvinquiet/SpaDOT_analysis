"""Run one S47 kernel-scale trial from the external frozen study manifest.

Set SPADOT_KERNEL_STUDY_ROOT to the workspace holding the original frozen
SpaDOT 1.0.5 source, kernel-scale configs and data manifest. It must point to
preprocessed real/simulation H5ADs available locally. Provide --trial, --gpu,
--output-dir and PYTHONHASHSEED matching the training seed before invoking Python.
The training grid uses 0.01, 0.03, 0.1, 0.3 and 1 with seeds 666, 1993 and 2026.
The default post-training OT from this original script is evaluated separately
from the paper-parameter sensitivity curve. See score_paper_ot.py.
"""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
os.environ['CUBLAS_WORKSPACE_CONFIG']=':4096:8'
from pathlib import Path
import sys,json,hashlib,time,traceback,argparse,importlib.metadata
from types import SimpleNamespace
ROOT=Path(os.environ['SPADOT_KERNEL_STUDY_ROOT']).resolve()
sys.path.insert(0,str(ROOT/'source'))
import numpy as np,pandas as pd,anndata as ad,torch,yaml
from scipy import sparse
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
from threadpoolctl import threadpool_limits
import SpaDOT
from SpaDOT.utils import _analyze_utils as A,_train_utils as U
import collaborator_metrics as M

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def write(p,x):
 q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(x,indent=2)+'\n');q.replace(p)
def dense(x):return x.toarray() if sparse.issparse(x) else np.asarray(x)
def raw_table(T,src,dst,rows,cols):
 return pd.DataFrame((src[:,None]==np.array(rows)[None,:]).astype(float).T@(T@(dst[:,None]==np.array(cols)[None,:]).astype(float)),index=rows,columns=cols)

def main():
 p=argparse.ArgumentParser();p.add_argument('--trial',required=True);p.add_argument('--gpu',type=int,required=True);p.add_argument('--output-dir',type=Path,required=True);args=p.parse_args()
 manifest=json.loads((ROOT/'manifest.json').read_text());trial=next(t for t in manifest['trials'] if t['id']==args.trial);case=next(c for c in manifest['cases'] if c['id']==trial['dataset'])
 out=args.output_dir/trial['dataset']/f"kernel_{trial['scale']:g}_seed{trial['seed']}";out.mkdir(parents=True,exist_ok=False)
 status={'id':trial['id'],'dataset':trial['dataset'],'scale':trial['scale'],'seed':trial['seed'],'gpu':args.gpu,'status':'running','stage':'preflight','started_unix':time.time()};write(out/'status.json',status)
 def stage(name):status['stage']=name;write(out/'status.json',status);print('STAGE',name,flush=True)
 try:
  assert Path(SpaDOT.__file__).resolve()==ROOT/'source/SpaDOT/__init__.py'
  assert importlib.metadata.version('SpaDOT')=='1.0.5'
  for rel,h in json.loads((ROOT/'source_sha256.json').read_text()).items():assert sha(ROOT/'source'/rel)==h,rel
  for name,h in manifest['code_sha256'].items():assert sha(ROOT/'code'/name)==h,name
  assert sha(Path(case['input']))==case['input_sha256']
  cfg_path=Path(trial['config']);assert sha(cfg_path)==manifest['config_sha256'][cfg_path.name]
  cfg=yaml.safe_load(cfg_path.read_text());reference=yaml.safe_load((ROOT/'source/SpaDOT/config.yaml').read_text())
  assert {k for k in cfg if cfg[k]!=reference.get(k)} <= {'kernel_scale','seed'}
  assert cfg['kernel_scale']==trial['scale'] and cfg['seed']==trial['seed']
  assert os.environ.get('PYTHONHASHSEED')==str(trial['seed'])
  torch.use_deterministic_algorithms(True);torch.set_num_threads(1);torch.set_num_interop_threads(1)
  torch.backends.cudnn.deterministic=True;torch.backends.cudnn.benchmark=False;torch.cuda.set_device(args.gpu)
  x=torch.eye(2,device=f'cuda:{args.gpu}',dtype=torch.float64);assert torch.equal(x@x,x)
  (out/'config_used.yaml').write_text(cfg_path.read_text())
  observer={}
  original_train=U.train_SpaDOT
  def observe(*a,**kw):
   model,loss=original_train(*a,**kw)
   observer['all_model_parameters_finite']=all(bool(torch.isfinite(v).all()) for v in model.parameters())
   observer['kernel_scales_after_training']={t:float(g.kernel.scale.item()) for t,g in model.svgp_dict.items()}
   observer['kernel_has_trainable_parameters']=any(len(list(g.kernel.parameters()))>0 for g in model.svgp_dict.values())
   assert observer['all_model_parameters_finite'] and not observer['kernel_has_trainable_parameters']
   assert all(v==trial['scale'] for v in observer['kernel_scales_after_training'].values())
   return model,loss
  U.train_SpaDOT=observe
  stage('training');start=time.monotonic()
  with threadpool_limits(limits=1):SpaDOT.train(SimpleNamespace(data=case['input'],output_dir=str(out),prefix='',config=str(cfg_path),save_model=True,device=f'cuda:{args.gpu}'))
  seconds=time.monotonic()-start;loss=pd.read_csv(out/'loss.csv',index_col=0)
  assert [int(i) for i in loss.index]==list(range(1,100)) and np.isfinite(loss.values).all()
  stage('analysis_and_metrics')
  lat=ad.read_h5ad(out/'latent.h5ad');assert lat.shape==(case['shape'][0],20) and np.isfinite(lat.X).all()
  input_obs=ad.read_h5ad(case['input'],backed='r');assert set(lat.obs_names)==set(input_obs.obs_names)
  assert np.array_equal(lat.obsm['spatial'],input_obs[lat.obs_names].obsm['spatial']);input_obs.file.close()
  lat.obs['timepoint']=lat.obs.timepoint.astype(str);lat.obs['annotation']=lat.obs.annotation.astype(str).map(M.norm)
  with threadpool_limits(limits=1):
   pred=A.KMeans_Clustering(lat,case['clusters'],seed=2026,n_init=10)
   A.OT_analysis(SimpleNamespace(output_dir=str(out),prefix=''),pred)
   A.plot_OT(SimpleNamespace(output_dir=str(out),prefix=''),pred)
   A.plot_domains(SimpleNamespace(output_dir=str(out),prefix=''),pred)
  pred.obs['kmeans']=pred.obs.kmeans.astype(str)
  pd.DataFrame({'spot':pred.obs_names,'timepoint':pred.obs.timepoint.values,'cluster':pred.obs.kmeans.values,'annotation':pred.obs.annotation.values}).to_csv(out/'predictions.csv',index=False)
  graphs=M.transition_graphs(case['id'],case['graph'],sorted(pred.obs.annotation.unique()))
  domain=[];transition=[];audit=[]
  for tp in case['timepoints']:
   mask=pred.obs.timepoint==tp;o=pred.obs.loc[mask];loc=pred.obsm['spatial'][mask]
   assert o.kmeans.nunique()==case['clusters'][case['timepoints'].index(tp)]
   domain.append(dict(dataset=case['id'],kind=case['kind'],kernel_scale=trial['scale'],seed=trial['seed'],timepoint=tp,n=len(o),ARI=adjusted_rand_score(o.annotation,o.kmeans),NMI=normalized_mutual_info_score(o.annotation,o.kmeans),PAS=M.compute_PAS(loc,o.kmeans.values),CHAOS=M.compute_CHAOS(loc,o.kmeans.values),PAS_annotation=M.compute_PAS(loc,o.annotation.values),CHAOS_annotation=M.compute_CHAOS(loc,o.annotation.values)))
  for i in range(len(case['timepoints'])-1):
   raw=ad.read_h5ad(out/f'OT_{i}_{i+1}.h5ad');T=dense(raw.X);assert np.isfinite(T).all() and (T>=0).all()
   src=pred.obs.loc[raw.obs_names];dst=pred.obs.loc[raw.var_names];s=src.annotation.to_numpy();d=dst.annotation.to_numpy()
   annotation_mass=raw_table(T,s,d,sorted(set(s)),sorted(set(d)))
   ctab=ad.read_h5ad(out/f'transition_table_{i}_{i+1}.h5ad')
   cluster_mass=raw_table(T,src.SpaDOT_pred_labels.to_numpy(),dst.SpaDOT_pred_labels.to_numpy(),list(ctab.obs_names),list(ctab.var_names))
   diff=float(np.abs(cluster_mass.values/cluster_mass.values.sum()-ctab.X/ctab.X.sum()).max());assert diff<1e-10
   majority=[]
   for tp in case['timepoints'][i:i+2]:majority.append(pred.obs.loc[pred.obs.timepoint==tp].groupby('kmeans')['annotation'].agg(lambda x:x.value_counts().idxmax()).to_dict())
   mapped=cluster_mass.copy();mapped.index=[majority[0][x.split('_')[-1]] for x in mapped.index];mapped.columns=[majority[1][x.split('_')[-1]] for x in mapped.columns]
   mapped=mapped.groupby(level=0).sum().T.groupby(level=0).sum().T
   for kind,mass in [('WOT_annotation',annotation_mass),('SpaDOT_domains_majority',mapped),('SpaDOT_domains_raw',cluster_mass)]:
    Q=M.min_norm(mass);assert np.isfinite(Q.values).all()
    mass.to_csv(out/f'mass_{kind}_{i}_{i+1}.csv');Q.to_csv(out/f'transition_{kind}_{i}_{i+1}.csv')
    if kind!='SpaDOT_domains_raw':
     acc,weight=M.evaluate_using_curated_transitions(Q,graphs[i])
     transition.append(dict(dataset=case['id'],kind=case['kind'],kernel_scale=trial['scale'],seed=trial['seed'],i=i,transition=f't{i}->t{i+1}',table=kind,naive_accuracy=M.accuracy(Q,graphs[i]),curated_accuracy=acc,total_weight=weight))
   audit.append(dict(pair=i,cluster_aggregation_max_error=diff,spot_shape=list(T.shape)))
  pd.DataFrame(domain).to_csv(out/'domain_metrics.csv',index=False);pd.DataFrame(transition).to_csv(out/'transition_metrics.csv',index=False)
  baseline={}
  if case['kind']=='real' and trial['scale']==.1 and trial['seed']==2026:
   original=ad.read_h5ad(case['baseline_latent']);by_index=lat[original.obs_names]
   baseline['latent_exact']=bool(np.array_equal(by_index.X,original.X));baseline['latent_max_abs_error']=float(np.abs(by_index.X-original.X).max())
   prior_loss=pd.read_csv(case['baseline_loss'],index_col=0);baseline['loss_exact']=bool(np.array_equal(prior_loss.values,loss.values));baseline['loss_max_abs_error']=float(np.abs(prior_loss.values-loss.values).max())
  assert sha(Path(case['input']))==case['input_sha256']
  report=dict(all_passed=True,trial=trial,case=case['id'],executed_epoch=99,training_seconds=seconds,model=observer,spot_aggregation=audit,baseline_reproduction=baseline,source_manifest_sha256=sha(ROOT/'source_sha256.json'),input_sha256=case['input_sha256'],analysis_ot=manifest['analysis_ot'])
  write(out/'verification.json',report)
  files={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name not in ['status.json','complete.json']}
  write(out/'complete.json',dict(all_passed=True,files=files,trial_id=trial['id']))
  status.update(status='complete',stage='complete',mean_ARI=float(np.mean([x['ARI'] for x in domain])),ARI={x['timepoint']:x['ARI'] for x in domain},training_seconds=seconds,finished_unix=time.time());write(out/'status.json',status)
  print('COMPLETE',json.dumps(status),flush=True)
 except Exception as e:
  status.update(status='failed',error=str(e),traceback=traceback.format_exc(),finished_unix=time.time());write(out/'status.json',status);print(traceback.format_exc(),flush=True);raise
if __name__=='__main__':main()
