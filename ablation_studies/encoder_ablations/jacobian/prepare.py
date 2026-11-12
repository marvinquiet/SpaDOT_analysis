"""Prepare S40b from the archived seed-2026 SpaDOT checkpoint.

Set SPADOT_ENCODER_STUDY_ROOT to the original parent study with manifest,
configs, checkpoint and frozen source. Set SPADOT_JACOBIAN_WORK_ROOT to a
new external work directory. Input and latent H5ADs are required separately.
"""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[name]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
import hashlib,json,shutil,time
R=Path(os.environ['SPADOT_JACOBIAN_WORK_ROOT']).resolve()
R.mkdir(parents=True, exist_ok=True)
P=Path(os.environ['SPADOT_ENCODER_STUDY_ROOT']).resolve()
P_SHA='1d6e5de97b593273ae6b70eed8a283b01dee55eb4537677ae2c9607fa84e6adc'
ID='AxolotlTelencephalon__dual_native__seed2026'

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');q.replace(p)
def farthest(coordinates,n=20):
 import numpy as np
 x=np.asarray(coordinates,dtype=float);x=(x-x.mean(0))/x.std(0)
 first=int(np.argmin(np.sum((x-x.mean(0))**2,axis=1)))
 chosen=[first];minimum=np.sum((x-x[first])**2,axis=1)
 for _ in range(1,min(n,len(x))):
  minimum[chosen]=-1;new=int(np.argmax(minimum));chosen.append(new)
  minimum=np.minimum(minimum,np.sum((x-x[new])**2,axis=1))
 return chosen

def main():
 import numpy as np,anndata as ad,importlib.metadata
 assert not(R/'manifest.json').exists() and not(R/'source').exists()
 assert sha(P/'manifest.json')==P_SHA
 parent=read(P/'manifest.json');trial=next(t for t in parent['trials'] if t['id']==ID)
 case=next(c for c in parent['cases'] if c['id']==trial['dataset']);O=P/'results'/ID
 v=read(O/'verification.json');done=read(O/'complete.json')
 assert done['training_passed'] and done['domain_passed'] and v['checkpoint_reconstruction_exact']
 consumed={str(P/'manifest.json'):P_SHA}
 for name in ['complete.json','verification.json','initialization.json','native_sampler.json','checkpoint_reconstruction.json']:
  consumed[str(O/name)]=sha(O/name)
 for name in ['SpaDOT_model.pth','latent.h5ad','config_used.yaml']:
  expected=done['files'][name];assert sha(O/name)==expected;consumed[str(O/name)]=expected
 consumed[case['input']]=case['input_sha256'];assert sha(case['input'])==case['input_sha256']
 for name,h in parent['source_sha256'].items():assert sha(P/'source'/name)==h;consumed[str(P/'source'/name)]=h
 shutil.copytree(P/'source',R/'source',ignore=shutil.ignore_patterns('__pycache__','*.pyc','._*'))
 for name,h in parent['source_sha256'].items():assert sha(R/'source'/name)==h
 packages={n:importlib.metadata.version(n) for n in parent['package_versions']};assert packages==parent['package_versions']
 data=ad.read_h5ad(case['input'],backed='r');latent=ad.read_h5ad(O/'latent.h5ad')
 assert tuple(case['shape'])==data.shape and latent.shape==(data.n_obs,20)
 assert data.obs_names.is_unique and latent.obs_names.is_unique and set(data.obs_names)==set(latent.obs_names)
 assert np.isfinite(latent.X).all()
 anchors={}
 for tp in case['timepoints']:
  mask=np.asarray(data.obs.timepoint.astype(str)==tp);spots=np.asarray(data.obs_names[mask],dtype=str)
  coords=np.asarray(data.obsm['spatial'][mask],dtype=float);ids=farthest(coords,20)
  anchors[tp]=dict(indices=ids,spots=spots[ids].tolist(),n_spots=len(spots),selection_uses_labels=False,
   rule='First closest standardized-coordinate centroid; then deterministic farthest point, ties source input order')
 data.file.close()
 reference=dict(input=case['input'],input_sha256=case['input_sha256'],checkpoint=str(O/'SpaDOT_model.pth'),checkpoint_sha256=consumed[str(O/'SpaDOT_model.pth')],
  latent=str(O/'latent.h5ad'),latent_sha256=consumed[str(O/'latent.h5ad')],config=str(O/'config_used.yaml'),config_sha256=consumed[str(O/'config_used.yaml')],
  initialization=str(O/'initialization.json'),native_sampler=str(O/'native_sampler.json'))
 manifest=dict(root=str(R),parent_root=str(P),parent_manifest_sha256=P_SHA,reference=reference,reference_result_dir=str(O),case=case,trial=trial,
  source_sha256=parent['source_sha256'],consumed_sha256=consumed,code_sha256={p.name:sha(p) for p in Path(__file__).resolve().parent.glob('*.py')},package_versions=packages,anchors=anchors,
  protocol=dict(backend='cpu',dtype='float64',numerical_threads=1,nice=10,training_seed=2026,no_training=True,no_OT=True,optimizer_steps=0,
   cross_device_atol=1e-9,cross_device_rtol=1e-10,anchors_per_timepoint=20,primary_influence='raw squared Frobenius energy normalized peranchor',
   secondary_influence='Fixed baseline covariance-ridge whitened squared energy',relative_ridge=1e-6,exact_basis_per_metric=10,
   native_graph_k={'t0':6,'t1':18,'t2':24},fixed_weights=True,coordinate_transforms=['identity','scale_0.5','scale_2','rigid_reflection'],
   coordinate_transform_inducing_locations_jointly=True,minimum_available_RAM_bytes=32*2**30,no_retry_or_fallback=True,
   source_context_GP='ALL sources within timepoint',GAT_context='Exact incoming3hop full-neighbor induced closure, verifiedagainst fullgraphforward',
   limitations=['Fixed-weight sensitivity is not a retraining influence function.','Three selectedtimepoints and one trainingseed are development diagnostics.','Geometrycounterfactual does not prove predictive usefulness.']))
 write(R/'manifest.json',manifest);write(R/'anchors.json',anchors)
 write(R/'preparation_verification.json',dict(all_passed=True,manifest_sha256=sha(R/'manifest.json'),consumed_files=len(consumed),anchors_per_timepoint=20,
  labels_not_used_for_selection=True,no_model_or_optimizer_executed=True,finished_unix=time.time()))
 print('PREPARED',sha(R/'manifest.json'),flush=True)
if __name__=='__main__':main()
