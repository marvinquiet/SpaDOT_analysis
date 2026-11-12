"""Compute S40b with exact CPU gradients of the archived Axolotl checkpoint.

Run prepare.py and mechanism_fixture.py first, passing the work root's frozen
source to --source and writing mechanism_fixture.json in that work root.
Use the same SPADOT_JACOBIAN_WORK_ROOT and SPADOT_ENCODER_STUDY_ROOT variables.
Requires Linux, SpaDOT, Torch, PyG, scikit-learn, anndata, NumPy, pandas and PyYAML.
"""
import os
os.environ['CUDA_VISIBLE_DEVICES']=''
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import json,sys,time,traceback,gc,hashlib
R=Path(os.environ['SPADOT_JACOBIAN_WORK_ROOT']).resolve()
sys.path[:0]=[str(R/'source'),str(R)]
from prepare import read,write,sha

def available_RAM():
 for line in Path('/proc/meminfo').read_text().splitlines():
  if line.startswith('MemAvailable:'):return int(line.split()[1])*1024
 raise RuntimeError('Cannot observe availableRAM')
def main():
 import numpy as np,pandas as pd,anndata as ad,torch,yaml
 from sklearn.preprocessing import StandardScaler
 from threadpoolctl import threadpool_limits
 from torch_geometric.utils import dense_to_sparse
 from SpaDOT.model.SpaDOT import SpaDOT
 from SpaDOT.utils import _train_utils as U,_utils
 from mechanism_core import array_digest,fit_whitener,exact_influence,summarize_influence,all_incoming_hops,coordinate_api
 m=read(R/'manifest.json');out=R/'analysis';out.mkdir(exist_ok=False);(out/'influence').mkdir();(out/'baseline').mkdir()
 status=dict(stage='initialization',started_unix=time.time(),pid=os.getpid(),backend='cpu',training=False,optimizer_steps=0,OT_solves=0,anchors_completed=0,anchors_planned=60)
 write(out/'status.json',status)
 def stage(name,**kwargs):status.update(stage=name,updated_unix=time.time(),**kwargs);write(out/'status.json',status);print(name,json.dumps(kwargs),flush=True)
 def verify():
  for p,h in m['consumed_sha256'].items():assert sha(p)==h,p
  for n,h in m['source_sha256'].items():assert sha(R/'source'/n)==h,n
  for n,h in m['code_sha256'].items():assert sha(Path(__file__).resolve().parent/n)==h,n
 try:
  assert os.environ.get('CUDA_VISIBLE_DEVICES')==''
  assert available_RAM()>=m['protocol']['minimum_available_RAM_bytes']
  assert read(R/'mechanism_fixture.json')['all_passed']
  verify();torch.set_num_threads(1);torch.set_num_interop_threads(1);torch.use_deterministic_algorithms(True)
  case=m['case'];data=ad.read_h5ad(m['reference']['input']);saved=ad.read_h5ad(m['reference']['latent'])
  cfg=yaml.safe_load(Path(m['reference']['config']).read_text());cfg.update(device=torch.device('cpu'),dtype=torch.float64,input_dim=data.n_vars,timepoints=case['timepoints'])
  assert cfg['seed']==2026 and cfg['z_dim']==20 and cfg['kernel_scale']==.1
  _utils.set_seed(2026);stage('native_dataloader_reconstruction')
  dl=U.prepare_dataloader(data.copy(),cfg)
  model=SpaDOT(cfg,dl).to('cpu');state=torch.load(m['reference']['checkpoint'],map_location='cpu',weights_only=True);model.load_state_dict(state,strict=True);del state
  model.eval()
  for parameter in model.parameters():parameter.requires_grad_(False)
  def state_digest():
   h=hashlib.sha256()
   for n,v in model.state_dict().items():h.update(n.encode());h.update(array_digest(v).encode())
   return h.hexdigest()
  before=state_digest()
  def forbidden(*a,**kw):raise AssertionError('Training/OT forbidden in read-only mechanism analysis')
  U.train_SpaDOT=forbidden;U._update_OT_matrix=forbidden
  stage('full_native_CPU_inference')
  with threadpool_limits(limits=1),torch.no_grad():baseline=U.get_latent(model,cfg,data,dl)
  native=saved[baseline.obs_names].X;new=np.asarray(baseline.X);difference=np.abs(new-native)
  tol=m['protocol'];assert np.allclose(new,native,atol=tol['cross_device_atol'],rtol=tol['cross_device_rtol']),('CPU/GPU checkpointmeans differ',float(difference.max()))
  reproduction=dict(all_passed=True,backend='cpu',reference_backend='savedGPU',bitwise_equal=bool(np.array_equal(new,native)),maximum_absolute_error=float(difference.max()),atol=tol['cross_device_atol'],rtol=tol['cross_device_rtol'],checkpoint_sha256=m['reference']['checkpoint_sha256'],model_weights_unchanged=True)
  write(out/'checkpoint_reproduction.json',reproduction)
  source_graph=read(m['reference']['native_sampler']);contexts={};whitening={};initial_U={}
  for tp in case['timepoints']:
   loc,y,indices=dl['datasets'][tp].dataset.tensors;edge,_=dense_to_sparse(dl['adjacency_matrices'][tp]);indices=indices.numpy();spots=np.asarray(data.obs_names[indices],dtype=str);raw=np.asarray(data.obsm['spatial'][indices],dtype=float)
   assert array_digest(edge)==source_graph[tp]['graph_sha256']
   graph_k=min(cfg['max_neighbors'],cfg['knn_cutoff']*round(len(spots)/1000));assert graph_k==m['protocol']['native_graph_k'][tp]
   z=np.asarray(baseline[spots].X);whiteners={};whitening[tp]={}
   for b,columns in [('GP',slice(0,10)),('GAT',slice(10,20))]:whiteners[b],whitening[tp][b]=fit_whitener(z[:,columns],m['protocol']['relative_ridge'])
   source_anchors=m['anchors'][tp];assert list(spots[source_anchors['indices']])==source_anchors['spots']
   scaler=StandardScaler().fit(raw);assert np.array_equal(scaler.transform(raw),loc.numpy())
   coords=loc.numpy();distance_edges=np.linspace(0,np.linalg.norm(coords.max(0)-coords.min(0))*(1+1e-12),16)
   contexts[tp]=dict(X=loc,Y=y,edge=edge,raw=raw,spots=spots,z=z,whiteners=whiteners,scaler=scaler,distance_edges=distance_edges,indices=indices)
   initial_U[tp]=array_digest(model.svgp_dict[tp].inducing_index_points)
   np.savez_compressed(out/'baseline'/(tp+'.npz'),spot=spots,latent=z,coordinates_raw=raw,coordinates_standardized=coords,edge_index=edge.numpy(),inducing_points=model.svgp_dict[tp].inducing_index_points.numpy())
  write(out/'whitening.json',whitening)
  index=[];summary=[];hops_rows=[];dist_rows=[];timings=[]
  # First preselectedcentroid anchor at EACH timepoint provides an honest pilot;
  # later anchors are processed in the same fixed round-robin order.
  for ordinal in range(20):
   for tp in case['timepoints']:
    context=contexts[tp];anchor=m['anchors'][tp]['indices'][ordinal];start=time.monotonic()
    stage('anchor_sensitivity',timepoint=tp,anchor_ordinal=ordinal,anchor_index=anchor,anchor_spot=str(context['spots'][anchor]))
    with threadpool_limits(limits=1):result=exact_influence(model,tp,context['X'],context['Y'],context['edge'],anchor,context['whiteners'],baseline_anchor=context['z'][anchor],atol=tol['cross_device_atol'],rtol=tol['cross_device_rtol'])
    hops=all_incoming_hops(context['edge'],len(context['spots']),anchor)
    filename='influence/'+tp+'_anchor%02d.npz'%ordinal
    payload=dict(hops=hops,coordinates_raw=context['raw'],coordinates_standardized=context['X'].numpy(),spot=context['spots'],anchor_index=np.array(anchor))
    for b in ['GP','GAT']:
     for metric,key in [('raw','raw_energy'),('whitened','whitened_energy')]:
      energy=result[b][key];payload[b+'_'+key]=energy
      stats=summarize_influence(energy,context['X'].numpy(),hops,anchor,context['distance_edges'])
      base=dict(timepoint=tp,anchor=anchor,anchor_ordinal=ordinal,anchor_spot=str(context['spots'][anchor]),branch=b,metric=metric)
      summary.append(dict(base,**{k:v for k,v in stats.items() if k not in ['hop_rows','distance_rows','normalized_mass','normalized_nonself_mass']}))
      hops_rows.extend(dict(base,**row) for row in stats['hop_rows']);dist_rows.extend(dict(base,**row) for row in stats['distance_rows'])
    np.savez_compressed(out/filename,**payload)
    metadata=dict(result['metadata'],timepoint=tp,anchor_ordinal=ordinal,elapsed_seconds=time.monotonic()-start)
    write(out/'influence'/(tp+'_anchor%02d.json'%ordinal),metadata)
    index.append(dict(timepoint=tp,anchor=anchor,anchor_ordinal=ordinal,anchor_spot=str(context['spots'][anchor]),npz=filename))
    timings.append(dict(timepoint=tp,anchor_ordinal=ordinal,seconds=metadata['elapsed_seconds']))
    write(out/'influence_index.json',index)
    for name,rows in [('influence_summary',summary),('hop_profiles',hops_rows),('distance_profiles',dist_rows),('anchor_timings',timings)]:pd.DataFrame(rows).to_csv(out/(name+'.csv'),index=False,float_format='%.17g')
    status['anchors_completed']=len(index);write(out/'status.json',status)
    del result,payload;gc.collect()
   if ordinal==0:write(out/'pilot_complete.json',dict(all_passed=True,anchors=3,each_timepoint_centroid=True,no_model_training=True,source_baseline_reproduction=reproduction,completed_unix=time.time()))
  coordinate_rows=[]
  for tp in case['timepoints']:
   c=contexts[tp];stage('coordinate_counterfactual',timepoint=tp)
   def rebuild(raw_coordinates):
    view=data[c['indices']].copy();view.obsm['spatial']=np.asarray(raw_coordinates)
    _utils._Cal_Spatial_Net(view,k_cutoff=min(cfg['max_neighbors'],cfg['knn_cutoff']*round(view.n_obs/1000)),max_neigh=cfg['max_neighbors'])
    matrix=torch.tensor(view.uns['adj'],dtype=torch.float64);edges,_=dense_to_sparse(matrix);return edges
   result=coordinate_api(model,tp,c['X'],c['Y'],c['edge'],c['raw'],c['scaler'].mean_,c['scaler'].scale_,rebuild,c['whiteners'],baseline_latent=c['z'],atol=tol['cross_device_atol'],rtol=tol['cross_device_rtol'])
   coordinate_rows.extend(dict(timepoint=tp,**row) for row in result['rows']);write(out/'coordinate_responses.json',coordinate_rows)
   for name,values in result['arrays'].items():np.savez_compressed(out/'baseline'/(tp+'_coordinates_'+name+'.npz'),**values)
   assert array_digest(model.svgp_dict[tp].inducing_index_points)==initial_U[tp]
  assert len(index)==60 and state_digest()==before;verify()
  status.update(stage='complete',finished_unix=time.time(),anchors_completed=60);write(out/'status.json',status)
  write(out/'verification.json',dict(all_passed=True,scope_timepoints=case['timepoints'],anchors_per_timepoint=20,anchors_completed=60,exact_basis_VJPs=True,backend='cpu',dtype='float64',optimizer_steps=0,OT_solves=0,learned_weights_unchanged=True,inducing_positions_restored_exactly=True,consumed_files_unchanged=True,checkpoint_reproduction=reproduction,original_source_unchanged=True,native_graphs_match_saved=True,raw_primary_whitened_secondary=True,labels_not_used_to_select_anchors=True,fullgraph_gradient_validation='Tiny native fixture only; real anchorforward checked against fullgraphbaseline',finished_unix=time.time()))
  write(out/'artifact_sha256.json',{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and p.name!='artifact_sha256.json'})
 except Exception as error:
  status.update(stage='failed',error=str(error),traceback=traceback.format_exc(),finished_unix=time.time());write(out/'status.json',status);raise
if __name__=='__main__':main()
