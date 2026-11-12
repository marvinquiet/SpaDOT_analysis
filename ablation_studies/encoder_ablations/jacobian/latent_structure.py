"""Read-only continuous-latent spatial structure; no clustering or neural forward.

Pairs are shared between branches. They are correlated observations, not
independent model replications. Primary normalization removes global latent scale;
whitening is a sensitivity analysis with transforms fitted once to baseline data.
"""
import os
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np

BRANCHES={'GP10':slice(0,10),'GAT10':slice(10,20)}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def array_sha(array):
    a=np.ascontiguousarray(array);h=hashlib.sha256()
    h.update(str(a.dtype).encode());h.update(str(a.shape).encode());h.update(a.tobytes())
    return h.hexdigest()


def unordered_pairs(n,count=100000,seed=20261104):
    assert n>=2 and count>0
    total=n*(n-1)//2
    ids=np.sort(np.random.default_rng(seed).choice(total,size=min(count,total),replace=False))
    row=np.floor(((2*n-1)-np.sqrt((2*n-1)**2-8*ids))/2).astype(np.int64)
    prefix=lambda k:k*(2*n-k-1)//2
    row-=prefix(row)>ids
    row+=prefix(row+1)<=ids
    col=row+1+ids-prefix(row)
    assert np.all((0<=row)&(row<col)&(col<n))
    assert len(set(zip(row.tolist(),col.tolist())))==len(ids)
    return row,col,total


def normalize_latent(z,ridge_relative=1e-6):
    z=np.asarray(z,dtype=np.float64)
    assert z.ndim==2 and z.shape[1]==10 and np.isfinite(z).all()
    mean=z.mean(0);center=z-mean
    population_covariance=center.T@center/len(center)
    trace=float(np.trace(population_covariance));assert trace>0 and np.isfinite(trace)
    covariance=center.T@center/(len(center)-1)
    values,vectors=np.linalg.eigh(covariance)
    ridge=float(np.trace(covariance)/z.shape[1])*ridge_relative
    whitener=(vectors*(1/np.sqrt(np.maximum(values,0)+ridge)))@vectors.T
    white=center@whitener
    white_trace=float(np.mean(np.sum(white*white,axis=1)))
    assert white_trace>0 and np.isfinite(white).all()
    return {'trace_normalized':center/np.sqrt(trace),
            'ridge_whitened':white/np.sqrt(white_trace)},dict(mean=mean.tolist(),
            covariance=covariance.tolist(),primary_population_covariance=population_covariance.tolist(),
            trace=trace,eigenvalues=values.tolist(),ridge=ridge,
            primary_covariance_denominator='N; trace normalization',
            secondary_covariance_denominator='N-1; fixed ridge is relative_ridge times mean sample variance',
            ridge_relative=ridge_relative,whitener=whitener.tolist(),whitened_trace=white_trace,
            fitted_on='unperturbed baseline only; not refit on perturbed data')


def variogram(x,z,pairs,bins=15):
    x=np.asarray(x,dtype=np.float64);i,j=pairs
    assert x.shape==(len(z),2)
    distance=np.linalg.norm(x[i]-x[j],axis=1)
    edges=np.unique(np.quantile(distance,np.linspace(0,1,bins+1)))
    assert len(edges)>=2 and edges[-1]>0
    group=np.minimum(np.searchsorted(edges,distance,side='right')-1,len(edges)-2)
    group=np.maximum(group,0)
    semivariance=.5*np.sum((z[i]-z[j])**2,axis=1)
    assert np.isfinite(semivariance).all()
    rows=[]
    for k in range(len(edges)-1):
        use=group==k
        if not use.any():continue
        vals=semivariance[use]
        rows.append(dict(bin=k,distance_left=float(edges[k]),distance_right=float(edges[k+1]),
            distance_mean=float(distance[use].mean()),pairs=int(use.sum()),
            semivariance_mean=float(vals.mean()),semivariance_pair_SD=float(vals.std(ddof=1)) if len(vals)>1 else None))
    return rows,dict(distance_bin_edges=edges.tolist(),pair_count=len(i),
        pair_distance_sha256=array_sha(distance),pair_index_sha256=array_sha(np.column_stack((i,j))),
        zero_distance_pairs=int(np.sum(distance==0)),binning='shared geometry-quantile bins',
        interpretation='Pair SD describes correlated spot pairs, not uncertainty across independent fits.')


def reference_paths(manifest):
    """Allow explicit diagnostic references or derive pinned parent result paths."""
    ref=manifest.get('reference',{})
    if 'latent' in ref:
        return Path(ref['latent']),ref.get('latent_sha256'),Path(ref.get('input',manifest.get('case',{}).get('input',''))),ref.get('input_sha256')
    case=manifest.get('case')
    if case is None:
        case=next(c for c in manifest['cases'] if c['id']=='AxolotlTelencephalon')
    result=manifest.get('reference_result_dir') or manifest.get('parent_result')
    if result is None:
        parent=manifest.get('parent_root') or manifest.get('reference_root')
        result=str(Path(parent)/'results/AxolotlTelencephalon__dual_native__seed2026')
    result=Path(result)
    complete=json.loads((result/'complete.json').read_text())
    return result/'latent.h5ad',complete['files']['latent.h5ad'],Path(case['input']),case['input_sha256']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--destination',type=Path,required=True)
    p.add_argument('--pairs',type=int,default=100000);p.add_argument('--bins',type=int,default=15)
    args=p.parse_args();root=args.root.resolve();destination=args.destination.resolve()
    manifest_path=root/'manifest.json';manifest=json.loads(manifest_path.read_text())
    latent_path,latent_hash,input_path,input_hash=reference_paths(manifest)
    assert latent_path.exists() and input_path.exists()
    assert latent_hash and input_hash,'Trusted input/latent SHA256 required'
    protected={str(manifest_path):sha(manifest_path),str(latent_path):latent_hash,str(input_path):input_hash}
    assert all(sha(k)==v for k,v in protected.items())
    assert not destination.exists();destination.mkdir(parents=True)
    import anndata as ad,pandas as pd
    source=ad.read_h5ad(input_path,backed='r');latent=ad.read_h5ad(latent_path)
    try:
        assert latent.n_vars==20 and latent.obs_names.is_unique and source.obs_names.is_unique
        assert set(latent.obs_names)==set(source.obs_names)
        aligned=source[latent.obs_names]
        assert np.array_equal(aligned.obsm['spatial'],latent.obsm['spatial'])
        assert aligned.obs.timepoint.astype(str).tolist()==latent.obs.timepoint.astype(str).tolist()
    finally:source.file.close()
    z=np.asarray(latent.X,dtype=np.float64);assert np.isfinite(z).all()
    labels=latent.obs.timepoint.astype(str).to_numpy();raw=np.asarray(latent.obsm['spatial'],dtype=np.float64)
    timepoints=['t0','t1','t2'];assert sorted(set(labels))==timepoints
    rows=[];summary=[];transforms={};pair_records={};standardized=np.empty_like(raw)
    for ix,tp in enumerate(timepoints):
        sel=labels==tp;coords=raw[sel];scale=coords.std(0);assert (scale>0).all()
        x=(coords-coords.mean(0))/scale;standardized[sel]=x
        i,j,total=unordered_pairs(len(x),args.pairs,seed=20261104+ix)
        for branch,columns in BRANCHES.items():
            normalized,transform=normalize_latent(z[sel,columns]);transforms[tp+'__'+branch]=transform
            for norm,values in normalized.items():
                records,provenance=variogram(x,values,(i,j),args.bins)
                rows.extend(dict(timepoint=tp,branch=branch,normalization=norm,**r) for r in records)
                pair_records[tp]=dict(provenance,total_possible_unordered_pairs=total,
                    sampling_seed=20261104+ix,shared_pair_ids_between_branches=True,
                    coordinate_mean=coords.mean(0).tolist(),coordinate_std_population=scale.tolist())
                summary.append(dict(timepoint=tp,branch=branch,normalization=norm,n_spots=int(sel.sum()),
                    sampled_pairs=len(i),global_pair_semivariance=float(np.mean(.5*np.sum((values[i]-values[j])**2,axis=1))),
                    first_bin_semivariance=records[0]['semivariance_mean'],
                    last_bin_semivariance=records[-1]['semivariance_mean']))
    pd.DataFrame(rows).to_csv(destination/'variogram.csv',index=False,float_format='%.17g')
    pd.DataFrame(summary).to_csv(destination/'structure_summary.csv',index=False,float_format='%.17g')
    np.savez_compressed(destination/'display_geometry.npz',spot=latent.obs_names.to_numpy(dtype=str),
        timepoint=labels.astype(str),coordinates_raw=raw,coordinates_standardized=standardized)
    (destination/'fixed_latent_transforms.json').write_text(json.dumps(transforms,indent=2,allow_nan=False)+'\n')
    (destination/'pair_sampling.json').write_text(json.dumps(pair_records,indent=2,allow_nan=False)+'\n')
    assert all(sha(k)==v for k,v in protected.items())
    verification=dict(all_passed=True,source_manifest_sha256=protected[str(manifest_path)],
        source_latent_sha256=latent_hash,input_sha256=input_hash,timepoints=timepoints,
        branches=list(BRANCHES),primary_normalization='trace_normalized',
        sensitivity_normalization='ridge_whitened',coordinates='population-standardized separately per timepoint, source preprocessing convention',
        no_neural_forward=True,no_model_training=True,no_clustering=True,no_OT=True,no_annotations_used=True,
        native_graph_energy_performed=False,native_graph_energy_note='No exact native edge file consumed; no substitute graph was invented.',
        pair_count_per_timepoint=args.pairs,pair_dependence_disclosed=True,
        interpretation='Continuous latent spatial structure of one fixed jointly trained model; not architecture necessity or predictive benefit.',
        visual_QA_performed=False,finished_unix=time.time())
    (destination/'verification.json').write_text(json.dumps(verification,indent=2,allow_nan=False)+'\n')
    artifacts={str(f.relative_to(destination)):sha(f) for f in destination.rglob('*') if f.is_file()}
    (destination/'artifact_sha256.json').write_text(json.dumps(artifacts,indent=2)+'\n')
    print(json.dumps(verification),flush=True)


if __name__=='__main__':main()
