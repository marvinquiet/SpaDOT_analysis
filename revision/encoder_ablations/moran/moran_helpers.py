"""Actual distance-bin Moran statistics, not relabelled semivariance.

Read-only CPU analysis of the same fixed latent, pairs and bins as the earlier
variogram. No neural model, optimizer, clustering or transport is imported.
"""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def weight_statistics(z,i,j):
    """Moran I for symmetric binary weights on selected unordered pairs.

    W_ij=W_ji=1 for the supplied pairs and S0=2m. The trace statistic is the
    variance-weighted mean of coordinate-wise Moran statistics, not a scalar
    spatially standardized Pearson coefficient and not generally bounded ±1.
    Input must already be globally centered within a timepoint.
    """
    z=np.asarray(z,dtype=np.float64);i=np.asarray(i,dtype=np.int64);j=np.asarray(j,dtype=np.int64)
    assert z.ndim==2 and np.isfinite(z).all() and len(z)>1
    assert i.shape==j.shape and i.ndim==1 and len(i)>0
    assert np.all((i>=0)&(i<j)&(j<len(z)))
    assert np.allclose(z.mean(0),0,atol=1e-10,rtol=0),'Center each branch globally before Moran I'
    n,m=len(z),len(i)
    denominator=np.sum(z*z,axis=0);total=float(denominator.sum())
    assert np.isfinite(total) and total>0
    products=np.sum(z[i]*z[j],axis=0)
    valid=denominator>0
    individual=np.full(z.shape[1],np.nan)
    individual[valid]=n*products[valid]/(m*denominator[valid])
    trace=float(n*products.sum()/(m*total))
    variances=denominator/total
    assert np.isclose(trace,np.sum(variances[valid]*individual[valid]),atol=1e-12,rtol=1e-12)
    scale=n/total
    gamma=float(scale*.5*np.mean(np.sum((z[i]-z[j])**2,axis=1)))
    endpoint=float(scale*.5*np.mean(np.sum(z[i]**2+z[j]**2,axis=1)))
    assert np.isclose(gamma,endpoint-trace,atol=2e-12,rtol=2e-12)
    return dict(Moran_I_trace=trace,Moran_I_unweighted_dimension_mean=float(individual[valid].mean()),
        normalized_semivariance=gamma,endpoint_energy=endpoint,
        one_minus_semivariance=1-gamma,
        difference_from_one_minus_semivariance=trace-(1-gamma),
        sampled_unordered_pairs=m,weight_sum_S0=2*m,n_spots=n,
        expected_randomization_I=-1/(n-1),
        dimension_I=individual,dimension_variance_weights=variances,
        valid_dimensions=valid,global_mean_squared_norm=total/n)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--structure',type=Path)
    parser.add_argument('--destination',type=Path,required=True)
    args=parser.parse_args();root=args.root.resolve();structure=(args.structure or root/'structure').resolve()
    destination=args.destination.resolve();assert not destination.exists()
    sys.path.insert(0,str(root))
    import latent_structure as original
    import pandas as pd
    manifest_path=root/'manifest.json';manifest=read(manifest_path)
    assert sha(root/'latent_structure.py')==manifest['code_sha256']['latent_structure.py']
    structure_hashes=read(structure/'artifact_sha256.json')
    for name,digest in structure_hashes.items():assert sha(structure/name)==digest,name
    previous=read(structure/'verification.json');assert previous['all_passed']
    latent_path,latent_hash,input_path,input_hash=original.reference_paths(manifest)
    assert latent_hash==previous['source_latent_sha256'] and input_hash==previous['input_sha256']
    protected={str(manifest_path):sha(manifest_path),str(root/'latent_structure.py'):sha(root/'latent_structure.py'),
        str(latent_path):latent_hash,str(input_path):input_hash,
        str(structure/'artifact_sha256.json'):sha(structure/'artifact_sha256.json')}
    for name,digest in protected.items():assert sha(name)==digest,name
    import anndata as ad
    latent=ad.read_h5ad(latent_path)
    assert latent.n_vars==20 and latent.obs_names.is_unique and np.isfinite(latent.X).all()
    z=np.asarray(latent.X,dtype=np.float64)
    raw=np.asarray(latent.obsm['spatial'],dtype=np.float64)
    labels=latent.obs.timepoint.astype(str).to_numpy()
    old_table=pd.read_csv(structure/'variogram.csv',float_precision='round_trip')
    old_pairs=read(structure/'pair_sampling.json')
    rows=[];dimensions=[];summaries=[];identity_errors=[];provenance={}
    for tp_index,tp in enumerate(('t0','t1','t2')):
        selected=labels==tp;coordinates=raw[selected]
        assert selected.sum()>1
        coordinates=(coordinates-coordinates.mean(0))/coordinates.std(0)
        i,j,total=original.unordered_pairs(len(coordinates),old_pairs[tp]['pair_count'],20261104+tp_index)
        distances=np.linalg.norm(coordinates[i]-coordinates[j],axis=1)
        assert original.array_sha(np.column_stack([i,j]))==old_pairs[tp]['pair_index_sha256']
        assert original.array_sha(distances)==old_pairs[tp]['pair_distance_sha256']
        edges=np.asarray(old_pairs[tp]['distance_bin_edges'])
        group=np.minimum(np.searchsorted(edges,distances,side='right')-1,len(edges)-2)
        group=np.maximum(group,0)
        provenance[tp]=dict(pair_index_sha256=original.array_sha(np.column_stack([i,j])),
            pair_distance_sha256=original.array_sha(distances),sampling_seed=20261104+tp_index,
            distance_bin_edges=edges.tolist(),unordered_pairs=len(i),total_possible_pairs=total,
            same_pairs_and_bins_as_original_variogram=True)
        for branch,columns in original.BRANCHES.items():
            normalized,_=original.normalize_latent(z[selected,columns])
            for normalization,values in normalized.items():
                assert np.isclose(np.mean(np.sum(values*values,axis=1)),1,atol=1e-12,rtol=1e-12)
                for bin_index in range(len(edges)-1):
                    use=group==bin_index
                    if not use.any():continue
                    statistics=weight_statistics(values,i[use],j[use])
                    old=old_table[(old_table.timepoint==tp)&(old_table.branch==branch)&
                        (old_table.normalization==normalization)&(old_table.bin==bin_index)]
                    assert len(old)==1
                    old=old.iloc[0]
                    assert int(old.pairs)==int(use.sum())
                    assert float(old.distance_mean)==float(distances[use].mean())
                    err=abs(statistics['normalized_semivariance']-float(old.semivariance_mean))
                    identity_errors.append(err)
                    assert err<2e-12,(tp,branch,normalization,bin_index,err)
                    base=dict(timepoint=tp,branch=branch,normalization=normalization,bin=bin_index,
                        distance_left=float(edges[bin_index]),distance_right=float(edges[bin_index+1]),
                        distance_mean=float(distances[use].mean()))
                    arrays={name:statistics.pop(name) for name in
                        ('dimension_I','dimension_variance_weights','valid_dimensions')}
                    rows.append(dict(base,**statistics))
                    for dim in range(10):
                        dimensions.append(dict(base,dimension=dim,
                            Moran_I=float(arrays['dimension_I'][dim]) if arrays['valid_dimensions'][dim] else None,
                            variance_weight=float(arrays['dimension_variance_weights'][dim]),
                            status='defined' if arrays['valid_dimensions'][dim] else 'undefined_constant_dimension'))
                whole=weight_statistics(values,i,j)
                summaries.append(dict(timepoint=tp,branch=branch,normalization=normalization,
                    n_spots=len(values),sampled_pairs=len(i),all_sampled_pairs_Moran_I_trace=whole['Moran_I_trace'],
                    first_bin_Moran_I_trace=next(r['Moran_I_trace'] for r in rows if r['timepoint']==tp and
                        r['branch']==branch and r['normalization']==normalization),
                    expected_randomization_I=-1/(len(values)-1)))
    destination.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(destination/'moran_by_distance.csv',index=False,float_format='%.17g')
    pd.DataFrame(dimensions).to_csv(destination/'moran_per_dimension.csv',index=False,float_format='%.17g')
    pd.DataFrame(summaries).to_csv(destination/'moran_summary.csv',index=False,float_format='%.17g')
    (destination/'pair_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    for name,digest in protected.items():assert sha(name)==digest,name
    verification=dict(all_passed=True,source_manifest_sha256=sha(manifest_path),
        source_structure_artifact_sha256=sha(structure/'artifact_sha256.json'),
        source_latent_sha256=latent_hash,input_sha256=input_hash,
        statistic='Variance-weighted average of coordinate-wise Moran I; equivalent trace cross-product statistic',
        formula='n/S0 * sum_ij(w_ij z_i dot z_j) / sum_i(z_i dot z_i)',
        weights='Symmetric binary weights for the SAME sampled unordered pairs in each distance bin; no row standardization',
        pair_count_per_timepoint=100000,distance_bins=15,
        primary_normalization='trace_normalized',secondary_normalization='ridge_whitened',
        per_dimension_unweighted_mean_secondary=True,computed_from_crossproducts=True,
        not_one_minus_semivariance=True,semivariance_match_max_error=max(identity_errors),
        empirical_identity='normalized_semivariance = endpoint_energy - Moran_I_trace',
        same_pairs_bins_and_normalizations_verified=True,
        stationarity_not_assumed_for_descriptive_statistic=True,
        no_permutation_pvalues=True,no_independent_pair_uncertainty_claim=True,
        no_training=True,no_neural_forward=True,no_OT=True,no_clustering=True,
        protected_files=protected,visual_QA_performed=False,finished_unix=time.time())
    (destination/'verification.json').write_text(json.dumps(verification,indent=2,allow_nan=False)+'\n')
    artifacts={p.name:sha(p) for p in destination.iterdir() if p.is_file()}
    (destination/'artifact_sha256.json').write_text(json.dumps(artifacts,indent=2)+'\n')
    print(json.dumps(verification),flush=True)


if __name__=='__main__':main()
