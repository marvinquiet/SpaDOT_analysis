"""Render verified continuous-latent diagnostics; missing mechanism sections stay pending.

No neural model, training, clustering, or transport solve is imported. PDFs must
be visually reviewed after rendering; this script does not certify visual QA.
"""
import os
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
from pathlib import Path
import argparse,hashlib,html,json,time
import numpy as np
import pandas as pd

COLORS={'GP10':'#2458a6','GAT10':'#d46a19'}
LABELS={'GP10':'SVGP','GAT10':'GAT'}
STAGES={'t0':'St.44','t1':'St.54','t2':'St.57'}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def checked_structure(folder):
    folder=Path(folder)
    v=json.loads((folder/'verification.json').read_text());assert v['all_passed']
    hashes=json.loads((folder/'artifact_sha256.json').read_text())
    for name,digest in hashes.items():assert sha(folder/name)==digest,name
    table=pd.read_csv(folder/'variogram.csv',float_precision='round_trip')
    summary=pd.read_csv(folder/'structure_summary.csv',float_precision='round_trip')
    expected={(tp,b,norm) for tp in STAGES for b in LABELS for norm in ['trace_normalized','ridge_whitened']}
    assert set(map(tuple,summary[['timepoint','branch','normalization']].to_numpy()))==expected
    assert not table.duplicated(['timepoint','branch','normalization','bin']).any()
    assert np.isfinite(table[['distance_mean','semivariance_mean']].to_numpy()).all()
    for tp in STAGES:
        for norm in ['trace_normalized','ridge_whitened']:
            parts=[table[(table.timepoint==tp)&(table.branch==b)&(table.normalization==norm)].sort_values('bin') for b in LABELS]
            assert np.array_equal(parts[0].bin,parts[1].bin)
            assert np.array_equal(parts[0].distance_mean,parts[1].distance_mean)
            assert np.array_equal(parts[0].pairs,parts[1].pairs)
    return table,summary,v,hashes


def render_variogram(table,summary,destination,normalization):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14,'axes.titlesize':18,
        'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':13,
        'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(15.8,5.2),sharey=True)
    for ax,tp in zip(axes,STAGES):
        for branch in LABELS:
            p=table[(table.timepoint==tp)&(table.branch==branch)&(table.normalization==normalization)].sort_values('bin')
            ax.plot(p.distance_mean,p.semivariance_mean,color=COLORS[branch],label=LABELS[branch],
                    marker='o',markersize=4.6,lw=2.4)
        n=int(summary[(summary.timepoint==tp)&(summary.normalization==normalization)].n_spots.iloc[0])
        ax.set_title(STAGES[tp]+'  |  n = '+f'{n:,}')
        ax.set_xlabel('Standardized spatial distance')
        ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.18,lw=.7)
    axes[0].set_ylabel('Normalized latent semivariance')
    axes[2].legend(frameon=False,loc='best')
    # Shared-axis limits must be set after every panel contributes its values;
    # setting only bottom inside the first panel otherwise freezes its upper
    # autoscale limit and can clip a later timepoint.
    maximum=float(table[table.normalization==normalization].semivariance_mean.max())
    axes[0].set_ylim(0,maximum*1.06)
    note=('Primary: centered branch representations normalized by total latent variance.' if normalization=='trace_normalized'
          else 'Sensitivity analysis: baseline covariance whitening with a fixed ridge, followed by variance normalization.')
    fig.text(.5,.065,note,ha='center',fontsize=12,color='#334155')
    fig.text(.5,.025,'The same 100,000 spot pairs are used for both branches at each timepoint; pairs are not independent model replicates.',
             ha='center',fontsize=11,color='#64748b')
    fig.subplots_adjust(left=.075,right=.985,top=.87,bottom=.24,wspace=.20)
    stem='latent_spatial_variogram' if normalization=='trace_normalized' else 'latent_spatial_variogram_whitened'
    paths=[]
    for suffix in ['png','pdf','svg']:
        p=destination/(stem+'.'+suffix)
        fig.savefig(p,dpi=200,facecolor='white');paths.append(p)
    plt.close(fig)
    return paths


def checked_mechanism(folder,structure):
    folder=Path(folder);v=json.loads((folder/'verification.json').read_text())
    assert v.get('all_passed') is True,'Only verified completed mechanism results can be rendered'
    hashes=json.loads((folder/'artifact_sha256.json').read_text())
    for name,digest in hashes.items():assert sha(folder/name)==digest,name
    index=json.loads((folder/'influence_index.json').read_text())
    assert isinstance(index,list) and index and {x['timepoint'] for x in index}==set(STAGES)
    assert len({(x['timepoint'],int(x['anchor'])) for x in index})==len(index)
    arrays=[]
    geometry=np.load(Path(structure)/'display_geometry.npz',allow_pickle=False)
    for entry in index:
        with np.load(folder/entry['npz'],allow_pickle=False) as z:data={k:z[k].copy() for k in z.files}
        needed={'GP_raw_energy','GAT_raw_energy','GP_whitened_energy','GAT_whitened_energy',
                'hops','coordinates_raw','coordinates_standardized','spot'}
        assert needed<=set(data)
        n=len(data['spot']);a=int(entry['anchor'])
        assert len(set(data['spot'].tolist()))==n and 0<=a<n
        assert str(data['spot'][a])==str(entry['anchor_spot'])
        for key in ['GP_raw_energy','GAT_raw_energy','GP_whitened_energy','GAT_whitened_energy']:
            assert data[key].shape==(n,) and np.isfinite(data[key]).all() and (data[key]>=0).all()
        assert data['coordinates_raw'].shape==data['coordinates_standardized'].shape==(n,2)
        assert data['hops'].shape==(n,) and data['hops'][a]==0
        assert np.all(data['GAT_raw_energy'][(data['hops']>3)|(data['hops']<0)]==0)
        sel=geometry['timepoint']==entry['timepoint']
        lookup={str(s):i for i,s in enumerate(geometry['spot'][sel])}
        positions=[lookup[str(s)] for s in data['spot']]
        assert np.array_equal(geometry['coordinates_raw'][sel][positions],data['coordinates_raw'])
        arrays.append((entry,data))
    geometry.close()
    hop=pd.read_csv(folder/'hop_profiles.csv',float_precision='round_trip')
    distance=pd.read_csv(folder/'distance_profiles.csv',float_precision='round_trip')
    summary=pd.read_csv(folder/'influence_summary.csv',float_precision='round_trip')
    for frame in [hop,distance,summary]:
        assert set(frame.timepoint)==set(STAGES) and set(frame.branch)=={'GP','GAT'}
        assert set(frame.metric)=={'raw','whitened'}
    coord=json.loads((folder/'coordinate_responses.json').read_text())
    assert coord and {row['timepoint'] for row in coord}==set(STAGES)
    return arrays,hop,distance,summary,coord,v


def save_figure(fig,dest,stem):
    import matplotlib.pyplot as plt
    paths=[]
    for ext in ['png','pdf','svg']:
        path=dest/(stem+'.'+ext);fig.savefig(path,dpi=200,facecolor='white');paths.append(path)
    plt.close(fig);return paths


def render_influence_maps(arrays,dest):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    selected=[]
    for tp in STAGES:
        # The first predeclared target is the geometry-centroid target. Never
        # select a displayed target by its influence result or anatomical label.
        candidates=[(e,d) for e,d in arrays if e['timepoint']==tp]
        selected.append(min(candidates,key=lambda x:float(np.linalg.norm(
            x[1]['coordinates_standardized'][int(x[0]['anchor'])]-
            x[1]['coordinates_standardized'].mean(axis=0)))))
    distributions=[]
    for entry,data in selected:
        for branch in ['GP','GAT']:
            energy=data[branch+'_raw_energy'].copy();energy[int(entry['anchor'])]=0
            distributions.append(energy/energy.sum() if energy.sum()>0 else energy)
    positive=np.concatenate([x[x>0] for x in distributions])
    vmax=float(positive.max()) if len(positive) else 1.
    vmin=min(1e-7,vmax/1000)
    cmap=plt.get_cmap('magma').copy();cmap.set_bad('#e5e7eb')
    norm=LogNorm(vmin=vmin,vmax=vmax)
    fig,axes=plt.subplots(2,3,figsize=(14.7,9.2))
    plot=None
    for col,(entry,data) in enumerate(selected):
        coords=data['coordinates_raw'];a=int(entry['anchor'])
        for row,branch in enumerate(['GP','GAT']):
            ax=axes[row,col];energy=data[branch+'_raw_energy'];total=energy.sum()
            selfmass=float(energy[a]/total) if total>0 else np.nan
            nonself=energy.copy();nonself[a]=0
            weights=nonself/nonself.sum() if nonself.sum()>0 else nonself
            # Scatter masks may omit zero-valued locations entirely. Draw the
            # full tissue in gray first, then overlay positive sensitivity.
            ax.scatter(coords[:,0],coords[:,1],c='#e5e7eb',s=5,
                       edgecolors='none',rasterized=True)
            plot=ax.scatter(coords[:,0],coords[:,1],c=np.ma.masked_where(weights<=0,weights),
                s=5,cmap=cmap,norm=norm,edgecolors='none',rasterized=True)
            ax.scatter(*coords[a],marker='*',s=95,c='#41b6c4',edgecolor='black',linewidth=.7,zorder=4)
            ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
            for spine in ax.spines.values():spine.set_visible(False)
            ax.set_title(('SVGP' if branch=='GP' else 'GAT')+' | '+STAGES[entry['timepoint']],fontsize=18)
            ax.text(.5,-.04,'Self-energy fraction: '+(f'{selfmass:.3f}' if np.isfinite(selfmass) else 'undefined'),
                    transform=ax.transAxes,ha='center',fontsize=12)
    cax=fig.add_axes([.93,.19,.018,.63]);bar=fig.colorbar(plot,cax=cax)
    bar.set_label('Normalized nonself sensitivity energy (log scale)',fontsize=12)
    fig.text(.48,.045,'Same preselected target for both branches; star marks the target. Gray = zero input sensitivity.',ha='center',fontsize=12)
    fig.text(.48,.018,'Primary metric: squared Frobenius Jacobian energy, normalized after removing the target itself.',ha='center',fontsize=11,color='#64748b')
    fig.subplots_adjust(left=.025,right=.90,bottom=.11,top=.95,hspace=.23,wspace=.08)
    return save_figure(fig,dest,'trained_encoder_influence_maps')


def render_hop_profiles(hop,dest):
    import matplotlib.pyplot as plt
    labels=['1','2','3','>3_or_unreachable'];names=['1 hop','2 hops','3 hops','>3 / unreachable']
    primary=hop[hop.metric=='raw'].copy();primary['bin']=primary['bin'].astype(str)
    fig,axes=plt.subplots(2,3,figsize=(15.8,8.8))
    for col,tp in enumerate(STAGES):
        part=primary[primary.timepoint==tp]
        n=int(part.anchor.nunique())
        for row,metric in enumerate(['nonself_mass','mean_nonself_mass_per_source']):
            ax=axes[row,col]
            for branch,color in [('GP',COLORS['GP10']),('GAT',COLORS['GAT10'])]:
                p=part[part.branch==branch]
                means=[pd.to_numeric(p[p['bin']==label][metric],errors='coerce').mean() for label in labels]
                ax.plot(np.arange(4),means,'o-',color=color,lw=2.1,label='SVGP' if branch=='GP' else 'GAT')
            ax.set_xticks(np.arange(4));ax.set_xticklabels(names,rotation=18,ha='right',fontsize=11)
            ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.18)
            if row==0:ax.set_title(STAGES[tp]+f' | {n} targets',fontsize=17)
            if col==0:ax.set_ylabel('Total nonself influence fraction' if row==0 else 'Mean influence fraction per source spot',fontsize=12)
    axes[0,2].legend(frameon=False)
    fig.text(.5,.025,'Each target is normalized separately. Curves average across the same spatially selected targets, not independent training seeds.',ha='center',fontsize=11,color='#475569')
    fig.subplots_adjust(left=.08,right=.985,bottom=.14,top=.91,wspace=.25,hspace=.43)
    return save_figure(fig,dest,'influence_by_directed_hop')


def render_range_profiles(distance,summary,dest):
    import matplotlib.pyplot as plt
    primary=distance[distance.metric=='raw'].copy();s=summary[summary.metric=='raw'].copy()
    fig,axes=plt.subplots(3,3,figsize=(15.8,11.4))
    for col,tp in enumerate(STAGES):
        for branch,color in [('GP',COLORS['GP10']),('GAT',COLORS['GAT10'])]:
            part=primary[(primary.timepoint==tp)&(primary.branch==branch)]
            curves=[]
            for _,a in part.groupby('anchor'):
                a=a.sort_values('lower');curves.append(a.nonself_mass.to_numpy().cumsum())
            edges=part[['lower','upper']].drop_duplicates().sort_values('lower')
            assert curves and all(len(c)==len(edges) for c in curves)
            axes[0,col].plot(edges.upper,np.nanmean(np.stack(curves),axis=0),color=color,lw=2.3,
                             label='SVGP' if branch=='GP' else 'GAT')
        axes[0,col].set_title(STAGES[tp],fontsize=18);axes[0,col].set_ylim(0,1.025)
        axes[0,col].set_xlabel('Standardized spatial distance',fontsize=12)
        for row,metric in [(1,'nonself_radius90'),(2,'self_mass')]:
            part=s[s.timepoint==tp]
            pivot=part.pivot(index='anchor',columns='branch',values=metric)
            for _,target in pivot.iterrows():
                if np.isfinite(target.get('GP',np.nan)) and np.isfinite(target.get('GAT',np.nan)):
                    axes[row,col].plot([0,1],[target.GP,target.GAT],color='#aeb7c4',lw=.8,alpha=.6,zorder=1)
            for ix,branch in enumerate(['GP','GAT']):
                values=part[part.branch==branch][metric].dropna().to_numpy()
                axes[row,col].scatter(np.full(len(values),ix),values,c=COLORS['GP10' if branch=='GP' else 'GAT10'],s=24,zorder=3)
                if len(values):axes[row,col].plot([ix-.13,ix+.13],[np.median(values)]*2,c='black',lw=2.2,zorder=4)
            axes[row,col].set_xticks([0,1]);axes[row,col].set_xticklabels(['SVGP','GAT'])
            axes[row,col].set_xlim(-.4,1.4);axes[row,col].set_ylim(bottom=0)
            if metric=='self_mass':axes[row,col].set_ylim(0,1.025)
        for ax in axes[:,col]:ax.grid(axis='y',alpha=.17)
    axes[0,0].set_ylabel('Cumulative nonself influence fraction',fontsize=12)
    axes[1,0].set_ylabel('90% nonself influence radius\n(standardized distance)',fontsize=12)
    axes[2,0].set_ylabel('Self-energy fraction',fontsize=12)
    axes[0,2].legend(frameon=False)
    fig.text(.5,.025,'Paired points are target locations from one fixed model. Black marks show medians; no seed-level confidence intervals are implied.',ha='center',fontsize=11,color='#475569')
    fig.subplots_adjust(left=.09,right=.985,bottom=.09,top=.94,wspace=.25,hspace=.44)
    return save_figure(fig,dest,'influence_distance_and_self_fraction')


def render_coordinate_response(rows,transforms,dest):
    import matplotlib.pyplot as plt
    order=['scale_0.5','identity','scale_2','rigid_reflection'];names=['0.5x','1x','2x','Reflection']
    fig,axes=plt.subplots(1,3,figsize=(15.8,5.4))
    for ax,tp in zip(axes,STAGES):
        part={r['transform']:r for r in rows if r['timepoint']==tp}
        assert set(part)==set(order)
        for branch,color in [('GP',COLORS['GP10']),('GAT',COLORS['GAT10'])]:
            divisor=np.sqrt(transforms[tp+'__'+branch+'10']['trace']/10)
            values=[part[name][branch]['rms_raw']/divisor if part[name]['status']=='evaluated' else np.nan for name in order]
            ax.plot(range(4),values,'o-',color=color,lw=2.3,label='SVGP' if branch=='GP' else 'GAT')
        for ix,name in enumerate(order):
            if part[name]['status']!='evaluated':
                ax.text(ix,.03,'Graph changed',rotation=90,transform=ax.get_xaxis_transform(),ha='center',va='bottom',color='#64748b',fontsize=10)
        ax.set_xticks(range(4));ax.set_xticklabels(names,rotation=15,ha='right')
        ax.set_title(STAGES[tp],fontsize=18);ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.18)
    axes[0].set_ylabel('RMS latent change / baseline RMS variation',fontsize=12)
    axes[2].legend(frameon=False)
    fig.text(.5,.07,'Input and inducing locations undergo the same transform; learned weights, original scaler and kernel scale remain fixed.',ha='center',fontsize=11)
    fig.text(.5,.025,'Only exactly unchanged native graphs are evaluated. Identity/reflection are invariance controls, not task-performance tests.',ha='center',fontsize=11,color='#475569')
    fig.subplots_adjust(left=.08,right=.985,bottom=.26,top=.87,wspace=.25)
    return save_figure(fig,dest,'coordinate_geometry_response')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--structure',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    parser.add_argument('--mechanism',type=Path)
    args=parser.parse_args();dest=args.destination.resolve();assert not dest.exists();dest.mkdir(parents=True)
    table,summary,verification,source_hashes=checked_structure(args.structure)
    paths=[]
    for norm in ['trace_normalized','ridge_whitened']:paths.extend(render_variogram(table,summary,dest,norm))
    mechanism_done=False;mechanism_provenance=None
    if args.mechanism is not None:
        arrays,hop,distance,influence_summary,coordinates,mechanism_verification=checked_mechanism(args.mechanism,args.structure)
        paths.extend(render_influence_maps(arrays,dest))
        paths.extend(render_hop_profiles(hop,dest))
        paths.extend(render_range_profiles(distance,influence_summary,dest))
        transforms=json.loads((args.structure/'fixed_latent_transforms.json').read_text())
        paths.extend(render_coordinate_response(coordinates,transforms,dest))
        for filename in ['influence_summary.csv','hop_profiles.csv','distance_profiles.csv','coordinate_responses.json']:
            (dest/filename).write_bytes((args.mechanism/filename).read_bytes())
        mechanism_done=True
        mechanism_provenance=sha(args.mechanism/'artifact_sha256.json')
    summary.to_csv(dest/'structure_summary.csv',index=False,float_format='%.17g')
    table.to_csv(dest/'variogram.csv',index=False,float_format='%.17g')
    small=summary[summary.normalization=='trace_normalized'][['timepoint','branch','n_spots','sampled_pairs','first_bin_semivariance','last_bin_semivariance']].copy()
    small['timepoint']=small.timepoint.map(STAGES);small['branch']=small.branch.map(LABELS)
    small.columns=['Stage','Branch','Spots','Shared pairs','Closest-distance bin','Farthest-distance bin']
    table_html=small.to_html(index=False,float_format=lambda x:f'{x:.4f}',border=0)
    pending='Pending: exact trained-model influence and coordinate-response diagnostics have not been consumed by this structure-only report.'
    if mechanism_done:
        pending='Completed: paired trained-model input sensitivities, unchanged-graph coordinate responses, and continuous-latent spatial structure. Numerical validation passed; figure visual QA is recorded separately.'
    page='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Axolotl encoder mechanisms</title><style>
body{margin:0;background:#f5f7fa;color:#172033;font:16px/1.55 -apple-system,BlinkMacSystemFont,Segoe UI,sans-serif}main{max-width:1380px;margin:32px auto;padding:0 26px}h1{font-size:30px;margin-bottom:8px}h2{font-size:22px}.card{background:white;border:1px solid #e1e6ee;border-radius:12px;padding:22px;margin:20px 0}img{width:100%;height:auto}.pending{border-left:5px solid #b7791f;background:#fff9ec;padding:13px 18px}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}th,td{text-align:left;padding:10px;border-bottom:1px solid #e5e7eb}small{color:#64748b}a{color:#2458a6}code{word-break:break-all}</style><main>
<h1>Axolotl telencephalon: SVGP and GAT spatial structure</h1><p>One fixed native SpaDOT model, training seed 2026, executed epoch 99. Both branches are jointly trained 10-dimensional posterior-mean representations.</p>
<div class="pending">'''+html.escape(pending)+'''</div>
<section class="card"><h2>Continuous latent spatial variation</h2><p>Each curve measures how different two latent representations are as their physical separation increases. Values are divided by total branch latent variance to remove a global scale difference. Lower short-distance semivariance indicates greater local smoothness; it does not establish more accurate anatomical recovery.</p><img src="latent_spatial_variogram.png" alt="Paired SVGP and GAT normalized semivariograms for three Axolotl stages"><p><a href="latent_spatial_variogram.pdf">PDF</a> · <a href="variogram.csv">Full values</a></p>'''+table_html+'''</section>
<section class="card"><h2>Whitening sensitivity analysis</h2><p>Whitening gives latent directions more equal weight using the baseline covariance and a fixed numerical ridge. It is supplementary to the primary total-variance normalization.</p><img src="latent_spatial_variogram_whitened.png" alt="Whitening sensitivity analysis"><p><a href="latent_spatial_variogram_whitened.pdf">PDF</a></p></section>
<section class="card"><h2>Scope and interpretation</h2><ul><li>No model training, neural forward, clustering or OT solve was performed for these curves.</li><li>Exactly the same sampled pairs and distance bins were used for both branches at each stage.</li><li>Spot pairs are correlated; the pair count is not a sample size for seed-level significance.</li><li>This is a description of one saved model on development data, not a test of universal encoder necessity.</li></ul><p><small>Source latent SHA256: <code>'''+html.escape(verification['source_latent_sha256'])+'''</code></small></p></section></main></html>'''
    if mechanism_done:
        sections='''<section class="card"><h2>Trained-model input influence</h2><p>Each map uses the same predeclared geometric target for both branches. Colors show squared Frobenius expression-to-latent Jacobian energy, normalized after removing self influence. The target is marked by a star; self influence remains reported separately. Maps describe one target, whereas the following profiles use every analyzed target.</p><img src="trained_encoder_influence_maps.png" alt="Paired SVGP/GAT trained-input sensitivity maps"><p><a href="trained_encoder_influence_maps.pdf">PDF</a></p></section>
<section class="card"><h2>Influence by directed graph distance</h2><p>Both total influence mass and mean mass per source spot are shown, so a larger far-neighbor population is not confused with stronger individual influence. Three GAT layers constrain fixed-weight expression sensitivity to incoming paths of at most three hops; they do not guarantee strongest one-hop influence.</p><img src="influence_by_directed_hop.png" alt="Directed hop influence profiles"><p><a href="influence_by_directed_hop.pdf">PDF</a> · <a href="hop_profiles.csv">Values</a></p></section>
<section class="card"><h2>Influence range and self sensitivity</h2><p>Profiles and paired target points describe one fixed model, not independent training-seed replications. Black marks indicate target medians. The 90% radius excludes the target itself.</p><img src="influence_distance_and_self_fraction.png" alt="Influence range and self sensitivity"><p><a href="influence_distance_and_self_fraction.pdf">PDF</a> · <a href="influence_summary.csv">Values</a></p></section>
<section class="card"><h2>Continuous spatial geometry response</h2><p>Observed and inducing locations receive the same prescribed transformation. Learned weights, the baseline coordinate standardization, and the kernel scale remain fixed. A transformation is evaluated only if the rebuilt native graph is exactly unchanged. This verifies geometry dependence; it is not evidence that a change improves task accuracy.</p><img src="coordinate_geometry_response.png" alt="Coordinate geometry response"><p><a href="coordinate_geometry_response.pdf">PDF</a> · <a href="coordinate_responses.json">Values and graph checks</a></p></section>'''
        page=page.replace('<section class="card"><h2>Scope and interpretation</h2>',sections+'<section class="card"><h2>Scope and interpretation</h2>')
    (dest/'report.html').write_text(page)
    record=dict(all_numeric_checks_passed=True,source_structure=str(args.structure.resolve()),
        source_artifact_manifest_sha256=sha(args.structure/'artifact_sha256.json'),
        plotted_timepoints=list(STAGES),plotted_branches=list(LABELS),
        same_pair_bins_verified=True,primary_normalization='trace_normalized',
        structure_complete=True,influence_complete=mechanism_done,coordinate_response_complete=mechanism_done,
        source_mechanism_artifact_manifest_sha256=mechanism_provenance,
        influence_primary_metric='Raw squared-Frobenius energy normalized per target; self excluded for spatial maps',
        visual_QA_pending=True,finished_unix=time.time())
    (dest/'render_verification.json').write_text(json.dumps(record,indent=2)+'\n')
    artifact={str(p.relative_to(dest)):sha(p) for p in dest.rglob('*') if p.is_file()}
    (dest/'artifact_sha256.json').write_text(json.dumps(artifact,indent=2)+'\n')
    print(json.dumps(record),flush=True)


if __name__=='__main__':main()
