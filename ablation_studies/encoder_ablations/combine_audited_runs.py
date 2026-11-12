"""Render hash-verified stable27 + D4 extension results; never fit a model.

Primary20D readouts only. No neural/model/OT imports, clustering fits, old writes
or missing-value imputation. Every case has one training seed, not seed SD.
"""
import os
import tempfile
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ.setdefault('MPLCONFIGDIR', os.path.join(tempfile.gettempdir(),'spadot_combined_D4_mpl'))
os.environ.setdefault('XDG_CACHE_HOME', os.path.join(tempfile.gettempdir(),'spadot_combined_D4_cache'))
from pathlib import Path
import argparse
import hashlib
import html
import json
import textwrap
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

CASES = ('ChickenHeart','AxolotlTelencephalon','MouseEmbryo')
PARENT_ARMS = ('dual_native','gat_d1','gat_d2','gat_d3','gat_d6','gat_d3_parammatched',
               'gat_d3_expanded','gat_gat_matched','mlp_gat')
ARMS = ('dual_native','gat_d1','gat_d2','gat_d3','gat_d4','gat_d6',
        'gat_d3_parammatched','gat_d3_expanded','gat_gat_matched','mlp_gat')
LABELS = dict(dual_native='SpaDOT GP+GAT',gat_d1='GAT20 depth1',gat_d2='GAT20 depth2',
    gat_d3='GAT20 depth3',gat_d4='GAT20 depth4 (extension)',gat_d6='GAT20 depth6',
    gat_d3_parammatched='GAT20 d3 matched',gat_d3_expanded='GAT20 d3 expanded',
    gat_gat_matched='GAT10+GAT10 matched',mlp_gat='MLP10+GAT10')
TITLES = dict(ChickenHeart='Chicken heart',AxolotlTelencephalon='Axolotl telencephalon',MouseEmbryo='Mouse embryo')
STAGES = dict(ChickenHeart=['Day4','Day7','Day10','Day14'],
    AxolotlTelencephalon=['St.44','St.54','St.57'],MouseEmbryo=['E12.5','E13.5','E15.5'])
METRICS = ('ARI','NMI','PAS','CHAOS')
SOLVER = '34d2c557045e7c2f0166e647684fa6b711cab26605ff80a6750ed43de9bec805'
COMMON = ('model_adapter.py','standalone_adapter.py','replacement_adapter.py','runtime.py',
          'worker.py','training_solver_guard.py','analysis_runtime.py','evaluate.py','collaborator_metrics.py')
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none'})


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def keep(protected,path,expected=None):
    path = Path(path).resolve(); value = sha(path)
    assert expected is None or value == expected, ('Hash mismatch',str(path))
    assert str(path) not in protected or protected[str(path)] == value
    protected[str(path)] = value
    return value


def safe(root,name):
    name = Path(name)
    assert not name.is_absolute() and '..' not in name.parts
    path = (Path(root)/name).resolve()
    assert Path(root).resolve() in path.parents
    return path


def csv(path):
    try:return pd.read_csv(path,float_precision='round_trip')
    except pd.errors.EmptyDataError:return pd.DataFrame()


def truth(value):
    assert str(value).lower() in ('true','false'), ('Invalid Boolean',value)
    return str(value).lower() == 'true'


def load_audit(audit,manifest_path,kind,protected):
    audit, manifest_path = Path(audit).resolve(), Path(manifest_path).resolve()
    manifest, summary = read(manifest_path), read(audit/'verification.json')
    artifacts = read(audit/'artifact_sha256.json')
    required = {'verification.json','preservation_hashes.json','trial_audit.csv',
                'domain_metrics.csv','regional_metrics.csv','endpoint_long.csv'}
    if kind == 'child':required |= {'parent_audit_provenance.json','parent_reference_verification.json'}
    assert required <= set(artifacts)
    for name,digest in artifacts.items():keep(protected,safe(audit,name),digest)
    keep(protected,audit/'artifact_sha256.json')
    original = read(audit/'preservation_hashes.json')
    keep(protected,manifest_path,original[str(Path(manifest['root'])/'manifest.json')])
    assert summary['manifest_sha256'] == sha(manifest_path)
    assert summary['core_integrity_passed'] and summary['saved_metrics_verified']
    assert summary['all_planned_attempted'] and summary['study_finished'] and not summary['dispatch_aborted']
    assert summary['training_solver_sha256'] == SOLVER and summary['training_epsilon'] == .05
    assert summary['training_returned_growth_index'] == 0
    assert manifest['training_seeds'] == [2026] and manifest['analysis_seed'] == 2026 and manifest['analysis_n_init'] == 10
    assert manifest['executed_epochs'] == list(range(1,100))
    assert manifest['analysis_ot']['epsilon'] == .01
    expected_arms = PARENT_ARMS if kind == 'parent' else ('gat_d4',)
    trials = {t['id']:t for t in manifest['trials']}
    expected = {(case,arm,2026) for case in CASES for arm in expected_arms}
    assert len(trials) == len(manifest['trials']) == len(expected) == summary['planned']
    assert {(t['dataset'],t['arm'],t['seed']) for t in trials.values()} == expected
    for t in trials.values():
        if kind == 'child':assert (t['family'],t['depth'],t['width'],t['graph_expansion']) == ('gat_only',4,512,False)
    tables = {name:csv(audit/(name+'.csv')) for name in ('trial_audit','domain_metrics','regional_metrics','endpoint_long')}
    statuses = tables['trial_audit']
    assert set(statuses.trial_id) == set(trials) and not statuses.trial_id.duplicated().any()
    for row in statuses.to_dict('records'):
        assert row['status'] in ('complete','complete_with_analysis_failure','failed')
        assert truth(row['valid_model']) == (row['status'] != 'failed')
    domains = tables['domain_metrics']
    if domains.empty:domains = pd.DataFrame(columns=['trial_id','dataset','arm','family','seed','representation','role','timepoint',*METRICS])
    primary = domains.loc[domains.role == 'primary'].copy()
    assert not primary.duplicated(['trial_id','timepoint']).any()
    for row in primary.to_dict('records'):
        trial = trials[row['trial_id']]
        assert (row['dataset'],row['arm'],int(row['seed'])) == (trial['dataset'],trial['arm'],2026)
        expected_representation = 'GATonly20' if trial['family'] == 'gat_only' else 'combined20'
        assert row['representation'] == expected_representation and np.isfinite([row[x] for x in METRICS]).all()
        assert truth(statuses.loc[statuses.trial_id==trial['id'],'valid_model'].iloc[0])
    regional = tables['regional_metrics']
    if not regional.empty:regional = regional.loc[regional.role == 'primary'].copy()
    tables.update(primary=primary,regional_primary=regional)
    return dict(manifest=manifest,summary=summary,trials=trials,tables=tables,original=original,
                manifest_path=manifest_path,audit=audit,kind=kind)


def bind_parent(parent,child):
    p,c = parent['manifest'],child['manifest']
    assert c['parent_root'] == p['root'] and c['parent_manifest_sha256'] == sha(parent['manifest_path'])
    assert child['summary']['parent_reference_contract_verified']
    assert child['summary']['parent_controls_counted_as_new_fits'] is False
    assert c['source_sha256'] == p['source_sha256'] and c['package_versions'] == p['package_versions']
    for name in COMMON:assert c['code_sha256'][name] == p['code_sha256'][name]
    assert c['cases'] == p['cases'] and c['protocol']['training_solver'] == p['protocol']['training_solver']
    assert c['protocol']['sampler'] == p['protocol']['sampler'] and c['analysis_ot'] == p['analysis_ot']
    provenance = read(child['audit']/'parent_audit_provenance.json')
    assert provenance['parent_manifest_sha256'] == sha(parent['manifest_path'])
    assert provenance['parent_audit_artifact_manifest_sha256'] == sha(parent['audit']/'artifact_sha256.json')
    references = read(child['audit']/'parent_reference_verification.json')
    assert references['parent_reference_contract_passed'] and references['child_fresh_fits'] == 3
    assert references['reused_comparison_controls'] == 3 and references['no_new_dual_fits']
    for case in CASES:
        t3 = next(t for t in parent['trials'].values() if t['dataset']==case and t['arm']=='gat_d3')
        t4 = next(t for t in child['trials'].values() if t['dataset']==case)
        assert t4['expected_parameters'] - t3['expected_parameters'] == 4200448


def aggregate(parent,child):
    primary = pd.concat([x['tables']['primary'].assign(study_origin=x['kind']) for x in (parent,child)],ignore_index=True)
    regional = pd.concat([x['tables']['regional_primary'].assign(study_origin=x['kind']) for x in (parent,child)],ignore_index=True)
    statuses = pd.concat([x['tables']['trial_audit'].assign(study_origin=x['kind']) for x in (parent,child)],ignore_index=True)
    assert len(statuses) == 30 and not statuses.trial_id.duplicated().any()
    cases = {c['id']:c for c in parent['manifest']['cases']}
    trials = {**parent['trials'],**child['trials']}
    rows = []
    for case in CASES:
        for arm in ARMS:
            trial = next(t for t in trials.values() if t['dataset']==case and t['arm']==arm)
            table = primary.loc[primary.trial_id==trial['id']]
            status = statuses.loc[statuses.trial_id==trial['id']].iloc[0]
            times = set(map(str,cases[case]['timepoints']))
            assert set(table.timepoint.astype(str)) <= times
            complete = set(table.timepoint.astype(str)) == times
            row = dict(trial_id=trial['id'],dataset=case,arm=arm,seed=2026,
                study_origin='child_extension' if arm=='gat_d4' else 'parent',
                status=status.status,parameters=trial['expected_parameters'],
                available_timepoints=len(table),expected_timepoints=len(times),complete=complete,
                training_seeds=1,sample_SD=None)
            for metric in METRICS:row[metric] = float(table[metric].mean()) if complete else np.nan
            rows.append(row)
    means = pd.DataFrame(rows)
    for metric in METRICS:
        ref = means.loc[means.arm=='dual_native'].set_index('dataset')[metric]
        means[metric+'_minus_SpaDOT'] = means.apply(lambda r:r[metric]-ref[r.dataset] if pd.notna(r[metric]) and pd.notna(ref[r.dataset]) else np.nan,axis=1)
    return primary,regional,statuses,means,cases,trials


def save_figure(fig,destination,name,width=2200):
    fig.savefig(destination/(name+'.png'),dpi=width/fig.get_figwidth(),facecolor='white')
    plt.close(fig)


def heatmaps(destination,primary,means,cases):
    fig,axes = plt.subplots(1,3,figsize=(19,8),gridspec_kw={'width_ratios':[5,4,4]})
    for index,(ax,case) in enumerate(zip(axes,CASES)):
        times = list(map(str,cases[case]['timepoints']))
        values = np.full((len(ARMS),len(times)+1),np.nan)
        for i,arm in enumerate(ARMS):
            q = primary.loc[(primary.dataset==case)&(primary.arm==arm)]
            for j,tp in enumerate(times):
                r=q.loc[q.timepoint.astype(str)==tp]
                if not r.empty:assert len(r)==1;values[i,j]=float(r.iloc[0].ARI)
            values[i,-1]=means.loc[(means.dataset==case)&(means.arm==arm),'ARI'].iloc[0]
        cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#e2e8f0')
        image=ax.imshow(np.ma.masked_invalid(values),cmap=cmap,vmin=0,vmax=1,aspect='auto')
        for i in range(len(ARMS)):
            for j in range(len(times)+1):
                v=values[i,j];ax.text(j,i,'—' if np.isnan(v) else f'{v:.3f}',ha='center',va='center',fontsize=9,color='#334155' if np.isnan(v) or v>.58 else 'white')
        ax.set_xticks(range(len(times)+1),[*STAGES[case],'Mean'],fontsize=10,rotation=25,ha='right')
        ax.set_yticks(range(len(ARMS)),[LABELS[a] for a in ARMS] if index==0 else ['']*len(ARMS),fontsize=10)
        ax.set_title(TITLES[case],fontsize=15,pad=12)
    fig.subplots_adjust(left=.18,right=.95,top=.88,bottom=.13,wspace=.12)
    fig.suptitle('Primary20D domain ARI | seed2026 | epoch99',fontsize=19,y=.97)
    fig.text(.55,.035,'D4 is a later extension. Gray = missing/failed, never zero; means require all timepoints. D6 outcomes remain recorded.',ha='center',fontsize=10)
    fig.colorbar(image,ax=list(axes),fraction=.018,pad=.025,label='ARI')
    save_figure(fig,destination,'primary_ARI_by_timepoint')
    fig,axes=plt.subplots(2,2,figsize=(18,14))
    for ax,metric in zip(axes.flat,METRICS):
        values=np.array([[means.loc[(means.dataset==case)&(means.arm==arm),metric].iloc[0] for case in CASES] for arm in ARMS])
        cmap=plt.get_cmap('viridis' if metric in ('ARI','NMI') else 'viridis_r').copy();cmap.set_bad('#e2e8f0')
        finite=values[np.isfinite(values)];low=0.;high=1. if metric in ('ARI','NMI') else float(finite.max()) if len(finite) else 1.
        ax.imshow(np.ma.masked_invalid(values),cmap=cmap,vmin=low,vmax=max(high,1e-12),aspect='auto')
        for i in range(len(ARMS)):
            for j in range(3):
                v=values[i,j];frac=(v-low)/max(high-low,1e-12) if np.isfinite(v) else .5
                light=(frac>.58 if metric in ('ARI','NMI') else frac<.42)
                ax.text(j,i,'—' if np.isnan(v) else f'{v:.3f}',ha='center',va='center',fontsize=10,color='#334155' if np.isnan(v) or light else 'white')
        ax.set_xticks(range(3),[TITLES[x] for x in CASES],rotation=15,fontsize=10)
        ax.set_yticks(range(len(ARMS)),[LABELS[a] for a in ARMS],fontsize=10)
        ax.set_title(metric+(' (higher better)' if metric in ('ARI','NMI') else ' (lower better)'),fontsize=14)
    fig.subplots_adjust(left=.15,right=.98,bottom=.06,top=.93,hspace=.28,wspace=.5)
    fig.suptitle('Equal-timepoint means | one training seed; no seed SD',fontsize=19,y=.98)
    save_figure(fig,destination,'primary_metric_means')


def display(value):return ' '.join(str(value).split())


def mapping(labels,groups,palette):
    clusters,annotations=sorted(set(groups)),sorted(set(labels))
    counts=np.array([[np.sum((groups==g)&(labels==a)) for a in annotations] for g in clusters])
    ri,ci=linear_sum_assignment(-counts); result={}
    for i,j in zip(ri,ci):
        if counts[i,j]>0:result[clusters[i]]=dict(annotation=annotations[j],color=palette[display(annotations[j])],overlap=int(counts[i,j]))
    unmatched=[g for g in clusters if g not in result]
    for i,g in enumerate(unmatched):
        gray=int(115+90*i/max(1,len(unmatched)-1));color=f'#{gray:02x}{gray:02x}{gray:02x}'
        result[g]=dict(annotation=None,color=color,overlap=0)
    return result


def load_spatial(directory,audited,protected):
    directory=Path(directory).resolve();bundle=read(directory/'bundle.json')
    index=read(directory/'artifact_sha256.json')
    for name,digest in index.items():keep(protected,safe(directory,name),digest)
    keep(protected,directory/'artifact_sha256.json')
    assert 'bundle.json' in index
    assert bundle['manifest_sha256']==sha(audited['manifest_path'])
    assert bundle['audit_verification_sha256']==sha(audited['audit']/'verification.json')
    assert bundle['audit_artifact_manifest_sha256']==sha(audited['audit']/'artifact_sha256.json')
    assert bundle['original_input_metadata_only'] and bundle['current_native_study_predictions_only']
    assert bundle['no_historical_predictions_or_spatial_cache']
    assert set(bundle['cases'])==set(CASES) and set(bundle['trials'])<=set(audited['trials'])
    meta,predictions={},{}
    for case in audited['manifest']['cases']:
        record=bundle['cases'][case['id']]
        assert record['input_sha256']==case['input_sha256']
        path=safe(directory,record['file']);keep(protected,path,record['sha256'])
        with np.load(path,allow_pickle=False) as archive:values={k:np.array(archive[k],copy=True) for k in archive.files}
        assert values['coordinates'].shape==(case['shape'][0],2) and np.isfinite(values['coordinates']).all()
        assert len(set(values['spot']))==len(values['spot'])==record['n']==case['shape'][0]
        assert np.array_equal(values['annotation'],np.array([x.replace('\n',' ') for x in values['annotation_source']]))
        assert set(values['timepoint'])==set(map(str,case['timepoints']))
        meta[case['id']]=values
    for ident,record in bundle['trials'].items():
        trial=audited['trials'][ident];case=trial['dataset'];values=meta[case]
        assert record['metadata_sha256']==bundle['cases'][case]['sha256']
        logical=str(Path(audited['manifest']['root'])/'results'/ident/'predictions.csv')
        assert record['source_prediction_sha256']==audited['original'][logical]
        path=safe(directory,record['file']);keep(protected,path,record['sha256'])
        with np.load(path,allow_pickle=False) as archive:groups=np.array(archive['cluster'],dtype=str,copy=True)
        assert len(groups)==len(values['spot'])==record['n']
        rows=audited['tables']['primary'].loc[audited['tables']['primary'].trial_id==ident]
        assert not rows.empty
        for row in rows.to_dict('records'):
            mask=values['timepoint']==str(row['timepoint'])
            assert np.isclose(adjusted_rand_score(values['annotation'][mask],groups[mask]),row['ARI'],atol=1e-12,rtol=1e-11)
        predictions[ident]=groups
    expected=set(audited['tables']['primary'].trial_id)
    assert set(predictions)==expected, 'An audited primary domain output is missing from spatial export'
    return meta,predictions,bundle


def spatial(destination,parent,child,parent_path,child_path,palette_path,protected):
    pm,pp,pb=load_spatial(parent_path,parent,protected)
    cm,cp,cb=load_spatial(child_path,child,protected)
    palette_source=read(palette_path);keep(protected,palette_path)
    palettes={case:{display(k):v for k,v in palette_source['palette_binding'][case]['annotation_colors'].items()} for case in CASES}
    checks,mappings,figures=[],[],[]
    for case in CASES:
        values=pm[case]
        for key in values:assert np.array_equal(values[key],cm[case][key]), ('Parent/child exported metadata differ',case,key)
        palette=palettes[case]
        assert {display(x) for x in values['annotation']}<=set(palette)
        times=list(map(str,next(c for c in parent['manifest']['cases'] if c['id']==case)['timepoints']))
        height,width=10.6,4*len(times)+.8
        fig,axes=plt.subplots(3,len(times),figsize=(width,height),squeeze=False)
        fig.subplots_adjust(left=.09,right=.99,top=.86,bottom=.17,wspace=.06,hspace=.20)
        fig.suptitle(TITLES[case]+' | standard SpaDOT and D4 extension',fontsize=19,y=.98)
        fig.text(.54,.935,'Fixed solver · original native two-hop sampler · epoch99 · seed2026 · final KMeans n_init=10',ha='center',fontsize=10)
        for j,tp in enumerate(times):axes[0,j].set_title(STAGES[case][j],fontsize=14,pad=10)
        for i,label in enumerate(('Original\nannotation','SpaDOT\nGP+GAT','GAT-only20\ndepth4')):
            pos=axes[i,0].get_position();fig.text(.035,(pos.y0+pos.y1)/2,label,ha='center',va='center',rotation=90,fontsize=13,fontweight='bold')
        for j,tp in enumerate(times):
            mask=values['timepoint']==tp;xy=values['coordinates'][mask];labels=values['annotation'][mask]
            lower,upper=xy.min(0),xy.max(0);span=np.maximum(upper-lower,1e-12)
            for i,arm in enumerate((None,'dual_native','gat_d4')):
                ax=axes[i,j];subtitle=f'{len(labels):,} spots';colors=[palette[display(x)] for x in labels]
                if arm:
                    origin=parent if arm=='dual_native' else child
                    lookup=pp if arm=='dual_native' else cp;ident=case+'__'+arm+'__seed2026'
                    rows=origin['tables']['primary'].loc[(origin['tables']['primary'].trial_id==ident)&(origin['tables']['primary'].timepoint.astype(str)==tp)]
                    if ident not in lookup or rows.empty:
                        status=origin['tables']['trial_audit'].loc[origin['tables']['trial_audit'].trial_id==ident,'status'].iloc[0]
                        colors=['#cbd5e1']*len(labels);subtitle='Training failure\nNo final domains' if status=='failed' else 'Missing domain endpoint\nNo imputation'
                        checks.append(dict(dataset=case,arm=arm,timepoint=tp,available=False,status=status))
                    else:
                        assert len(rows)==1
                        groups=lookup[ident][mask];score=float(adjusted_rand_score(labels,groups))
                        assert np.isclose(score,rows.iloc[0].ARI,atol=1e-12,rtol=1e-11)
                        bind=mapping(labels,groups,palette);colors=[bind[x]['color'] for x in groups];subtitle=f'ARI {score:.3f}'
                        mappings.append(dict(dataset=case,arm=arm,timepoint=tp,mapping=bind))
                        checks.append(dict(dataset=case,arm=arm,timepoint=tp,available=True,ARI=score,independent_ARI_matches=True))
                size=8 if case=='ChickenHeart' else 3.2 if case=='AxolotlTelencephalon' else 5.5
                ax.scatter(xy[:,0],xy[:,1],c=colors,s=size,linewidths=0)
                ax.set_aspect('equal');ax.set_xlim(lower[0]-.035*span[0],upper[0]+.035*span[0]);ax.set_ylim(lower[1]-.035*span[1],upper[1]+.035*span[1]);ax.axis('off')
                ax.text(.5,-.025,subtitle,transform=ax.transAxes,ha='center',va='top',fontsize=11,color='#334155')
        legend=[Patch(facecolor=color,label=textwrap.fill(label,24)) for label,color in palette.items() if label in {display(x) for x in values['annotation']}]
        legend.append(Patch(facecolor='#8c96a5',label='Unmatched / unavailable'))
        fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.54,.035),ncol=5,fontsize=8.5,frameon=False)
        fig.text(.54,.012,'Maximum-overlap color matching uses annotations for display only; predictions and ARIs are unchanged.',ha='center',fontsize=9,color='#475569')
        name=case+'_SpaDOT_D4_spatial';save_figure(fig,destination,name)
        figures.append(dict(dataset=case,file=name+'.png'))
    write(destination/'spatial_display_mappings.json',mappings)
    write(destination/'spatial_panel_verification.json',checks)
    return figures,checks


def table(frame,columns):
    return frame.loc[:,[x for x in columns if x in frame]].to_html(index=False,na_rep='—',float_format=lambda x:f'{x:.4f}',escape=True)


def render(args):
    protected={};keep(protected,__file__)
    parent=load_audit(args.parent_audit,args.parent_manifest,'parent',protected)
    child=load_audit(args.child_audit,args.child_manifest,'child',protected)
    bind_parent(parent,child)
    destination=args.destination.resolve()
    assert not destination.exists()
    for source in (args.parent_audit,args.child_audit,args.parent_spatial,args.child_spatial):
        assert Path(source).resolve() not in destination.parents and destination!=Path(source).resolve()
    destination.mkdir(parents=True,exist_ok=False)
    primary,regional,statuses,means,cases,trials=aggregate(parent,child)
    for name,frame in [('primary_per_timepoint',primary),('primary_equal_timepoint_means',means),('primary_regional_metrics',regional),('all_trial_outcomes',statuses)]:frame.to_csv(destination/(name+'.csv'),index=False)
    heatmaps(destination,primary,means,cases)
    figures,checks=spatial(destination,parent,child,args.parent_spatial,args.child_spatial,args.palette_file,protected)
    failure=statuses.loc[statuses.status=='failed']
    comparison=means.loc[means.arm!='dual_native',['dataset','arm',*[x+'_minus_SpaDOT' for x in METRICS]]]
    comparison.to_csv(destination/'paired_primary_mean_differences.csv',index=False)
    counts={x:sum(statuses.status==x) for x in ('complete','complete_with_analysis_failure','failed')}
    assert sum(counts.values())==30
    images=''.join(f'<h3>{html.escape(TITLES[x["dataset"]])}</h3><a href="{x["file"]}"><img src="{x["file"]}" alt="Annotation, SpaDOT and D4 spatial comparison"></a>' for x in figures)
    parent_link=f'<p><a href="{html.escape(str(Path(args.parent_report).resolve()))}">Full parent nine-arm report and spatial panels</a></p>' if args.parent_report else ''
    if args.parent_report:assert Path(args.parent_report).is_file();keep(protected,args.parent_report)
    region_columns=['dataset','arm','timepoint','pair','pair_ARI','dominant_clusters_different','cross_region_same_cluster_fraction',
        'Heart_precision','Heart_recall','Heart_contamination','Liver_precision','Liver_recall','Liver_contamination',
        'MP_precision','MP_recall','MP_contamination','DP_precision','DP_recall','DP_contamination','VZ_precision','VZ_recall','VZ_contamination']
    body=f'''<!doctype html><html><head><meta charset="utf-8"><title>Fixed-solver native study plus D4 extension</title><style>
body{{font:15px/1.5 system-ui,sans-serif;color:#163047;margin:32px auto;max-width:1500px;padding:0 20px}}h1,h2{{line-height:1.2}}img{{width:100%;height:auto}}table{{border-collapse:collapse;font-size:12px;white-space:nowrap}}th,td{{padding:7px 9px;border-bottom:1px solid #dce4eb;text-align:right}}th{{background:#eff4f8}}.scroll{{overflow-x:auto}}.note{{background:#edf4fa;padding:16px;border-left:4px solid #31577f}}a{{color:#235e8e}}</style></head><body>
<h1>Native-sampler encoder comparison with fixed OT</h1><p>Three datasets · seed2026 only · executed epoch99 · original native two-hop finite sampler · final KMeans n_init10.</p>
<div class="note"><b>30 fresh attempts:</b> parent27 plus three later D4 fits. {int(counts['complete'])} complete, {int(counts['complete_with_analysis_failure'])} valid fitted models with analysis failure, {int(counts['failed'])} training failures. Failed values remain missing. D4 was added after reviewing D6 failures; it does not replace or retry those recorded outcomes. No seed SD, confidence interval or generalization estimate is available.</div>
<h2>Primary20D performance</h2><p>All metrics use the whole20D representation. Equal-timepoint means require every stage; no partial means or failure-as-zero scoring. Higher ARI/NMI and lower PAS/CHAOS are preferable. PAS/CHAOS measure continuity, not proof that one encoder is exclusively local or global.</p>
<img src="primary_ARI_by_timepoint.png" alt="Primary ARI by stage and arm"><img src="primary_metric_means.png" alt="Primary ARI NMI PAS CHAOS means">
<details open><summary>Full means and coverage</summary><div class="scroll">{table(means,['dataset','arm','status','parameters','available_timepoints','expected_timepoints',*METRICS])}</div></details>
<details><summary>Per-timepoint ARI/NMI/PAS/CHAOS</summary><div class="scroll">{table(primary,['dataset','arm','timepoint','n',*METRICS])}</div></details>
<h2>Regional separation and purity</h2><p>Heart/Liver and MP/DP restricted ARI does not establish purity against other tissues. Precision, recall, cross-region mixing and contamination use the original annotations and all regions. VZ is reported at every Axolotl stage. No anatomy-selected epoch or model was introduced.</p>
<div class="scroll">{table(regional,region_columns)}</div>
<h2>SpaDOT and D4 spatial patterns</h2><p>Original input coordinates, unchanged point order and independently audited predictions. Hungarian one-to-one maximum-overlap annotation colors assist display only. Failed Mouse D4 panels show gray original coordinates and explicitly contain no predicted domains.</p>{images}{parent_link}
<h2>Recorded training failures</h2><div class="scroll">{table(failure,['dataset','arm','status','error_category','last_recorded_loss_epoch','last_finite_prefix_epoch','error'])}</div>
<h2>What this comparison can establish</h2><p>Compare architectures case by case and retain improvements as well as counterexamples. Adding a GP branch is supported only where the matched empirical comparison improves relevant outcomes; a difference in branch patterns or model capacity does not prove universal necessity. Standalone GAT20 arms are actual separately fitted models; secondary10D readouts from the two-branch models are not standalone ablations.</p>
<p>D4 has4,200,448more GAT parameters than natural D3. It is not capacity matched. Four graph layers still use native two-hop sampled unions, whereas final inference uses the full graph. The expanded arm changes graph cutoff/cap and can increase native fanout. All arms share the fixed NumPy solver: training epsilon0.05 with the original first-growth return; spot analysis epsilon0.01 with WOT's third map. The fix does not guarantee neural or OT numerical stability.</p>
<p>Primary comparisons retain original saved endpoints. Posthoc-derived OT scores without original pipeline score records are excluded. Current audit integrity does not imply all endpoints exist.</p>
<p>Downloads: <a href="primary_per_timepoint.csv">per-stage scores</a> · <a href="primary_equal_timepoint_means.csv">means/coverage</a> · <a href="primary_regional_metrics.csv">regional metrics</a> · <a href="paired_primary_mean_differences.csv">same-seed differences</a> · <a href="all_trial_outcomes.csv">all outcomes</a> · <a href="verification.json">verification</a>.</p></body></html>'''
    (destination/'report.html').write_text(body)
    assert all(sha(path)==digest for path,digest in protected.items()), 'Protected render inputs changed'
    write(destination/'source_sha256.json',protected)
    write(destination/'verification.json',dict(all_passed=True,planned_fresh_attempts=30,parent_fresh=27,child_fresh=3,
        counts={k:int(v) for k,v in counts.items()},training_seeds=[2026],sample_SD=None,
        primary20D_only=True,failed_scores_not_zero=True,no_partial_timepoint_means=True,
        D6_preserved=True,D4_later_extension=True,additional_D4_parameters_over_D3=4200448,
        parent_manifest_sha256=sha(args.parent_manifest),child_manifest_sha256=sha(args.child_manifest),
        parent_audit_artifact_sha256=sha(args.parent_audit/'artifact_sha256.json'),
        child_audit_artifact_sha256=sha(args.child_audit/'artifact_sha256.json'),
        parent_reference_provenance_bound=True,spatial_panels=len(checks),
        available_spatial_panel_ARIs_verified=sum(x['available'] for x in checks),
        annotation_assisted_display_only=True,no_model_or_clustering_fit=True,no_neural_forward=True,no_OT_solve=True,
        input_files_verified=len(protected),visual_inspection_performed=False))
    write(destination/'artifact_sha256.json',{str(p.relative_to(destination)):sha(p) for p in sorted(destination.rglob('*')) if p.is_file() and p.name!='artifact_sha256.json'})
    print(json.dumps(dict(report=str(destination/'report.html'),all_passed=True,counts={k:int(v) for k,v in counts.items()})))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('parent-audit','child-audit','parent-manifest','child-manifest','parent-spatial','child-spatial','palette-file','destination'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--parent-report',type=Path)
    render(parser.parse_args())


if __name__=='__main__':main()
