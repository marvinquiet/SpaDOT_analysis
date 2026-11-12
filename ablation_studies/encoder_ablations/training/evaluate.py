"""Evaluate dual-branch replacements; preserve domains before scoped epsilon.01 OT."""

from contextlib import contextmanager
import hashlib
import inspect
import json
from pathlib import Path
import time
import traceback
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics.pairwise import pairwise_distances
from threadpoolctl import threadpool_limits

from SpaDOT.utils import _analyze_utils as A
from analysis_runtime import installed_solver, verify_saved_maps
import collaborator_metrics as M


def _dense(value):
    return value.toarray() if sparse.issparse(value) else np.asarray(value)


def _write(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def _raw_table(transport, source, target, rows, columns):
    source_membership = (source[:, None] == np.asarray(rows)[None, :]).astype(float)
    target_membership = (target[:, None] == np.asarray(columns)[None, :]).astype(float)
    values = source_membership.T @ (transport @ target_membership)
    return pd.DataFrame(values, index=rows, columns=columns)


@contextmanager
def _serial_wot_cost():
    """Preserve WOT cost arithmetic and solver; prevent joblib subprocesses."""
    model_class = A.wot.ot.OTModel
    descriptor = inspect.getattr_static(model_class, 'compute_default_cost_matrix')
    if not isinstance(descriptor, staticmethod):
        raise TypeError('Expected WOT compute_default_cost_matrix to be static')
    if list(inspect.signature(descriptor.__func__).parameters) != ['a', 'b', 'eigenvals']:
        raise TypeError('Unexpected WOT cost API')
    calls = []
    source = inspect.getsource(descriptor.__func__)
    record = {'arithmetic': 'sqeuclidean / median', 'n_jobs': 1,
              'solver_changed': False, 'calls': calls,
              'original_cost_source_sha256': hashlib.sha256(source.encode()).hexdigest()}

    def serial_cost(a, b, eigenvals=None):
        before_shape = [list(a.shape), list(b.shape)]
        if eigenvals is not None:
            a = a.dot(eigenvals)
            b = b.dot(eigenvals)
        matrix = pairwise_distances(_dense(a), _dense(b), metric='sqeuclidean', n_jobs=1)
        median = float(np.median(matrix))
        if not np.isfinite(matrix).all() or median <= 0 or not np.isfinite(median):
            raise ValueError('Nonfinite or zero-median WOT cost')
        calls.append({'input_shapes': before_shape, 'cost_shape': list(matrix.shape),
                      'median_cost': median, 'eigenvalue_scaling': eigenvals is not None})
        return matrix / median

    model_class.compute_default_cost_matrix = staticmethod(serial_cost)
    try:
        yield record
    finally:
        model_class.compute_default_cost_matrix = descriptor


def _metadata(case, trial):
    return {'dataset': case['id'], 'kind': case['kind'], 'trial_id': trial['id'],
            'arm': trial['arm'], 'depth': int(trial['depth']),
            'width': int(trial['width']), 'seed': int(trial['seed']),
            'family': trial['family'], 'kernel_scale': 0.1 if trial['family'] == 'dual' else None, 'clustering_seed': 2026, 'clustering_n_init': 10}


def _domain_scores(prediction, case, metadata, representation):
    rows = []
    for i, timepoint in enumerate(case['timepoints']):
        timepoint = str(timepoint)
        mask = prediction.obs['timepoint'].to_numpy() == timepoint
        observations = prediction.obs.loc[mask]
        coordinates = np.asarray(prediction.obsm['spatial'])[mask]
        if observations.empty or observations['kmeans'].nunique() != case['clusters'][i]:
            raise ValueError('Incorrect number of nonempty domains at ' + timepoint)
        annotation = observations['annotation'].to_numpy()
        clusters = observations['kmeans'].to_numpy()
        row = dict(metadata, representation=representation, timepoint=timepoint,
                   n=len(observations), n_clusters=int(case['clusters'][i]),
                   ARI=float(adjusted_rand_score(annotation, clusters)),
                   NMI=float(normalized_mutual_info_score(annotation, clusters)),
                   PAS=float(M.compute_PAS(coordinates, clusters)),
                   CHAOS=float(M.compute_CHAOS(coordinates, clusters)),
                   PAS_annotation=float(M.compute_PAS(coordinates, annotation)),
                   CHAOS_annotation=float(M.compute_CHAOS(coordinates, annotation)))
        if not np.isfinite([row[k] for k in ('ARI', 'NMI', 'PAS', 'CHAOS',
                                             'PAS_annotation', 'CHAOS_annotation')]).all():
            raise ValueError('Nonfinite domain metric at ' + timepoint)
        rows.append(row)
    return rows


def _save_predictions(prediction, path):
    pd.DataFrame({'spot': prediction.obs_names,
                  'timepoint': prediction.obs['timepoint'].to_numpy(),
                  'cluster': prediction.obs['kmeans'].to_numpy(),
                  'annotation': prediction.obs['annotation'].to_numpy()}).to_csv(path, index=False)


def _regional_scores(prediction, case, metadata, representation):
    """Prespecified posthoc biological checks; annotations never choose a fit."""
    results = []
    timepoints = [str(t) for t in case['timepoints']]
    pairs = ([(tp, 'Heart', 'Liver') for tp in timepoints[:2]]
             if case['id'] == 'MouseEmbryo' else
             [('t2', 'MP', 'DP')] if case['id'] == 'AxolotlTelencephalon' else [])
    for tp, left, right in pairs:
        q = prediction.obs[prediction.obs['timepoint'] == tp]
        labels, clusters = q['annotation'].to_numpy(), q['kmeans'].to_numpy()
        if not (labels == left).any() or not (labels == right).any():
            raise ValueError('Missing prespecified regional label: ' + tp + ':' + left + '/' + right)
        mask = np.isin(labels, [left, right])
        row = dict(metadata, representation=representation, timepoint=tp, pair=left + '/' + right,
                   pair_ARI=float(adjusted_rand_score(labels[mask], clusters[mask])))
        dominants = {}
        for tissue in [left, right]:
            truth = labels == tissue
            ks, ns = np.unique(clusters[truth], return_counts=True)
            k = ks[np.argmax(ns)]
            dominants[tissue] = str(k)
            selected = clusters == k
            hits = int(np.sum(selected & truth))
            row.update({tissue + '_dominant_cluster': str(k),
                        tissue + '_dominant_tie_count': int(np.sum(ns == ns.max())),
                        tissue + '_precision': hits / int(selected.sum()),
                        tissue + '_recall': hits / int(truth.sum()),
                        tissue + '_contamination': 1 - hits / int(selected.sum())})
        row['dominant_clusters_different'] = dominants[left] != dominants[right]
        row['cross_region_same_cluster_fraction'] = float(np.mean(
            clusters[labels == left][:, None] == clusters[labels == right][None, :]))
        results.append(row)
    if case['id'] == 'AxolotlTelencephalon':
        for tp in timepoints:
            q = prediction.obs[prediction.obs['timepoint'] == tp]
            labels, clusters = q['annotation'].to_numpy(), q['kmeans'].to_numpy()
            truth = labels == 'VZ'
            if not truth.any():
                raise ValueError('Missing prespecified VZ label at ' + tp)
            ks, ns = np.unique(clusters[truth], return_counts=True)
            k = ks[np.argmax(ns)]
            selected = clusters == k
            hits = int(np.sum(selected & truth))
            results.append(dict(metadata, representation=representation, timepoint=tp, pair='VZ',
                                VZ_dominant_cluster=str(k), VZ_precision=hits / int(selected.sum()),
                                VZ_recall=hits / int(truth.sum()), VZ_contamination=1 - hits / int(selected.sum())))
    return results


def evaluate(root, out, case, trial, latent):
    """Domain metrics are retained even if subsequent association analysis fails.

    The caller installs the frozen source and collaborator_metrics on sys.path.
    This function evaluates each latent once and never retries an OT failure.
    """
    root, out = Path(root), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    start = out / 'evaluation_started.json'
    if start.exists():
        raise FileExistsError('Evaluation already attempted: ' + str(start))
    metadata = _metadata(case, trial)
    with start.open('x') as handle:
        json.dump(dict(metadata, root=str(root), started_unix=time.time()), handle, indent=2)
        handle.write('\n')
    lat = ad.read_h5ad(latent) if isinstance(latent, (str, Path)) else latent.copy()
    if lat.n_vars != 20 or not np.isfinite(_dense(lat.X)).all():
        raise ValueError('Expected a finite 20-dimensional latent')
    if not lat.obs_names.is_unique:
        raise ValueError('Spot identifiers must be unique')
    for column in ('timepoint', 'annotation'):
        if lat.obs[column].isna().any():
            raise ValueError('Missing ' + column)
        lat.obs[column] = lat.obs[column].astype(str)
    lat.obs['annotation'] = lat.obs['annotation'].map(M.norm)
    timepoints = [str(x) for x in case['timepoints']]
    if sorted(lat.obs['timepoint'].unique()) != timepoints:
        raise ValueError('Case timepoints must match the sorted latent timepoints')
    if len(case['clusters']) != len(timepoints):
        raise ValueError('One final cluster count is required per timepoint')

    result = dict(metadata, domain_passed=False, domain_plots_passed=False,
                  association_passed=False, domain_metrics=[], combined_domain_metrics=[],
                  transition_metrics=[], aggregation_audit=[], domain_failures=[])
    combined = None
    primary = 'GATonly20' if trial['family'] == 'gat_only' else 'combined20'
    representations = ([(primary, slice(None))] if trial['family'] == 'gat_only' else
                       [(primary, slice(None)), (('GPVAE10' if trial['family']=='dual' else 'MLP10' if trial['family']=='mlp_gat' else 'GAT_A10'), slice(0, 10)), ('GAT10', slice(10, 20))])
    result['primary_representation'] = primary
    result['regional_metrics'] = []
    for representation, dimensions in representations:
        try:
            with threadpool_limits(limits=1):
                predicted = A.KMeans_Clustering(lat[:, dimensions].copy(), case['clusters'],
                                                seed=2026, n_init=10)
                predicted.obs['kmeans'] = predicted.obs['kmeans'].astype(str)
                rows = _domain_scores(predicted, case, metadata, representation)
            suffix = '' if representation == primary else '_' + representation
            _save_predictions(predicted, out / ('predictions' + suffix + '.csv'))
            pd.DataFrame(rows).to_csv(out / ('domain_metrics' + suffix + '.csv'), index=False)
            result['domain_metrics'].extend(rows)
            pd.DataFrame(result['domain_metrics']).to_csv(
                out / 'domain_metrics_all_representations.csv', index=False)
            result['regional_metrics'].extend(_regional_scores(predicted, case, metadata, representation))
            if result['regional_metrics']:
                pd.DataFrame(result['regional_metrics']).to_csv(out / 'regional_metrics.csv', index=False)
            if representation == primary:
                combined = predicted
                result['combined_domain_metrics'] = rows
        except Exception as error:
            failure = {'representation': representation, 'error': str(error),
                       'traceback': traceback.format_exc()}
            result['domain_failures'].append(failure)
            _write(out / ('domain_failure_' + representation + '.json'), failure)
    result['domain_passed'] = len(result['domain_metrics']) == len(representations) * len(timepoints)
    _write(out / 'domain_evaluation.json', {
        'all_representations_passed': result['domain_passed'],
        'metrics': result['domain_metrics'], 'failures': result['domain_failures']})

    if combined is None:
        failure = {'stage': 'association', 'error': 'Primary 20D domain evaluation did not complete',
                   'domain_metrics_preserved': True, 'automatic_retry': False}
        _write(out / 'analysis_failure.json', failure)
        result['analysis_failure'] = failure
        _write(out / 'evaluation.json', result)
        return result
    args = SimpleNamespace(output_dir=str(out), prefix='')
    try:
        A.plot_domains(args, combined)
        result['domain_plots_passed'] = True
    except Exception as error:
        failure = {'error': str(error), 'traceback': traceback.format_exc()}
        _write(out / 'domain_plot_failure.json', failure)
        result['domain_plot_failure'] = failure

    try:
        graph = M.transition_graphs(case['id'], case['graph'],
                                    sorted(combined.obs['annotation'].unique()))
        if len(graph) != len(timepoints) - 1:
            raise ValueError('One correspondence graph is required per adjacent pair')
        with threadpool_limits(limits=1), _serial_wot_cost() as cost_record:
            with installed_solver(A.wot.ot.OTModel, root, out) as solver_record:
                A.OT_analysis(args, combined)
            verify_saved_maps(solver_record, out, ad.read_h5ad, len(timepoints)-1)
        result['analysis_epsilon'] = .01
        result['analysis_solver_patch_sha256'] = solver_record['patch_file_sha256']
        result['wot_cost'] = cost_record
        _write(out / 'wot_cost_execution.json', cost_record)
        if len(cost_record['calls']) != len(timepoints) - 1:
            raise ValueError('Unexpected number of WOT cost calculations')
        if any(call['eigenvalue_scaling'] or any(shape[1] != 20 for shape in call['input_shapes'])
               for call in cost_record['calls']):
            raise ValueError('Expected unprojected 20-dimensional association costs')

        for i in range(len(timepoints) - 1):
            coupling = ad.read_h5ad(out / f'OT_{i}_{i + 1}.h5ad')
            transport = _dense(coupling.X)
            if (not np.isfinite(transport).all() or (transport < 0).any()
                    or float(transport.sum()) <= 0):
                raise ValueError('Invalid spot transport matrix')
            source = combined.obs.loc[coupling.obs_names]
            target = combined.obs.loc[coupling.var_names]
            if (set(coupling.obs_names) != set(combined.obs_names[
                    combined.obs['timepoint'] == timepoints[i]])
                    or set(coupling.var_names) != set(combined.obs_names[
                        combined.obs['timepoint'] == timepoints[i + 1]])):
                raise ValueError('OT spot/timepoint alignment mismatch')
            source_labels = source['annotation'].to_numpy()
            target_labels = target['annotation'].to_numpy()
            rows, columns = sorted(set(source_labels)), sorted(set(target_labels))
            annotation_mass = _raw_table(transport, source_labels, target_labels, rows, columns)
            direct = np.asarray([[transport[np.ix_(source_labels == row,
                                                   target_labels == column)].sum()
                                  for column in columns] for row in rows])
            annotation_error = float(np.abs(direct - annotation_mass.to_numpy()).max())
            if not np.allclose(direct, annotation_mass.to_numpy(), atol=1e-12, rtol=1e-10):
                raise ValueError('Direct annotation aggregation verification failed')
            if not np.isclose(annotation_mass.to_numpy().sum(), transport.sum(),
                              atol=1e-12, rtol=1e-10):
                raise ValueError('Annotation aggregation did not preserve total mass')

            table = ad.read_h5ad(out / f'transition_table_{i}_{i + 1}.h5ad')
            cluster_mass = _raw_table(transport, source['SpaDOT_pred_labels'].to_numpy(),
                                      target['SpaDOT_pred_labels'].to_numpy(),
                                      list(table.obs_names), list(table.var_names))
            original = _dense(table.X)
            if not np.isfinite(original).all() or original.sum() <= 0:
                raise ValueError('Invalid saved predicted-domain transport table')
            difference = float(np.abs(cluster_mass.to_numpy() / cluster_mass.to_numpy().sum()
                                      - original / original.sum()).max())
            if difference >= 1e-10:
                raise ValueError('Direct predicted-domain aggregation verification failed')
            majority = [frame.groupby('SpaDOT_pred_labels')['annotation'].agg(
                lambda values: values.value_counts().idxmax()).to_dict()
                for frame in (source, target)]
            mapped = cluster_mass.copy()
            mapped.index = [majority[0][label] for label in mapped.index]
            mapped.columns = [majority[1][label] for label in mapped.columns]
            mapped = mapped.groupby(level=0).sum().T.groupby(level=0).sum().T
            for table_kind, mass in [('WOT_annotation', annotation_mass),
                                     ('SpaDOT_domains_majority', mapped),
                                     ('SpaDOT_domains_raw', cluster_mass)]:
                mass.to_csv(out / f'mass_{table_kind}_{i}_{i + 1}.csv')
                probability = M.min_norm(mass)
                if not np.isfinite(probability.to_numpy()).all():
                    raise ValueError('Nonfinite conditional transition probability')
                probability.to_csv(out / f'transition_{table_kind}_{i}_{i + 1}.csv')
                if table_kind != 'SpaDOT_domains_raw':
                    curated, weight = M.evaluate_using_curated_transitions(probability, graph[i])
                    naive = float(M.accuracy(probability, graph[i]))
                    if not np.isfinite([naive, curated, weight]).all():
                        raise ValueError('Nonfinite transition accuracy')
                    result['transition_metrics'].append(dict(
                        metadata, i=i, transition=timepoints[i] + '->' + timepoints[i + 1],
                        table=table_kind, naive_accuracy=naive,
                        curated_accuracy=float(curated), total_weight=float(weight)))
                    pd.DataFrame(result['transition_metrics']).to_csv(
                        out / 'transition_metrics.csv', index=False)
            result['aggregation_audit'].append({
                'pair': i, 'spot_shape': list(transport.shape),
                'annotation_aggregation_max_error': annotation_error,
                'cluster_aggregation_max_error': difference,
                'annotation_total_mass_preserved': True,
                'graph_source_labels_without_outgoing_edges': [label for label in rows
                    if label not in graph[i] or graph[i].out_degree(label) == 0],
                'graph_target_labels_absent_from_data': sorted({str(b) for _, b in graph[i].edges
                                                               if b not in columns})})
            _write(out / 'aggregation_audit.json', result['aggregation_audit'])
        A.plot_OT(args, combined)
        result['association_passed'] = True
    except Exception as error:
        failure = {'stage': 'association', 'error': str(error),
                   'traceback': traceback.format_exc(), 'domain_metrics_preserved': True,
                   'automatic_retry': False, 'completed_pairs': len(result['aggregation_audit'])}
        if 'cost_record' in locals():
            result['wot_cost'] = cost_record
            _write(out / 'wot_cost_execution.json', cost_record)
        _write(out / 'analysis_failure.json', failure)
        result['analysis_failure'] = failure
    result['evaluation_finished_unix'] = time.time()
    _write(out / 'evaluation.json', result)
    return result
