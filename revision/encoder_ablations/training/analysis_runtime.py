"""Scoped CPU-only OT analysis solver replacement; frozen WOT growth is retained."""
from contextlib import contextmanager
from pathlib import Path
import hashlib
import importlib.util
import inspect
import json
import sys
import time

import numpy as np

PATCH_SHA256 = '34d2c557045e7c2f0166e647684fa6b711cab26605ff80a6750ed43de9bec805'
EXPECTED = dict(epsilon=.01, epsilon0=1., lambda1=.1, lambda2=5., tau=1000., growth_iters=3)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def array_sha(value):
    x = np.ascontiguousarray(value)
    h = hashlib.sha256()
    h.update(str(x.dtype).encode())
    h.update(json.dumps(list(x.shape)).encode())
    h.update(memoryview(x).cast('B'))
    return h.hexdigest()


def load_patch(patch_root):
    path = Path(patch_root) / 'model_patch_source/SpaDOT/utils/OT_loss/ot_solvers.py'
    assert sha(path) == PATCH_SHA256, 'Unexpected numerical patch source'
    spec = importlib.util.spec_from_file_location('_spadot_analysis_stable_ot_patch', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _pair_frame(instance):
    frame = sys._getframe(1)
    while frame is not None:
        if frame.f_code.co_name == 'compute_single_transport_map' and frame.f_locals.get('self') is instance:
            return float(frame.f_locals['t0']), float(frame.f_locals['t1'])
        frame = frame.f_back
    raise RuntimeError('Solver called outside the original WOT transport-map entrypoint')


@contextmanager
def installed_solver(model_class, patch_root, destination):
    """Replace each newly constructed instance's solver, restoring __init__ on exit.

    No WOT growth-loop or cost-method replacement occurs here. The frozen
    evaluator remains responsible for its own serial cost context.
    """
    dest = Path(destination)
    vectors = dest / 'solver_vectors'
    vectors.mkdir(exist_ok=False)
    patch = load_patch(patch_root)
    function = patch.optimal_transport_duality_gap
    original_init = model_class.__init__
    record = dict(solver_replaced=True, patch_file_sha256=PATCH_SHA256,
                  parameters=EXPECTED, instances=[], calls=[], pairs=[], automatic_retry=False,
                  constructor_restored=False, diagnostics_only=True,
                  algorithm_source_unchanged=True, vectors_kept_remote=True)

    def new_init(instance, *args, **kwargs):
        # User-authorized analysis-only epsilon change; frozen source stays intact.
        requested_epsilon = kwargs['epsilon']
        assert requested_epsilon in (.02, .01), requested_epsilon
        effective_kwargs = dict(kwargs, epsilon=EXPECTED['epsilon'])
        assert all(effective_kwargs[k] is v for k, v in kwargs.items() if k != 'epsilon')
        original_init(instance, *args, **effective_kwargs)
        assert kwargs.get('solver', 'duality_gap') == 'duality_gap'
        for key, expected in EXPECTED.items():
            assert instance.ot_config[key] == expected, (key, instance.ot_config[key], expected)
        assert instance.ot_config['local_pca'] == 0, 'The 20D cost must remain unprojected'
        identifier = len(record['instances'])
        info = dict(instance=identifier, constructor_requested_epsilon=requested_epsilon,
                    effective_epsilon=EXPECTED['epsilon'], epsilon_only_constructor_override=True,
                    timepoints=list(map(float, instance.timepoints)),
                    config=dict(instance.ot_config), original_solver=instance.solver.__module__ + '.' + instance.solver.__name__,
                    calls=0)
        record['instances'].append(info)
        state = dict(pair=None, cost_hash=None, calls_in_pair=0, previous_rows=None)

        def solver(**params):
            started = time.monotonic()
            pair = _pair_frame(instance)
            for key, expected in EXPECTED.items():
                assert params[key] == expected, (key, params[key], expected)
            C, G = np.asarray(params['C']), np.asarray(params['G'])
            assert C.ndim == 2 and G.shape == (C.shape[0],)
            assert np.isfinite(C).all() and (C >= 0).all() and np.isfinite(G).all() and (G > 0).all()
            cost_hash, growth_hash = array_sha(C), array_sha(G)
            if state['pair'] != pair:
                assert state['pair'] is None or state['calls_in_pair'] == 3
                assert np.array_equal(G, np.ones(C.shape[0])), 'Original unit growth changed'
                state.update(pair=pair, cost_hash=cost_hash, calls_in_pair=0, previous_rows=None)
                record['pairs'].append(dict(instance=identifier, source_day=pair[0], target_day=pair[1],
                                            calls=[], cost_sha256=cost_hash, shape=list(C.shape)))
            else:
                assert state['calls_in_pair'] < 3 and state['cost_hash'] == cost_hash
                assert np.array_equal(G, state['previous_rows']), 'WOT growth is not the previous returned marginal'
            index = len(record['calls'])
            growth_call = state['calls_in_pair'] + 1
            item = dict(index=index, instance=identifier, source_day=pair[0], target_day=pair[1],
                        growth_call=growth_call, C_sha256=cost_hash, G_sha256=growth_hash,
                        shape=list(C.shape), parameters={k: params[k] for k in EXPECTED},
                        tolerance=float(params['tolerance']), max_iter=int(params['max_iter']),
                        batch_size=int(params['batch_size']), previous_marginal_exact=True,
                        vector_file=str(Path('solver_vectors') / ('call_%04d.npz' % index)))
            record['calls'].append(item)
            record['pairs'][-1]['calls'].append(index)
            info['calls'] += 1
            captured = {}
            captured_vectors = {}
            previous_profile = sys.getprofile()
            if previous_profile is not None:
                raise RuntimeError('Refusing to replace an existing Python profiler')

            def profile(frame, event, returned):
                if event != 'return' or frame.f_code is not function.__code__ or returned is None:
                    return
                values = frame.f_locals
                for key in ['epsilon_i', 'duality_gap', 'threshold', 'pri', 'dua']:
                    captured[key] = float(values[key])
                captured['current_iter'] = int(values['current_iter'])
                captured['final_scale'] = int(values['e'])
                captured['C_solver_sha256'] = array_sha(values['C'])
                captured['G_solver_sha256'] = array_sha(values['G'])
                captured['R_sha256'] = array_sha(values['R'])
                captured['returned_R_over_J_exact'] = bool(np.array_equal(returned, values['R'] / values['J']))
                for key in ['log_a', 'log_b', 'p', 'q']:
                    captured_vectors[key] = np.array(values[key], copy=True)

            try:
                sys.setprofile(profile)
                returned = function(**params)
            except Exception as error:
                item.update(status='failed', error=str(error), exception_type=type(error).__name__,
                            seconds=time.monotonic() - started)
                raise
            finally:
                sys.setprofile(previous_profile)
            assert captured, 'Stable solver return diagnostics were not captured'
            finite = bool(np.isfinite(list(captured[k] for k in ['epsilon_i', 'duality_gap', 'threshold', 'pri', 'dua'])).all())
            assert finite and captured['duality_gap'] <= captured['threshold'] == params['tolerance']
            assert np.isclose(captured['epsilon_i'], params['epsilon'], atol=1e-15, rtol=1e-14)
            assert captured['final_scale'] == 5 and captured['returned_R_over_J_exact']
            assert captured['C_solver_sha256'] == cost_hash and captured['G_solver_sha256'] == growth_hash
            assert np.isfinite(returned).all() and (returned >= 0).all() and returned.sum() > 0
            rows = np.asarray(returned.sum(axis=1))
            columns = np.asarray(returned.sum(axis=0))
            assert all(np.isfinite(x).all() for x in captured_vectors.values())
            np.savez(vectors / ('call_%04d.npz' % index), input_G=G, returned_row_sums=rows,
                     returned_column_sums=columns, **captured_vectors)
            item.update(status='complete', **captured, primal_dual_finite=finite, tolerance_met=True,
                        returned_sha256=array_sha(returned), returned_row_sums_sha256=array_sha(rows),
                        returned_column_sums_sha256=array_sha(columns),
                        returned_total=float(returned.sum()), input_G_min=float(G.min()), input_G_max=float(G.max()),
                        returned_rows_min=float(rows.min()), returned_rows_max=float(rows.max()),
                        vector_file_sha256=sha(vectors / ('call_%04d.npz' % index)), seconds=time.monotonic() - started)
            state.update(calls_in_pair=growth_call, previous_rows=rows.copy())
            return returned

        instance.solver = solver

    model_class.__init__ = new_init
    try:
        yield record
    finally:
        model_class.__init__ = original_init
        record['constructor_restored'] = model_class.__init__ is original_init
        write(dest / 'solver_diagnostics.json', record)


def verify_saved_maps(record, out, read_h5ad, expected_pairs):
    """Check the final saved map is growth call3's R/J, including growth marginals."""
    out = Path(out)
    assert record['constructor_restored'] and len(record['instances']) == 1
    assert len(record['pairs']) == expected_pairs and len(record['calls']) == expected_pairs * 3
    for i, pair in enumerate(record['pairs']):
        assert pair['source_day'] == float(i) and pair['target_day'] == float(i + 1)
        assert len(pair['calls']) == 3
        calls = [record['calls'][j] for j in pair['calls']]
        assert [c['growth_call'] for c in calls] == [1, 2, 3]
        assert all(c['status'] == 'complete' and c['tolerance_met'] for c in calls)
        for previous, current in zip(calls, calls[1:]):
            assert previous['returned_row_sums_sha256'] == current['G_sha256']
            with np.load(out / previous['vector_file']) as a, np.load(out / current['vector_file']) as b:
                assert np.array_equal(a['returned_row_sums'], b['input_G'])
        path = out / ('OT_%d_%d.h5ad' % (i, i + 1))
        saved = read_h5ad(path)
        matrix = saved.X.toarray() if hasattr(saved.X, 'toarray') else np.asarray(saved.X)
        assert array_sha(matrix) == calls[-1]['returned_sha256'], 'Saved WOT map is not the third growth return'
        delta_days = pair['target_day'] - pair['source_day']
        assert delta_days > 0
        for j, call in enumerate(calls):
            with np.load(out / call['vector_file']) as values:
                expected_metadata = np.power(values['input_G'], 1.0 / delta_days)
                assert np.array_equal(np.asarray(saved.obs['g' + str(j)]), expected_metadata)
        with np.load(out / calls[-1]['vector_file']) as values:
            expected_metadata = np.power(values['returned_row_sums'], 1.0 / delta_days)
            assert np.array_equal(np.asarray(saved.obs['g3']), expected_metadata)
        pair.update(final_map=str(path), final_map_file_sha256=sha(path),
                    final_map_array_sha256=array_sha(matrix), final_map_is_third_return_R_over_J=True,
                    learned_growth_vectors_exact=True,
                    learned_growth_metadata_transform='np.power(raw_growth, 1.0 / delta_days)',
                    delta_days=delta_days)
    record['all_final_maps_verified'] = True
    write(out / 'solver_diagnostics.json', record)
