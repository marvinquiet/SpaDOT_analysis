"""Fit an original loss/architecture ablation with its archived OT solver.

Supply the original study manifest, source, historical model modules, configs,
evaluation neighbor graphs and preprocessed H5ADs separately. Set
SPADOT_LOSS_STUDY_ROOT to that workspace and pass its job ID as the first
argument. Use the matching CUDA environment. See ../training/worker.py for
the separate numerical-recovery run used for the Chicken alignment-off result.
"""
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[key] = '1'
os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'
import sys
import copy
import inspect
import json
import time
import traceback
from pathlib import Path
from common import ROOT, sha, check_frozen, write_json
sys.path.insert(0, str(ROOT / 'source'))
import torch
import numpy as np
import pandas as pd
import anndata as ad
import yaml
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits
from torch_geometric.utils import dense_to_sparse
from SpaDOT.utils import _train_utils as trainer, _utils, _analyze_utils as analysis
from SpaDOT.model.SpaDOT import SpaDOT as CurrentSpaDOT
from historical_interface import compatible_class, original_class, raw_inference

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'patches'))
from skip_inactive_ot import install
install(trainer)

NATIVE_UPDATE = trainer._update_Kmeans
NATIVE_OT = trainer._update_OT_matrix


def loaders_for(original, config, times):
    full = copy.deepcopy(config)
    full['timepoints'] = sorted(original.obs.timepoint.unique().tolist())
    full['input_dim'] = original.n_vars
    full.update(dtype=torch.float64, device=torch.device('cuda:0'))
    _utils.set_seed(2026)
    data = trainer.prepare_dataloader(original, full)
    scoped = {key: {tp: values[tp] for tp in times} for key, values in data.items()}
    full['timepoints'] = list(times)
    return full, scoped


def score_case(latent, original, info, times):
    with threadpool_limits(limits=1):
        ks = [info['final_clusters'][info['timepoints'].index(t)] for t in times]
        pred = analysis.KMeans_Clustering(latent, ks, seed=2026, n_init=10)
    graph = np.load(ROOT / 'evaluation/fixed_spatial_neighbors.npz')
    scores = {}
    for tp in times:
        p = pred[pred.obs.timepoint.astype(str) == tp]
        source = original[original.obs.timepoint.astype(str) == tp]
        assert np.array_equal(p.obs_names, source.obs_names)
        assert np.array_equal(p.obsm['spatial'], source.obsm['spatial'])
        labels = p.obs.kmeans.astype(str).to_numpy()
        key = info['dataset'] + '__' + tp
        neighbors = graph[key]
        assert neighbors.shape == (p.n_obs, 10)
        assert not np.any(neighbors == np.arange(p.n_obs)[:, None])
        scores[tp] = dict(ARI=float(adjusted_rand_score(p.obs.annotation.astype(str), labels)),
                         PAS=float(np.mean(np.sum(labels[neighbors] != labels[:, None], axis=1) >= 6)),
                         spots=p.n_obs, actual_clusters=int(p.obs.kmeans.nunique()))
    return pred, scores


def main():
    manifest = json.loads((ROOT / 'manifest.json').read_text())
    job = next(j for j in manifest['jobs'] if j['id'] == sys.argv[1])
    info = dict(manifest['datasets'][job['dataset']], dataset=job['dataset'])
    out = ROOT / 'runs' / job['id']; out.mkdir(exist_ok=False)
    write_json(out / 'started.json', dict(job=job, pid=os.getpid(), started=time.time(), physical_gpu=0))
    stage = 'setup'
    try:
        torch.set_num_threads(1); torch.set_num_interop_threads(1)
        torch.use_deterministic_algorithms(True)
        assert torch.cuda.is_available()
        check_frozen(manifest)
        original = ad.read_h5ad(info['input'])
        assert original.obs_names.is_unique and original.shape == (info['spots'], info['features'])
        base = yaml.safe_load(Path(job['config']).read_text())
        groups = [info['timepoints']] if job['mode'] == 'combined' else [[t] for t in info['timepoints']]
        pieces = []; fit_results = []
        for times in groups:
            tag = 'combined' if len(times) > 1 else times[0]
            directory = out / tag; directory.mkdir()
            cfg, data = loaders_for(original, base, times)
            trainer.SpaDOT.SpaDOT = CurrentSpaDOT if job['method']['family'] == 'full' else compatible_class(job['method'])
            updates = []; ot_events = []; skipped_ot_events = []

            def observe(model, config, loaders):
                NATIVE_UPDATE(model, config, loaders)
                frame = inspect.currentframe().f_back
                epoch = int(frame.f_locals['epoch'])
                values = {k: float(v) for k, v in frame.f_locals['loss_dict'][epoch].items()}
                assert np.isfinite(list(values.values())).all()
                steps = max(int(v['step'].item()) for v in frame.f_locals['optimizer'].state.values() if 'step' in v)
                updates.append(dict(epoch=epoch, optimizer_steps=steps))
                write_json(out / 'progress.json', dict(fit=tag, timepoints=times, epoch=epoch,
                    completed_fits=len(fit_results), total_fits=len(groups), losses=values, updated=time.time()))

            def observe_ot(model, config):
                if config['omiga3'] == 0:
                    assert NATIVE_OT(model, config) is None
                    skipped_ot_events.append(updates[-1]['epoch'])
                    return
                try:
                    NATIVE_OT(model, config)
                    assert all(np.isfinite(v).all() and (v >= 0).all() for v in model.gammas.values())
                    ot_events.append(updates[-1]['epoch'])
                except BaseException:
                    torch.save(model.state_dict(), directory / 'failure_checkpoint.pth')
                    np.savez_compressed(directory / 'failure_ot_centers.npz', **model.kmeans_center_dict)
                    raise

            trainer._update_Kmeans = observe; trainer._update_OT_matrix = observe_ot
            stage = 'training'
            with threadpool_limits(limits=1):
                model, losses = trainer.train_SpaDOT(data, cfg)
            stage = 'verification'
            assert [v['epoch'] for v in updates] == list(range(1, 100))
            expected_steps = 99 * sum(len(data['dataloaders'][t]) for t in times)
            assert updates[-1]['optimizer_steps'] == expected_steps
            expected_refresh = list(range(10, 100, 10))
            if cfg['omiga3'] == 0:
                assert not ot_events and skipped_ot_events == expected_refresh
                assert not model.gammas
            else:
                assert ot_events == expected_refresh and not skipped_ot_events
            loss = losses.T
            assert loss.index.tolist() == list(range(1, 100)) and np.isfinite(loss.values).all()
            beta = trainer._beta_cycle_linear_period(99, 30, stop=base['beta1'])
            rebuilt_loss = (base['lambda1'] * loss.Recon - beta * loss.SVGP_KL + base['beta2'] * loss.GAT_KL
                + base['omiga1'] * loss.alignment + base['omiga2'] * loss.KMeans + base['omiga3'] * loss.OT)
            assert np.allclose(rebuilt_loss, loss.elbo, rtol=1e-10, atol=1e-8)
            if cfg['omiga3'] == 0:
                assert np.array_equal(loss.OT.to_numpy(), np.zeros(99))
            loss.to_csv(directory / 'loss.csv')
            torch.save(model.state_dict(), directory / 'SpaDOT_model.pth')
            assert all(torch.isfinite(v).all() for v in model.state_dict().values() if v.is_floating_point())
            latent = trainer.get_latent(model, cfg, original, data)
            latent.write_h5ad(directory / 'latent.h5ad')
            assert latent.n_vars == 20 and np.isfinite(latent.X).all()
            # Direct original-class checkpoint readout, not the compatibility subclass.
            cls = CurrentSpaDOT if job['method']['family'] == 'full' else original_class(job['method'])
            fresh_cfg, fresh_data = loaders_for(original, base, times)
            fresh = cls(fresh_cfg, fresh_data).to('cuda:0')
            state = torch.load(directory / 'SpaDOT_model.pth', map_location='cpu', weights_only=True)
            fresh.load_state_dict(state); fresh.eval(); arrays = []
            with torch.no_grad():
                for tp in times:
                    X, Y, _ = fresh_data['datasets'][tp].dataset.tensors
                    edge, _ = dense_to_sparse(fresh_data['adjacency_matrices'][tp])
                    z = fresh.all_latent_samples(X, Y, edge, tp) if job['method']['family'] == 'full' else raw_inference(fresh, job['method'], X, Y, edge, tp)
                    arrays.append(np.asarray(z))
            reproduced = np.concatenate(arrays)
            assert np.array_equal(reproduced, latent.X)
            pred, scores = score_case(latent, original, info, times)
            pred.obs[['timepoint', 'annotation', 'kmeans']].to_csv(directory / 'predictions.csv')
            fit = dict(all_passed=True, timepoints=times, executed_epoch=99, selected_epoch=99,
                scores=scores, optimizer_steps=expected_steps, source_commit=manifest['commit'],
                original_class_checkpoint_exact=True, checkpoint_max_abs_error=0,
                loss_objective_max_abs_error=float(np.max(np.abs(rebuilt_loss-loss.elbo))),
                trainable_parameters=sum(p.numel() for p in model.parameters() if p.requires_grad),
                inactive_ot_guard=bool(cfg['omiga3'] == 0), inactive_ot_skipped_epochs=skipped_ot_events,
                ot_refresh_epochs=ot_events, output_hashes={p.name:sha(p) for p in directory.iterdir() if p.is_file()})
            write_json(directory / 'result.json', fit)
            fit_results.append(fit); pieces.append(latent)
            del model, fresh, data, fresh_data, state
            trainer._update_Kmeans = NATIVE_UPDATE; trainer._update_OT_matrix = NATIVE_OT
            torch.cuda.empty_cache()
        stage = 'evaluation'
        latent = ad.concat(pieces)
        pred, scores = score_case(latent, original, info, info['timepoints'])
        pred.obs[['timepoint', 'annotation', 'kmeans']].to_csv(out / 'predictions.csv')
        pd.DataFrame(pred.obsm['spatial'], index=pred.obs_names, columns=['pixel_x', 'pixel_y']).to_csv(out / 'coordinates.csv')
        stage = 'verification'; check_frozen(manifest)
        result = dict(all_passed=True, job=job, source_commit=manifest['commit'], executed_epoch=99, selected_epoch=99,
            mode=job['mode'], neural_parameter_sets=len(fit_results), scores=scores,
            mean_ARI=float(np.mean([s['ARI'] for s in scores.values()])),
            mean_PAS=float(np.mean([s['PAS'] for s in scores.values()])),
            original_class_checkpoint_exact=True, all_fits_verified=True, fit_results=fit_results,
            retries=0, study_version='skipinactiveot_v1', original_trial_id=job['id'],
            inactive_ot_guard=bool(base['omiga3'] == 0), output_hashes={p.name:sha(p) for p in out.iterdir() if p.is_file()})
        write_json(out / 'result.json', result)
        print('VERIFIED_RESULT', json.dumps(result), flush=True)
    except BaseException:
        write_json(out / 'failure.json', dict(stage=stage, error=traceback.format_exc(), updated=time.time()))
        raise


if __name__ == '__main__':
    main()
