"""Run one S39 native-sampler ablation fit with an external frozen study.

Set SPADOT_ENCODER_STUDY_ROOT to the separately prepared parent or depth-four
study workspace. Supply its frozen source, manifest, configs, provenance
files and numerical solver patch separately. Inputs are the preprocessed H5ADs
recorded in the manifest. Use the matching SpaDOT/Torch/PyG CUDA environment.
Set PYTHONHASHSEED to the trial seed before launching Python.
The script writes a new trial directory under --output-dir.
"""
import os
for _key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[_key] = '1'
os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'
from pathlib import Path
import argparse, importlib.metadata, json, sys, time, traceback, warnings
from types import SimpleNamespace

ROOT = Path(os.environ['SPADOT_ENCODER_STUDY_ROOT']).resolve()
sys.path.insert(0, str(ROOT / 'source'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
import anndata as ad
import torch
import yaml
from threadpoolctl import threadpool_limits
import SpaDOT
from SpaDOT.utils import _train_utils as U, _utils
from runtime import installed, sha, read, write, verify_native_loaders
from evaluate import evaluate


def _finite_parameters(model):
    return all(bool(torch.isfinite(p).all()) for p in model.parameters())


def _verify_losses(loss, config, family):
    """Independently replay scalar coefficients on the saved raw loss columns."""
    assert list(map(int, loss.index)) == list(range(1, 100))
    assert np.isfinite(loss.to_numpy()).all()
    beta = U._beta_cycle_linear_period(99, config['kl_cycle_epochs'], stop=config['beta1'])
    terms = pd.DataFrame(index=loss.index)
    terms['weighted_reconstruction'] = config['lambda1'] * loss.Recon
    terms['weighted_SVGP_KL'] = -beta * loss.SVGP_KL.to_numpy()
    terms['weighted_GAT_KL'] = config['beta2'] * loss.GAT_KL
    terms['weighted_cross_encoder_alignment'] = config['omiga1'] * loss.alignment
    terms['weighted_KMeans'] = config['omiga2'] * loss.KMeans
    terms['weighted_OT'] = config['omiga3'] * loss.OT
    terms['replayed_objective'] = terms.sum(axis=1)
    terms['saved_objective'] = loss.elbo
    terms['objective_absolute_error'] = np.abs(terms.replayed_objective - loss.elbo)
    terms['SVGP_beta1'] = beta
    assert np.allclose(terms.replayed_objective, loss.elbo, rtol=1e-11, atol=1e-7)
    if family in ('gat_only', 'mlp_gat', 'gat_gat_matched'):
        assert (loss.SVGP_KL == 0).all()
    if family == 'gat_only':
        assert (loss.alignment == 0).all()
    assert float(loss.KMeans.iloc[0]) == 0.
    assert (loss.OT.iloc[:config['ot_epoch']] == 0).all()
    return terms


def _verify_inputs(case, latent):
    assert latent.shape == (case['shape'][0], 20)
    assert latent.obs_names.is_unique and np.isfinite(latent.X).all()
    original = ad.read_h5ad(case['input'], backed='r')
    try:
        assert set(original.obs_names) == set(latent.obs_names)
        reference = original[latent.obs_names]
        assert np.array_equal(reference.obsm['spatial'], latent.obsm['spatial'])
        for key in ['timepoint', 'annotation']:
            assert reference.obs[key].astype(str).tolist() == latent.obs[key].astype(str).tolist()
    finally:
        original.file.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--trial', required=True)
    parser.add_argument('--gpu', type=int, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    manifest = read(ROOT / 'manifest.json')
    trial = next(t for t in manifest['trials'] if t['id'] == args.trial)
    case = next(c for c in manifest['cases'] if c['id'] == trial['dataset'])
    out = args.output_dir / trial['id']
    out.mkdir(parents=True, exist_ok=False)
    status = dict(id=trial['id'], dataset=trial['dataset'], arm=trial['arm'],
                  family=trial['family'], depth=trial['depth'], width=trial['width'],
                  seed=trial['seed'], gpu=args.gpu, pid=os.getpid(), status='running',
                  stage='preflight', started_unix=time.time())
    write(out / 'status.json', status)
    def stage(name):
        status['stage'] = name
        write(out / 'status.json', status)
        print('STAGE', name, flush=True)
    try:
        assert Path(SpaDOT.__file__).resolve() == ROOT / 'source/SpaDOT/__init__.py'
        protected = {str(ROOT / 'source' / n): h for n, h in manifest['source_sha256'].items()}
        protected.update({str(ROOT / n): h for n, h in manifest['code_sha256'].items()})
        cfgp = Path(trial['config'])
        assert cfgp.parent == ROOT / 'configs'
        protected[str(cfgp)] = manifest['config_sha256'][cfgp.name]
        protected[str(ROOT / 'manifest.json')] = sha(ROOT / 'manifest.json')
        protected[str(ROOT / 'native_reference_inventory.json')] = manifest['native_reference_inventory_sha256']
        protected[case['input']] = case['input_sha256']
        protected.update({str(ROOT/n): h for n,h in manifest['patch_sha256'].items()})
        for path, expected in protected.items():
            assert sha(path) == expected, path
        packages = {name: importlib.metadata.version(name) for name in manifest['package_versions']}
        assert packages == manifest['package_versions']
        assert os.environ.get('PYTHONHASHSEED') == str(trial['seed'])
        config = yaml.safe_load(cfgp.read_text())
        reference_config = yaml.safe_load((ROOT / 'source/SpaDOT/config.yaml').read_text())
        allowed = {'seed', 'gat_depth', 'gat_encoder_hidden', 'gat_capacity_control',
                   'study_family', 'study_sampling_hops', 'replacement_width', 'knn_cutoff', 'max_neighbors'}
        assert {k for k in config if config[k] != reference_config.get(k)} <= allowed
        assert (config['epoch_start'], config['epoch_stop_exclusive'], config['kl_cycle_epochs'],
                config['z_dim']) == (1, 100, 30, 20)
        assert (config['gat_depth'], config['gat_encoder_hidden'], config['study_family'],
                config['study_sampling_hops']) == (trial['depth'], trial['width'], trial['family'], 2)
        assert config['replacement_width']==trial['replacement_width']
        assert (config['knn_cutoff'],config['max_neighbors']) == ((12,60) if trial['graph_expansion'] else (6,30))
        assert config['kernel_scale'] == .1 and not config['gat_capacity_control']
        warnings.filterwarnings('ignore', message="Can't initialize NVML")
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
        torch.use_deterministic_algorithms(True)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        torch.cuda.set_device(args.gpu)
        device = torch.device(f'cuda:{args.gpu}')
        eye = torch.eye(2, dtype=torch.float64, device=device)
        assert torch.equal(eye @ eye, eye)
        free, total = torch.cuda.mem_get_info(args.gpu)
        del eye
        torch.cuda.reset_peak_memory_stats(device)
        write(out / 'worker_preflight.json', dict(packages=packages, actual_CUDA_passed=True,
              gpu=args.gpu, free_bytes_before_model=free, total_gpu_bytes=total,
              protected_files=len(protected), source_and_input_hashes_verified=True,
              config_sha256=sha(cfgp), manifest_sha256=sha(ROOT / 'manifest.json')))
        (out / 'config_used.yaml').write_text(cfgp.read_text())
        observer, holder = {'capture_batch_hashes':True}, {}
        original_train, original_update = U.train_SpaDOT, U._update_Kmeans
        branch_history = []
        def observed_train(dl, cfg):
            holder.update(dl=dl, config=cfg)
            write(out/'native_sampler.json',verify_native_loaders(dl,cfg))
            model, loss = original_train(dl, cfg)
            holder['model'] = model
            assert sum(p.numel() for p in model.parameters() if p.requires_grad) == trial['expected_parameters']
            assert _finite_parameters(model)
            return model, loss
        def observed_update(model, cfg, dl):
            result = original_update(model, cfg, dl)
            parent = sys._getframe(1).f_locals
            if 'epoch' in parent and 'loss_dict' in parent:
                epoch = int(parent['epoch'])
                values = {key: float(v) for key, v in parent['loss_dict'][epoch].items()}
                if not np.isfinite(list(values.values())).all():
                    raise FloatingPointError('Non-finite recorded training loss')
                branch = getattr(model, '_branch_epoch_diagnostics', None)
                if branch is not None:
                    record = dict(epoch=epoch, **branch)
                    branch_history.append(record)
                    pd.DataFrame(branch_history).to_csv(out / 'branch_KL_progress.csv', index=False, float_format='%.17g')
                    model._branch_epoch_diagnostics = None
                write(out / 'progress.json', dict(epoch=epoch, losses=values,
                      updated_unix=time.time(), fixed_executed_epoch=99,
                      selected_epoch=99, no_early_stopping=True))
                # This is a diagnostic scalar history, not a resumable checkpoint.
                partial = pd.DataFrame.from_dict({int(e): row for e, row in parent['loss_dict'].items()
                                                  if int(e) <= epoch}, orient='index')
                partial.to_csv(out / 'loss_progress.csv', float_format='%.17g')
            return result
        from training_solver_guard import observe_training_solver
        U.train_SpaDOT, U._update_Kmeans = observed_train, observed_update
        start = time.monotonic()
        try:
            with threadpool_limits(limits=1), installed(U, trial, observer), observe_training_solver(ROOT,out/'training_solver_trace') as training_ot:
                stage('training')
                SpaDOT.train(SimpleNamespace(data=case['input'], output_dir=str(out), prefix='',
                             config=str(cfgp), save_model=True, device=str(device)))
        finally:
            U.train_SpaDOT, U._update_Kmeans = original_train, original_update
        training_seconds = time.monotonic() - start
        assert training_ot['all_recorded_calls_returned'] and training_ot['original_function_restored'] and training_ot['calls']
        model, used, dl = holder['model'], holder['config'], holder['dl']
        assert _finite_parameters(model)
        actual_parameters = sum(p.numel() for p in model.parameters() if p.requires_grad)
        assert actual_parameters == trial['expected_parameters']
        if trial['family']=='gat_only':
            assert not hasattr(model,'svgp_dict') and not hasattr(model,'SVGPEncoder')
            assert model.GAT_z_dim==20
        elif trial['family'] in ('mlp_gat','gat_gat_matched'):
            assert not hasattr(model,'svgp_dict')
            assert model.SVGP_z_dim==0 and model.replacement_z_dim==model.GAT_z_dim==10
        else:
            assert trial['family']=='dual' and model.SVGP_z_dim==model.GAT_z_dim==10
        write(out/'first_training_batches.json',observer['first_training_batches'])
        initialization = observer['initialization']
        write(out / 'initialization.json', initialization)
        latent = ad.read_h5ad(out / 'latent.h5ad')
        _verify_inputs(case, latent)
        loss = pd.read_csv(out / 'loss.csv', index_col=0, float_precision='round_trip')
        objective_terms = _verify_losses(loss, config, trial['family'])
        if trial['family'] in ('mlp_gat','gat_gat_matched'):
            assert len(branch_history)==99
            for record in branch_history:
                assert np.isclose(record['KL_A_epoch_loss']+record['KL_B_epoch_loss'],loss.loc[record['epoch'],'GAT_KL'],atol=1e-8,rtol=1e-11)
            pd.DataFrame(branch_history).to_csv(out/'branch_KL.csv',index=False,float_format='%.17g')
        objective_terms.to_csv(out / 'weighted_objective_replay.csv', float_format='%.17g')
        checkpoint = out / 'SpaDOT_model.pth'
        state = torch.load(checkpoint, map_location=device, weights_only=True)
        assert all(bool(torch.isfinite(value).all()) for value in state.values())
        stage('checkpoint_reconstruction')
        # Release the training graph before independently rebuilding the model.
        del model, dl
        holder.clear()
        torch.cuda.empty_cache()
        original_input = ad.read_h5ad(case['input'])
        rebuild_observer = {}
        _utils.set_seed(trial['seed'])
        with threadpool_limits(limits=1), installed(U, trial, rebuild_observer):
            rebuilt_dl = U.prepare_dataloader(original_input, used)
            reconstructed_sampler=verify_native_loaders(rebuilt_dl,used)
            assert reconstructed_sampler==read(out/'native_sampler.json')
            rebuilt = U.SpaDOT.SpaDOT(used, rebuilt_dl).to(device)
            rebuilt.load_state_dict(state, strict=True)
            reconstructed = U.get_latent(rebuilt, used, original_input, rebuilt_dl)
        assert np.array_equal(reconstructed[latent.obs_names].X, latent.X), 'Saved latent cannot be reproduced exactly'
        assert np.array_equal(reconstructed[latent.obs_names].obsm['spatial'], latent.obsm['spatial'])
        assert _finite_parameters(rebuilt)
        reconstruction_parameters = sum(p.numel() for p in rebuilt.parameters() if p.requires_grad)
        assert reconstruction_parameters == actual_parameters
        write(out / 'checkpoint_reconstruction.json', dict(all_passed=True,
              saved_latent_exact=True, checkpoint_sha256=sha(checkpoint), latent_sha256=sha(out / 'latent.h5ad'),
              actual_parameters=actual_parameters, family=trial['family'], latent_dimensions=20,
              native_sampler=reconstructed_sampler))
        native_reference_verified=None
        if trial['family']=='dual':
            reference=trial['native_reference']
            stage('native_reference_verification')
            for key in ['checkpoint','latent','loss','config']:
                assert sha(reference[key])==reference[key+'_sha256']
            historical=ad.read_h5ad(reference['latent'])[latent.obs_names].copy()
            oldloss=pd.read_csv(reference['loss'],index_col=0,float_precision='round_trip')
            assert list(oldloss.index)==list(loss.index)==list(range(1,100))
            assert set(oldloss.columns)==set(loss.columns)
            oldloss=oldloss[loss.columns]
            delta=np.abs(oldloss.to_numpy()-loss.to_numpy())
            loss_exact=bool(np.array_equal(oldloss.to_numpy(),loss.to_numpy()))
            latent_exact=bool(np.array_equal(historical.X,latent.X))
            historical_state=torch.load(reference['checkpoint'],map_location='cpu',weights_only=True)
            state_exact=set(historical_state)==set(state) and all(torch.equal(v.detach().cpu(),historical_state[k].detach().cpu()) for k,v in state.items())
            native_reference_verified=dict(all_passed=True,
                 comparison_completed=True,historical_exact_required=False,
                 historical_all_exact=bool(loss_exact and latent_exact and state_exact),
                 numerical_patch_only=True,training_solver_sha256=manifest['protocol']['training_solver']['sha256'],
                loss_all99_rows_exact=loss_exact,loss_max_absolute_error=float(delta.max()),
                latent_exact=latent_exact,state_dict_tensor_values_exact=bool(state_exact),
                original_reference=reference,new_checkpoint_sha256=sha(checkpoint),
                new_latent_sha256=sha(out/'latent.h5ad'))
            write(out/'native_reference_verification.json',native_reference_verified)
            assert np.isfinite(delta).all(), 'Nonfinite historical comparison'
            assert native_reference_verified['comparison_completed']
            del historical,historical_state,oldloss,delta
        del rebuilt, rebuilt_dl, reconstructed, original_input, state
        torch.cuda.empty_cache()
        stage('domain_and_association_evaluation')
        with threadpool_limits(limits=1):
            evaluation = evaluate(ROOT, out, case, trial, latent)
        for path, expected in protected.items():
            assert sha(path) == expected, path
        report = dict(all_training_checks_passed=True, all_domain_checks_passed=evaluation['domain_passed'],
                      association_passed=evaluation['association_passed'], executed_epoch=99, selected_epoch=99,
                      no_early_stopping=True, checkpoint_reconstruction_exact=True, trial=trial,
                      training_solver_sha256=training_ot['solver_sha256'],training_solver_calls=len(training_ot['calls']),
                      input_sha256=case['input_sha256'], actual_parameters=actual_parameters,
                      initialization=initialization,native_sampler=read(out/'native_sampler.json'),
                      native_reference_verified=native_reference_verified,checkpoint_sha256=sha(checkpoint),
                      protected_files_unchanged=len(protected), protected_sha256=protected,
                      objective_replay_max_error=float(objective_terms.objective_absolute_error.max()),
                      max_GPU_allocated_bytes=int(torch.cuda.max_memory_allocated(device)),
                      max_GPU_reserved_bytes=int(torch.cuda.max_memory_reserved(device)),
                      training_seconds=training_seconds,
                      objective_notes=manifest['protocol']['objectives'][trial['family']])
        write(out / 'verification.json', report)
        primary = evaluation['combined_domain_metrics']
        analysis_ok = evaluation['domain_passed'] and evaluation['association_passed'] and evaluation['domain_plots_passed']
        status.update(status='complete' if analysis_ok else 'complete_with_analysis_failure', stage='terminal',
                      executed_epoch=99, selected_epoch=99,
                      ARI={row['timepoint']: row['ARI'] for row in primary},
                      mean_ARI=float(np.mean([row['ARI'] for row in primary])) if primary else None,
                      domain_passed=evaluation['domain_passed'], association_passed=evaluation['association_passed'],
                      finished_unix=time.time())
        write(out / 'status.json', status)
        files = {p.name: sha(p) for p in out.iterdir() if p.is_file() and p.name not in ['status.json', 'complete.json']}
        write(out / 'complete.json', dict(training_passed=True, domain_passed=evaluation['domain_passed'],
              association_passed=evaluation['association_passed'], files=files))
        print('COMPLETE', json.dumps(status), flush=True)
    except Exception as error:
        text = traceback.format_exc()
        from training_solver_guard import numerical_failure
        failure=dict(status, error=str(error), exception_type=type(error).__name__, traceback=text)
        recognized_numeric=numerical_failure(failure,ROOT,trial['family'])
        status.update(status='failed', error=str(error), exception_type=type(error).__name__,
                      error_category='documented_training_numerical' if recognized_numeric else 'unexpected',
                      recognized_training_numerical_failure=recognized_numeric, traceback=text, finished_unix=time.time(),
                      automatic_retry=False)
        write(out / 'status.json', status)
        print(text, flush=True)
        raise


if __name__ == '__main__':
    main()
