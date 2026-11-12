"""Proposed study-only guard; do not install into an active frozen worker."""
def install(trainer):
    native_update = trainer._update_OT_matrix
    native_loss = trainer._compute_OT_loss

    def update(model, config):
        if config['omiga3'] == 0:
            return None
        return native_update(model, config)

    def loss(model, config, cur_tp, tp_ix, latent, prev_tp):
        if config['omiga3'] == 0:
            return latent.new_zeros(())
        return native_loss(model, config, cur_tp, tp_ix, latent, prev_tp)

    trainer._update_OT_matrix = update
    trainer._compute_OT_loss = loss
    return native_update, native_loss


if __name__ == '__main__':
    import sys
    import json
    from pathlib import Path
    from types import SimpleNamespace
    import numpy as np
    import torch
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / 'source'))
    from SpaDOT.utils import _train_utils as trainer
    torch.set_num_threads(1); torch.manual_seed(2026)
    points = torch.randn(6, 20, dtype=torch.float64)
    cfg = dict(dtype=torch.float64, device=torch.device('cpu'), omiga3=0.)
    rng = np.random.default_rng(2026)
    model = SimpleNamespace(kmeans_index_dict={'t1': {i:i%3 for i in range(6)}},
        kmeans_cluster_dict={'t1':[0,1,2,0,1,2]},
        kmeans_center_dict={'t0':rng.normal(size=(3,20)), 't1':rng.normal(size=(3,20))},
        gammas={'t0_t1':np.full((3,3),1/3)})
    actual = points.clone().requires_grad_()
    cost = trainer._compute_OT_loss(model,cfg,'t1',torch.arange(6),actual,'t0')
    native_objective = actual.square().sum() + 0*cost
    native_objective.backward()
    native_update, native_loss = install(trainer)
    guarded = points.clone().requires_grad_()
    guarded_objective = guarded.square().sum() + 0*trainer._compute_OT_loss(model,cfg,'t1',torch.arange(6),guarded,'t0')
    guarded_objective.backward()
    assert torch.equal(native_objective, guarded_objective)
    assert torch.equal(actual.grad, guarded.grad)
    assert trainer._update_OT_matrix(model,cfg) is None
    nonzero=dict(cfg,omiga3=1.)
    assert torch.equal(native_loss(model,nonzero,'t1',torch.arange(6),guarded,'t0'),
                       trainer._compute_OT_loss(model,nonzero,'t1',torch.arange(6),guarded,'t0'))
    result=dict(all_passed=True,finite_zero_weight_fixture_loss_exact=True,
        finite_zero_weight_fixture_gradients_exact=True,nonzero_ot_loss_delegated_unchanged=True,
        real_training_steps=0,scope='Proposed helper only; no active frozen code modified')
    (root/'proposed/guard_fixture_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
