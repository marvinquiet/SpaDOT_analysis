"""Isolated two-branch replacements; no frozen SpaDOT source is modified.

Both replacement families concatenate Gaussian branch A10 and the original
GAT B10 into the original decoder20. Their two standard-Gaussian KL terms are
separately normalized by ten, summed into the trainer's GAT_KL output, and
weighted by its unchanged beta2. There is no GP object or SVGP loss. The
original sampled-norm alignment, reconstruction, KMeans and temporal OT remain.
"""
import hashlib
from types import SimpleNamespace

import torch
from torch import nn
from torch.nn import functional as F
from torch_geometric.nn import GATConv

from SpaDOT.model.decoder import Decoder
from SpaDOT.model.encoder import SVGPEncoder as GaussianMLPEncoder


def _positive_int(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(name + ' must be a positive integer')
    return value


def parameter_count(input_dim, width, depth, GAT_z_dim=20, num_heads=4):
    """Exact active GAT-encoder count, excluding the shared decoder."""
    for value, name in ((input_dim, 'input_dim'), (width, 'width'),
                        (depth, 'depth'), (GAT_z_dim, 'GAT_z_dim'),
                        (num_heads, 'num_heads')):
        _positive_int(value, name)
    if depth > 6:
        raise ValueError('depth must be at most six')
    hw = num_heads * width
    return (input_dim * hw + (depth - 1) * hw ** 2 + depth * 2 * hw
            + (depth - 1) * hw + width + 2 * GAT_z_dim * width + 2 * GAT_z_dim)


def mlp_parameter_count(input_dim, hidden_dims, latent_dim=10):
    """Source Gaussian network, including active BatchNorm affine parameters."""
    _positive_int(input_dim, 'input_dim')
    _positive_int(latent_dim, 'latent_dim')
    hidden_dims = list(hidden_dims)
    if not hidden_dims:
        raise ValueError('hidden_dims cannot be empty')
    total, previous = 0, input_dim
    for width in hidden_dims:
        _positive_int(width, 'hidden width')
        total += previous * width + width + 2 * width
        previous = width
    return total + previous * 2 * latent_dim + 2 * latent_dim


def state_digest(state):
    """Hash named tensor shapes, dtypes and values, including BN buffers."""
    digest = hashlib.sha256()
    for name, value in sorted(state.items()):
        value = value.detach().cpu().contiguous()
        digest.update(name.encode())
        digest.update(str(value.dtype).encode())
        digest.update(str(tuple(value.shape)).encode())
        digest.update(value.numpy().tobytes())
    return digest.hexdigest()


def _rng_digest(state):
    return hashlib.sha256(state.detach().cpu().numpy().tobytes()).hexdigest()


def _count(module):
    return sum(p.numel() for p in module.parameters() if p.requires_grad)


class FlexibleGATEncoder(nn.Module):
    """Original graph stages and Gaussian head, without inactive parameters."""
    def __init__(self, input_dim, GAT_z_dim, hidden_dim=512, num_heads=4, depth=3):
        super().__init__()
        for value, name in ((input_dim, 'input_dim'), (GAT_z_dim, 'GAT_z_dim'),
                            (hidden_dim, 'hidden_dim'), (num_heads, 'num_heads')):
            _positive_int(value, name)
        if type(depth) is not int or depth not in range(1, 7):
            raise ValueError('depth must be an integer from one through six')
        self.depth = self.stages = depth
        self.capacity_control = False
        self.hidden_dim, self.num_heads = hidden_dim, num_heads
        for index in range(1, depth + 1):
            in_features = input_dim if index == 1 else hidden_dim * num_heads
            layer = GATConv(in_features, hidden_dim, heads=num_heads, concat=index < depth)
            nn.init.xavier_uniform_(layer.lin.weight)
            setattr(self, 'gat' + str(index), layer)
        self.GAT_fc = nn.Linear(hidden_dim, GAT_z_dim * 2)
        nn.init.xavier_uniform_(self.GAT_fc.weight)

    def forward(self, x, edge_index):
        h = self.gat1(x, edge_index)
        for index in range(2, self.depth + 1):
            h = getattr(self, 'gat' + str(index))(F.leaky_relu(h), edge_index)
        mean, log_variance = torch.chunk(self.GAT_fc(h), 2, dim=1)
        return mean, torch.exp(log_variance)


class DualBranchReplacement(nn.Module):
    """MLP10+GAT10 or parameter-matched GAT10+GAT10 trainer-compatible model."""
    def __init__(self, model_config, dataloader_dict, family, replacement_width=None):
        super().__init__()
        if family not in ('mlp_gat', 'gat_gat_matched'):
            raise ValueError('Unknown replacement family: ' + str(family))
        if model_config['z_dim'] != 20 or model_config['gat_attention_heads'] != 4:
            raise ValueError('This study fixes decoder20 and four graph attention heads')
        if model_config.get('gat_depth', 3) != 3 or model_config['gat_encoder_hidden'] != 512:
            raise ValueError('The original branch B must remain depth3/width512')
        if model_config.get('gat_capacity_control', False):
            raise ValueError('Nodewise controls are not part of these replacements')
        if model_config['beta2'] != 1e-4:
            raise ValueError('The frozen primary Gaussian-KL coefficient is beta2=1e-4')
        self.family = family
        self.input_dim = model_config['input_dim']
        self.dtype, self.device = model_config['dtype'], model_config['device']
        self.timepoints = tuple(model_config['timepoints'])
        self.replacement_z_dim, self.GAT_z_dim, self.SVGP_z_dim = 10, 10, 0
        self.replacement_width = replacement_width
        if family == 'mlp_gat':
            self.ReplacementEncoder = GaussianMLPEncoder(
                input_dim=self.input_dim, SVGP_z_dim=10,
                hidden_dims=model_config['svgp_encoder_layers']).to(dtype=self.dtype)
        else:
            _positive_int(replacement_width, 'replacement_width')
            self.ReplacementEncoder = FlexibleGATEncoder(
                self.input_dim, 10, hidden_dim=replacement_width, num_heads=4, depth=3).to(dtype=self.dtype)
        self.GATEncoder = FlexibleGATEncoder(
            self.input_dim, 10, hidden_dim=512, num_heads=4, depth=3).to(dtype=self.dtype)
        self.decoder = Decoder(input_dim=self.input_dim, z_dim=20,
                               decoder_layers=model_config['decoder_layers']).to(dtype=self.dtype)
        self.gammas, self.kmeans_center_dict = {}, {}
        self.kmeans_cluster_dict, self.kmeans_index_dict = {}, {}
        # Diagnostics are plain Python data, never parameters/buffers/state_dict.
        self._branch_epoch_diagnostics = None
        self.last_branch_losses = None

    def _branch_A(self, y, edge_index, seed_count=None):
        if self.family == 'mlp_gat':
            # Preserve source GP-network training BatchNorm on SEED nodes only.
            return self.ReplacementEncoder(y if seed_count is None else y[:seed_count])
        mean, variance = self.ReplacementEncoder(y, edge_index)
        return (mean, variance) if seed_count is None else (mean[:seed_count], variance[:seed_count])

    @staticmethod
    def _gaussian_kl(mean, variance):
        return -0.5 * torch.sum(1 + torch.log(variance) - mean.pow(2) - variance) / 10

    def _record_branch_losses(self, tp, kl_a, kl_b):
        if not self.training:
            return
        a, b = float(kl_a.detach().cpu().item()), float(kl_b.detach().cpu().item())
        self.last_branch_losses = {'KL_A': a, 'KL_B': b, 'timepoint': str(tp)}
        if self._branch_epoch_diagnostics is None:
            self._branch_epoch_diagnostics = {'KL_A': 0., 'KL_B': 0., 'batch_count': 0,
                'per_timepoint': {}, 'KL_A_epoch_loss': 0., 'KL_B_epoch_loss': 0.}
        record = self._branch_epoch_diagnostics
        record['KL_A'] += a
        record['KL_B'] += b
        record['batch_count'] += 1
        part = record['per_timepoint'].setdefault(str(tp), {'KL_A': 0., 'KL_B': 0., 'batch_count': 0})
        part['KL_A'] += a
        part['KL_B'] += b
        part['batch_count'] += 1
        # Source losses SUM the per-timepoint batch means. Counts differ by TP.
        for name in ('KL_A', 'KL_B'):
            record[name + '_epoch_loss'] = sum(p[name] / p['batch_count'] for p in record['per_timepoint'].values())

    def forward(self, x, y, edge_index, tp, batch_size):
        if type(batch_size) is not int or not 0 < batch_size <= len(y):
            raise ValueError('batch_size must identify a nonempty seed-node prefix')
        if tp not in self.timepoints:
            raise KeyError('Unknown timepoint: ' + str(tp))
        mean_a, var_a = self._branch_A(y, edge_index, batch_size)
        # Keep original sampling order: branch A noise, then branch B noise.
        sample_a = mean_a + torch.randn_like(mean_a) * torch.sqrt(var_a)
        mean_b, var_b = self.GATEncoder(y, edge_index)
        mean_b, var_b = mean_b[:batch_size], var_b[:batch_size]
        sample_b = mean_b + torch.randn_like(mean_b) * torch.sqrt(var_b)
        latent = torch.cat([sample_a, sample_b], dim=1)
        reconstruction = torch.sum((y[:batch_size] - self.decoder(latent)) ** 2) / self.input_dim
        kl_a, kl_b = self._gaussian_kl(mean_a, var_a), self._gaussian_kl(mean_b, var_b)
        alignment = F.mse_loss(sample_a.norm(dim=1) / 10, sample_b.norm(dim=1) / 10, reduction='sum')
        self._record_branch_losses(tp, kl_a, kl_b)
        return reconstruction, reconstruction.new_zeros(()), kl_a + kl_b, alignment, latent

    def all_latent_samples(self, X, Y, edge_index, tp):
        if tp not in self.timepoints:
            raise KeyError('Unknown timepoint: ' + str(tp))
        Y = torch.as_tensor(Y, dtype=self.dtype, device=self.device)
        edge_index = torch.as_tensor(edge_index, dtype=torch.long, device=self.device)
        mean_a, _ = self._branch_A(Y, edge_index)
        mean_b, _ = self.GATEncoder(Y, edge_index)
        return torch.cat([mean_a, mean_b], dim=1).detach().cpu().numpy()


def make_model_factory(original_class, trial, observer=None):
    """Pair shared initial tensors/RNG with the exact original full control."""
    family = trial['family']
    if family not in ('mlp_gat', 'gat_gat_matched', 'dual'):
        raise ValueError('Unsupported family: ' + str(family))
    width = trial.get('replacement_width', trial.get('width')) if family == 'gat_gat_matched' else None
    if family == 'gat_gat_matched':
        _positive_int(width, 'replacement_width')

    def construct(cfg, dl):
        reference_cfg = dict(cfg, gat_depth=3, gat_encoder_hidden=512, gat_capacity_control=False)
        baseline = original_class(reference_cfg, dl)
        cpu_rng = torch.get_rng_state().clone()
        device = torch.device(cfg['device'])
        cuda_rng = torch.cuda.get_rng_state(device).clone() if device.type == 'cuda' else None
        try:
            reference_count = _count(baseline)
            ref_decoder_sha = state_digest(baseline.decoder.state_dict())
            ref_gat_sha = state_digest(baseline.GATEncoder.state_dict())
            ref_mlp_sha = state_digest(baseline.SVGPEncoder.state_dict())
            gp_count = _count(baseline.svgp_dict)
            assert gp_count == 0, 'GP module has trainable parameters; exact MLP matching needs re-audit'
            expected_mlp = mlp_parameter_count(cfg['input_dim'], cfg['svgp_encoder_layers'], 10)
            assert _count(baseline.SVGPEncoder) == expected_mlp
            assert _count(baseline.GATEncoder) == parameter_count(cfg['input_dim'], 512, 3, 10, 4)
            if family == 'dual':
                model = baseline
                replacement_count = expected_mlp
                replacement_sha = ref_mlp_sha
            else:
                model = DualBranchReplacement(reference_cfg, dl, family=family, replacement_width=width)
                model.GATEncoder.load_state_dict(baseline.GATEncoder.state_dict(), strict=True)
                model.decoder.load_state_dict(baseline.decoder.state_dict(), strict=True)
                if family == 'mlp_gat':
                    model.ReplacementEncoder.load_state_dict(baseline.SVGPEncoder.state_dict(), strict=True)
                    assert state_digest(model.ReplacementEncoder.state_dict()) == ref_mlp_sha
                assert not hasattr(model, 'svgp_dict') and not hasattr(model, 'SVGPEncoder')
                replacement_count = _count(model.ReplacementEncoder)
                expected_replacement = expected_mlp if family == 'mlp_gat' else parameter_count(cfg['input_dim'], width, 3, 10, 4)
                assert replacement_count == expected_replacement
                replacement_sha = state_digest(model.ReplacementEncoder.state_dict())
            decoder_sha, gat_sha = state_digest(model.decoder.state_dict()), state_digest(model.GATEncoder.state_dict())
            assert decoder_sha == ref_decoder_sha and gat_sha == ref_gat_sha
            count = _count(model)
            assert count == replacement_count + _count(model.GATEncoder) + _count(model.decoder)
            if family == 'mlp_gat':
                assert count == reference_count
            if 'expected_parameters' in trial:
                assert count == trial['expected_parameters']
            relative = (count - reference_count) / reference_count
            if family == 'gat_gat_matched':
                assert abs(relative) < .005, ('Parameter match outside0.5%', relative)
            metadata = {
                'family': family, 'depth': 3, 'width': 512, 'replacement_width': width,
                'latent_dimensions': 20, 'branch_A_latent_dimensions': 10, 'branch_B_latent_dimensions': 10,
                'GAT_latent_dimensions': 10, 'shared_encoder_across_timepoints': True,
                'model_contains_GP_modules': family == 'dual',
                'reference_decoder_sha256': ref_decoder_sha, 'common_decoder_sha256': decoder_sha,
                'reference_GAT_sha256': ref_gat_sha, 'common_GAT_sha256': gat_sha,
                'reference_SVGP_encoder_sha256': ref_mlp_sha, 'replacement_encoder_sha256': replacement_sha,
                'common_weights_sha256': state_digest({'GAT.' + k: v for k, v in model.GATEncoder.state_dict().items()}
                    | {'decoder.' + k: v for k, v in model.decoder.state_dict().items()}),
                'reference_postinit_cpu_rng_sha256': _rng_digest(cpu_rng),
                'reference_postinit_cuda_rng_sha256': None if cuda_rng is None else _rng_digest(cuda_rng),
                'cuda_rng_device': str(device) if cuda_rng is not None else None,
                'active_model_parameters': count, 'reference_dual_parameters': reference_count,
                'active_replacement_parameters': replacement_count, 'reference_SVGP_encoder_parameters': expected_mlp,
                'active_GAT_parameters': _count(model.GATEncoder), 'active_decoder_parameters': _count(model.decoder),
                'GP_module_trainable_parameters': gp_count, 'parameter_count_formula_verified': True,
                'exact_MLP_parameter_match': family == 'mlp_gat' and count == reference_count,
                'relative_parameter_difference': relative,
                'replacement_initialized_from_reference_MLP': family == 'mlp_gat',
                'replacement_GAT_independent_initialization': family == 'gat_gat_matched',
                'shared_weights_are_trainable': True,
                'additional_weights_rng_policy': 'Fresh replacement-GAT draws after full reference construction; restore reference post-init CPU and selected-device CUDA RNG',
                'latent_order': ['GPVAE10', 'GAT10'] if family == 'dual' else
                    ['MLP10' if family == 'mlp_gat' else 'ReplacementGAT10', 'OriginalGAT10'],
                'replacement_objective': None if family == 'dual' else {
                    'SVGP_KL': 0, 'KL_A': '-0.5*sum(1+log(var_A)-mean_A^2-var_A)/10',
                    'KL_B': '-0.5*sum(1+log(var_B)-mean_B^2-var_B)/10',
                    'trainer_GAT_KL': 'KL_A+KL_B', 'Gaussian_KL_coefficient': cfg['beta2'],
                    'cross_encoder_norm_alignment': 'sum MSE(norm(sample_A)/10,norm(sample_B)/10)',
                    'reconstruction': 'sum squared error/input_dim',
                    'temporal_KMeans_and_OT': 'Original trainer terms and coefficients unchanged',
                    'epoch_KL_replay': 'Sum of per-timepoint batch means, not global batch mean'},
            }
            model._replacement_initialization = metadata
            if observer is not None:
                if callable(observer):
                    observer(metadata)
                else:
                    observer['initialization'] = metadata
            return model
        finally:
            torch.set_rng_state(cpu_rng)
            if cuda_rng is not None:
                torch.cuda.set_rng_state(cuda_rng, device)
    return construct


def install_training_adapter(train_utils, trial, observer=None):
    """Rebind a namespace; never replace the original class's module global."""
    original_module = train_utils.SpaDOT
    original_class = original_module.SpaDOT
    if original_class.__init__.__globals__.get('SpaDOT') is not original_class:
        raise AssertionError('Original SpaDOT class global is not intact')
    train_utils.SpaDOT = SimpleNamespace(SpaDOT=make_model_factory(original_class, trial, observer))
    assert original_module.SpaDOT is original_class
    return original_module
