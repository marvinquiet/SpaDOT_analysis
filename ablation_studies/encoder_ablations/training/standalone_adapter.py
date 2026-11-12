"""Isolated GAT-VAE replacement; the frozen SpaDOT source is never modified.

The standalone model has one shared 20D GAT encoder and the original 20D
decoder. It contains no GP modules. Its zero ``alignment`` output is the
removed cross-encoder norm-matching term, not the temporal OT objective.
"""

import hashlib
from types import SimpleNamespace

import torch
from torch import nn
from torch.nn import functional as F
from torch_geometric.nn import GATConv

from SpaDOT.model.decoder import Decoder


def _positive_int(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(name + " must be a positive integer")
    return value


def parameter_count(input_dim, width, depth, GAT_z_dim=20, num_heads=4):
    """Exact GAT-encoder parameter count, including its mean/logvar head.

    This excludes the shared decoder and assumes the reference GATConv options:
    no residual projection or edge features, concatenated intermediate heads,
    averaged final heads, and a bias on every graph layer.
    """
    for value, name in ((input_dim, "input_dim"), (width, "width"),
                        (depth, "depth"), (GAT_z_dim, "GAT_z_dim"),
                        (num_heads, "num_heads")):
        _positive_int(value, name)
    if depth > 6:
        raise ValueError("depth must be at most 6")
    heads_width = num_heads * width
    return (input_dim * heads_width + (depth - 1) * heads_width ** 2
            + depth * 2 * heads_width + (depth - 1) * heads_width + width
            + 2 * GAT_z_dim * width + 2 * GAT_z_dim)


def state_digest(state):
    """Digest named tensor values, including their shapes and dtypes."""
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


class FlexibleGATEncoder(nn.Module):
    """Original GAT stages extended to depths 1..6, without nodewise controls.

    Depth three preserves the native construction order, tensor names,
    initialization calls and forward operations exactly. The other depths are
    isolated experimental extensions. Intermediate layers concatenate heads;
    the final graph layer averages its heads.
    """

    def __init__(self, input_dim, GAT_z_dim, hidden_dim=512, num_heads=4, depth=3):
        super().__init__()
        _positive_int(input_dim, "input_dim")
        _positive_int(GAT_z_dim, "GAT_z_dim")
        _positive_int(hidden_dim, "hidden_dim")
        _positive_int(num_heads, "num_heads")
        if type(depth) is not int or depth not in range(1, 7):
            raise ValueError("depth must be an integer from 1 through 6")
        self.depth = self.stages = depth
        self.capacity_control = False
        self.hidden_dim, self.num_heads = hidden_dim, num_heads
        for index in range(1, depth + 1):
            in_features = input_dim if index == 1 else hidden_dim * num_heads
            layer = GATConv(in_features, hidden_dim, heads=num_heads,
                            concat=index < depth)
            nn.init.xavier_uniform_(layer.lin.weight)
            setattr(self, "gat" + str(index), layer)
        self.GAT_fc = nn.Linear(hidden_dim, GAT_z_dim * 2)
        nn.init.xavier_uniform_(self.GAT_fc.weight)

    def forward(self, x, edge_index):
        h = self.gat1(x, edge_index)
        for index in range(2, self.depth + 1):
            h = getattr(self, "gat" + str(index))(F.leaky_relu(h), edge_index)
        mean, log_variance = torch.chunk(self.GAT_fc(h), 2, dim=1)
        return mean, torch.exp(log_variance)


class StandaloneGATVAE(nn.Module):
    """20D GAT-VAE with the reference trainer's forward/latent API."""

    def __init__(self, model_config, dataloader_dict, depth=None, width=None):
        super().__init__()
        if model_config["z_dim"] != 20:
            raise ValueError("This replacement study requires a 20D latent space")
        if model_config["gat_attention_heads"] != 4:
            raise ValueError("This replacement study fixes four GAT attention heads")
        if model_config.get("gat_capacity_control", False):
            raise ValueError("Standalone replacement does not use nodewise controls")
        self.input_dim = model_config["input_dim"]
        self.GAT_z_dim, self.SVGP_z_dim = 20, 0
        self.dtype, self.device = model_config["dtype"], model_config["device"]
        self.timepoints = tuple(model_config["timepoints"])
        self.GATEncoder = FlexibleGATEncoder(
            input_dim=self.input_dim, GAT_z_dim=20,
            hidden_dim=model_config["gat_encoder_hidden"] if width is None else width,
            num_heads=4,
            depth=model_config.get("gat_depth", 3) if depth is None else depth,
        ).to(dtype=self.dtype)
        self.decoder = Decoder(
            input_dim=self.input_dim, z_dim=20,
            decoder_layers=model_config["decoder_layers"],
        ).to(dtype=self.dtype)
        # The unchanged temporal KMeans/OT trainer populates these dictionaries.
        self.gammas = {}
        self.kmeans_center_dict = {}
        self.kmeans_cluster_dict = {}
        self.kmeans_index_dict = {}

    def forward(self, x, y, edge_index, tp, batch_size):
        if type(batch_size) is not int or not 0 < batch_size <= len(y):
            raise ValueError("batch_size must identify a nonempty seed-node prefix")
        if tp not in self.timepoints:
            raise KeyError("Unknown timepoint: " + str(tp))
        mean, variance = self.GATEncoder(y, edge_index)
        mean, variance = mean[:batch_size], variance[:batch_size]
        latent = mean + torch.randn_like(mean) * torch.sqrt(variance)
        reconstruction = torch.sum((y[:batch_size] - self.decoder(latent)) ** 2) / self.input_dim
        gat_kl = -0.5 * torch.sum(
            1 + torch.log(variance) - mean.pow(2) - variance
        ) / self.GAT_z_dim
        zero = reconstruction.new_zeros(())
        return reconstruction, zero, gat_kl, zero, latent

    def all_latent_samples(self, X, Y, edge_index, tp):
        if tp not in self.timepoints:
            raise KeyError("Unknown timepoint: " + str(tp))
        Y = torch.as_tensor(Y, dtype=self.dtype, device=self.device)
        edge_index = torch.as_tensor(edge_index, dtype=torch.long, device=self.device)
        mean, _ = self.GATEncoder(Y, edge_index)
        return mean.detach().cpu().numpy()


def _copy_standalone_initialization(model, baseline):
    """Copy exact shared tensors and preserve the original decoder coordinates."""
    model.decoder.load_state_dict(baseline.decoder.state_dict(), strict=True)
    target = model.GATEncoder.state_dict()
    source = baseline.GATEncoder.state_dict()
    copied = []
    for name, value in source.items():
        if name.startswith("GAT_fc."):
            continue
        if name in target and target[name].shape == value.shape:
            target[name].copy_(value)
            copied.append(name)
    # The dual latent's original GAT dimensions occupy coordinates 10..19.
    # The new variance logits occupy 20..39; the matching original GAT logits
    # therefore occupy 30..39. Additional coordinates retain fresh draws.
    head_rows = []
    for name in ("GAT_fc.weight", "GAT_fc.bias"):
        src, dst = source[name], target[name]
        if src.shape[0] != 20 or dst.shape[0] != 40:
            raise AssertionError("Unexpected reference/replacement head dimensions")
        if src.shape[1:] == dst.shape[1:]:
            dst[10:20].copy_(src[:10])
            dst[30:40].copy_(src[10:20])
            head_rows.append(name)
    model.GATEncoder.load_state_dict(target, strict=True)
    return copied, head_rows


def make_model_factory(original_class, trial, observer=None):
    """Return a trainer constructor with a paired, original-D3 initialization.

    Only the selected CUDA device's RNG state is observed/restored. Additional
    standalone tensors use fresh construction draws; common tensors are copied
    exactly, not merely assigned equal seeds. No GP object is attached to the
    returned standalone model.
    """
    family = trial["family"]
    if family not in ("dual", "gat_only"):
        raise ValueError("family must be 'dual' or 'gat_only'")
    depth = trial.get("depth", 3)
    width = trial.get("width", 512)
    if family == "dual" and (depth != 3 or width != 512):
        raise ValueError("The dual reference is the original depth3/width512 model")

    def construct(cfg, dl):
        reference_cfg = dict(cfg, gat_depth=3, gat_encoder_hidden=512,
                             gat_capacity_control=False)
        baseline = original_class(reference_cfg, dl)
        cpu_rng = torch.get_rng_state().clone()
        device = torch.device(cfg["device"])
        cuda_rng = torch.cuda.get_rng_state(device).clone() if device.type == "cuda" else None
        reference_decoder_sha = state_digest(baseline.decoder.state_dict())
        copied, head_rows = [], []
        try:
            if family == "dual":
                model = baseline
            else:
                standalone_cfg = dict(cfg, gat_depth=depth, gat_encoder_hidden=width,
                                      gat_capacity_control=False)
                model = StandaloneGATVAE(standalone_cfg, dl)
                with torch.no_grad():
                    copied, head_rows = _copy_standalone_initialization(model, baseline)
                if any(name.startswith(("SVGPEncoder.", "svgp_dict."))
                       for name in model.state_dict()):
                    raise AssertionError("Standalone model unexpectedly contains GP state")
            decoder_sha = state_digest(model.decoder.state_dict())
            if decoder_sha != reference_decoder_sha:
                raise AssertionError("Paired decoder initialization differs")
            active_gat_parameters = sum(p.numel() for p in model.GATEncoder.parameters()
                                        if p.requires_grad)
            expected_gat_parameters = parameter_count(
                cfg["input_dim"], width, depth, 10 if family == "dual" else 20,
                cfg["gat_attention_heads"],
            )
            if active_gat_parameters != expected_gat_parameters:
                raise AssertionError("GAT parameter count disagrees with its architecture")
            metadata = {
                "family": family, "depth": depth, "width": width,
                "latent_dimensions": 20,
                "GAT_latent_dimensions": 10 if family == "dual" else 20,
                "reference_decoder_sha256": reference_decoder_sha,
                "common_decoder_sha256": decoder_sha,
                "common_weights_sha256": decoder_sha,
                "matching_shape_GAT_tensors_copied": copied,
                "expanded_head_tensors_with_copied_rows": head_rows,
                "expanded_head_mapping": None if family == "dual" else {
                    "reference_mean_rows": [0, 10], "replacement_mean_rows": [10, 20],
                    "reference_logvar_rows": [10, 20], "replacement_logvar_rows": [30, 40],
                    "index_convention": "zero-based half-open",
                    "new_mean_rows": [0, 10], "new_logvar_rows": [20, 30],
                    "weight_rows_require_matching_input_width": True,
                },
                "reference_postinit_cpu_rng_sha256": _rng_digest(cpu_rng),
                "reference_postinit_cuda_rng_sha256": None if cuda_rng is None else _rng_digest(cuda_rng),
                "cuda_rng_device": str(device) if cuda_rng is not None else None,
                "active_model_parameters": sum(p.numel() for p in model.parameters() if p.requires_grad),
                "active_GAT_parameters": active_gat_parameters,
                "parameter_count_formula_verified": True,
                "active_decoder_parameters": sum(p.numel() for p in model.decoder.parameters() if p.requires_grad),
                "reference_dual_parameters": sum(p.numel() for p in baseline.parameters() if p.requires_grad),
                "reference_GAT_parameters": sum(p.numel() for p in baseline.GATEncoder.parameters() if p.requires_grad),
                "shared_encoder_across_timepoints": True,
                "model_contains_GP_modules": family == "dual",
                "standalone_objective": None if family == "dual" else {
                    "reconstruction": "sum squared error / input_dim",
                    "GAT_KL": "-0.5 * sum(1 + log(var) - mean^2 - var) / 20",
                    "SVGP_KL": 0, "cross_encoder_norm_alignment": 0,
                    "temporal_KMeans_and_OT": "unchanged trainer terms and coefficients",
                },
            }
            model._replacement_initialization = metadata
            if observer is not None:
                if callable(observer):
                    observer(metadata)
                else:
                    observer["initialization"] = metadata
            return model
        finally:
            torch.set_rng_state(cpu_rng)
            if cuda_rng is not None:
                torch.cuda.set_rng_state(cuda_rng, device)

    return construct


def install_training_adapter(train_utils, trial, observer=None):
    """Rebind the trainer namespace, never the original class's module globals.

    The caller must restore the returned original namespace after its run.
    """
    original_module = train_utils.SpaDOT
    original_class = original_module.SpaDOT
    if original_class.__init__.__globals__.get("SpaDOT") is not original_class:
        raise AssertionError("Original SpaDOT class global is not intact")
    train_utils.SpaDOT = SimpleNamespace(
        SpaDOT=make_model_factory(original_class, trial, observer)
    )
    if original_module.SpaDOT is not original_class:
        raise AssertionError("Original class namespace was modified")
    return original_module
