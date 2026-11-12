"""Read-only fixed-weight GP/GAT mechanisms; no trainer or solver import.

The native network and native posterior methods are called without changing
their implementation.  Source expression includes every spot for the GP.
GAT differentiation uses the exact incoming three-hop induced closure only
after its anchor output agrees with native full-graph inference.
"""
from contextlib import contextmanager
import hashlib
import numpy as np
import torch


def array_digest(value):
    if torch.is_tensor(value):
        value = value.detach().cpu().numpy()
    value = np.ascontiguousarray(value)
    h = hashlib.sha256()
    h.update(str(value.dtype).encode())
    h.update(str(value.shape).encode())
    h.update(value.tobytes())
    return h.hexdigest()


def fit_whitener(z, relative_ridge=1e-6):
    """Fit ONCE on unperturbed branch means; never refit after perturbing.

    Ridge is proportional to mean baseline variance.  Whitened and raw maps
    are both returned downstream; the ridge prevents near-zero axes from
    dominating.  This does not imply arbitrary reparameterization invariance.
    """
    z = np.asarray(z, dtype=np.float64)
    assert z.ndim == 2 and z.shape[1] == 10 and np.isfinite(z).all()
    centered = z - z.mean(axis=0, keepdims=True)
    covariance = centered.T @ centered / (len(z) - 1)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    mean_variance = float(np.trace(covariance) / z.shape[1])
    assert mean_variance > 0 and eigenvalues.min() >= -1e-10 * eigenvalues.max()
    ridge = float(relative_ridge * mean_variance)
    used = np.maximum(eigenvalues, 0.) + ridge
    whitener = (eigenvectors * (1. / np.sqrt(used))) @ eigenvectors.T
    metadata = dict(relative_ridge=float(relative_ridge), ridge=ridge,
                    baseline_mean=z.mean(axis=0).tolist(),
                    covariance=covariance.tolist(), eigenvalues=eigenvalues.tolist(),
                    whitener_sha256=array_digest(whitener),
                    fitted_on_unperturbed_branch_only=True,
                    interpretation='Fixed ridge-Mahalanobis output metric; raw sensitivity retained')
    return whitener, metadata


def incoming_closure(edge_index, n_nodes, anchor, depth=3):
    """Exact source-to-target incoming support, including the anchor itself."""
    edge = edge_index.detach().cpu().numpy() if torch.is_tensor(edge_index) else np.asarray(edge_index)
    assert edge.shape[0] == 2 and 0 <= anchor < n_nodes
    distance = np.full(n_nodes, -1, dtype=np.int64)
    distance[anchor] = 0
    frontier = np.array([anchor], dtype=np.int64)
    for hop in range(1, depth + 1):
        sources = np.unique(edge[0, np.isin(edge[1], frontier)])
        frontier = sources[distance[sources] < 0]
        distance[frontier] = hop
        if not len(frontier):
            break
    nodes = np.flatnonzero(distance >= 0)
    return nodes, distance


def all_incoming_hops(edge_index, n_nodes, anchor):
    edge = edge_index.detach().cpu().numpy() if torch.is_tensor(edge_index) else np.asarray(edge_index)
    distance = np.full(n_nodes, -1, dtype=np.int64)
    distance[anchor] = 0
    frontier = np.array([anchor], dtype=np.int64)
    hop = 0
    while len(frontier):
        hop += 1
        sources = np.unique(edge[0, np.isin(edge[1], frontier)])
        frontier = sources[distance[sources] < 0]
        distance[frontier] = hop
    return distance


def induced_edges(edge_index, nodes, n_nodes):
    """Preserve original edge ordering; GATConv handles native self loops."""
    nodes = np.asarray(nodes, dtype=np.int64)
    lookup = torch.full((n_nodes,), -1, dtype=torch.long, device=edge_index.device)
    keep_nodes = torch.as_tensor(nodes, dtype=torch.long, device=edge_index.device)
    lookup[keep_nodes] = torch.arange(len(nodes), device=edge_index.device)
    keep = (lookup[edge_index[0]] >= 0) & (lookup[edge_index[1]] >= 0)
    return lookup[edge_index[:, keep]]


def gp_anchor_mean(model, tp, X, Y, anchors):
    """Native posterior, full expression-derived pseudo-mean AND uncertainty.

    Only the query rows are reduced. index_points_train and its b remain ALL
    sources, preserving N_train/b.  No RBF-prior-as-influence substitution.
    """
    anchors = torch.as_tensor(anchors, dtype=torch.long, device=X.device)
    pseudo_mean, noise = model.SVGPEncoder(Y)
    columns = []
    for dim in range(model.SVGP_z_dim):
        mean, _, _, _ = model.svgp_dict[str(tp)].approximate_posterior_params(
            X[anchors], X, pseudo_mean[:, dim], noise[:, dim])
        columns.append(mean)
    return torch.stack(columns, dim=1)


def gat_anchor_mean(model, Y, edge_index, anchor, nodes=None):
    if nodes is None:
        nodes, _ = incoming_closure(edge_index, len(Y), anchor, depth=3)
    nodes = np.asarray(nodes, dtype=np.int64)
    position = np.flatnonzero(nodes == anchor)
    assert len(position) == 1
    subedge = induced_edges(edge_index, nodes, len(Y))
    selected = torch.as_tensor(nodes, dtype=torch.long, device=Y.device)
    means, _ = model.GATEncoder(Y[selected], subedge)
    return means[int(position[0])], nodes


def _basis_jacobian_energy(output, expression, whitener):
    """Ten exact basis VJPs per output metric, accumulated without storing J.

    Primary raw squared-Frobenius energy is invariant under output rotations
    and normalization removes a uniform branch scale.  Secondary frozen-ridge
    whitening requires ten additional VJPs, not additional network forwards.
    """
    assert tuple(output.shape) == (10,) and expression.requires_grad
    raw_energy = torch.zeros(len(expression),dtype=expression.dtype,device=expression.device)
    white_energy = torch.zeros_like(raw_energy)
    for dim in range(10):
        gradient, = torch.autograd.grad(output[dim], expression,
                                        retain_graph=True, create_graph=False)
        assert torch.isfinite(gradient).all()
        raw_energy += gradient.detach().square().sum(dim=1)
        del gradient
    W = torch.as_tensor(whitener, dtype=expression.dtype, device=expression.device)
    whitened_output = output @ W
    for dim in range(10):
        gradient, = torch.autograd.grad(whitened_output[dim], expression,
                                        retain_graph=dim < 9, create_graph=False)
        assert torch.isfinite(gradient).all()
        white_energy += gradient.detach().square().sum(dim=1)
        del gradient
    result = dict(raw_energy=raw_energy.cpu().numpy(),
                  whitened_energy=white_energy.cpu().numpy(),
                  raw_norm=raw_energy.sqrt().cpu().numpy(),
                  whitened_norm=white_energy.sqrt().cpu().numpy())
    del whitened_output, raw_energy, white_energy
    return result


def exact_influence(model, tp, X, Y, edge_index, anchor, whiteners,
                    baseline_anchor=None, atol=1e-9, rtol=1e-10):
    """Paired fixed-weight input sensitivities for one unlabeled target.

    baseline_anchor, if supplied, is a saved/reproduced 20D row: GP then GAT.
    Full-graph no-grad baseline equality is an external prerequisite.  The
    tight reduced-query/subgraph agreement is separately recorded here.
    """
    assert not model.training and not model.SVGPEncoder.training and not model.GATEncoder.training
    assert X.dtype == Y.dtype == torch.float64
    nodes, support_distance = incoming_closure(edge_index, len(Y), anchor, 3)
    result = {}
    values = {}
    for branch in ('GP', 'GAT'):
        expression = Y.detach().clone().requires_grad_(True)
        if branch == 'GP':
            output = gp_anchor_mean(model, tp, X, expression, [anchor])[0]
        else:
            output, _ = gat_anchor_mean(model, expression, edge_index, anchor, nodes)
        assert torch.isfinite(output).all()
        values[branch] = output.detach().cpu().numpy()
        if baseline_anchor is not None:
            reference = np.asarray(baseline_anchor)[0:10 if branch == 'GP' else 20]
            if branch == 'GAT':
                reference = np.asarray(baseline_anchor)[10:20]
            assert np.allclose(values[branch], reference, atol=atol, rtol=rtol), (
                branch, anchor, float(np.max(np.abs(values[branch] - reference))))
        result[branch] = _basis_jacobian_energy(output, expression, whiteners[branch])
        del output, expression
    assert np.all(result['GAT']['raw_energy'][support_distance < 0] == 0.)
    assert np.all(result['GAT']['whitened_energy'][support_distance < 0] == 0.)
    result['metadata'] = dict(anchor=int(anchor), n_sources=len(Y), n_genes=int(Y.shape[1]),
        exact_output_basis_vectors=10, raw_basis_VJPs=10, whitened_basis_VJPs=10,
        total_basis_VJPs_per_branch=20, stochastic_probes=0,
        primary_metric='raw squared-Frobenius energy, normalized per anchor',
        secondary_metric='fixed-baseline ridge-Mahalanobis squared energy',
        GP_sources='ALL timepoint spots; derivative through pseudo-mean and noise',
        GAT_sources='Exact incoming3-hop induced closure; zeros outside are structural',
        incoming3hop_nodes=int(len(nodes)), support_nodes=nodes.tolist(),
        baseline_anchor_atol=atol, baseline_anchor_rtol=rtol,
        GP_query_and_GAT_closure_mean={key:value.tolist() for key,value in values.items()},
        fixed_model_weights=True, deterministic_inference=True)
    if baseline_anchor is not None:
        result['metadata']['baseline_anchor_max_abs_errors'] = {
            'GP': float(np.max(np.abs(values['GP'] - np.asarray(baseline_anchor)[:10]))),
            'GAT': float(np.max(np.abs(values['GAT'] - np.asarray(baseline_anchor)[10:]))) }
    return result


def summarize_influence(energy, coordinates, hops, anchor, distance_edges=None):
    """Keep self influence visible; compare mass AND per-source influence."""
    energy = np.asarray(energy, dtype=np.float64)
    coordinates, hops = np.asarray(coordinates), np.asarray(hops)
    assert energy.ndim == 1 and len(energy) == len(coordinates) == len(hops)
    assert np.isfinite(energy).all() and np.min(energy) >= 0
    nonself = np.arange(len(energy)) != anchor
    distance = np.linalg.norm(coordinates - coordinates[anchor], axis=1)
    total, nonself_total = float(energy.sum()), float(energy[nonself].sum())
    valid_total = total > 0
    valid_nonself = nonself_total > 0
    mass = energy / total if valid_total else np.zeros_like(energy)
    nonself_mass = np.zeros_like(energy)
    if nonself_total > 0:
        nonself_mass[nonself] = energy[nonself] / nonself_total
    def radius(weights):
        if weights.sum() == 0:
            return None
        order = np.argsort(distance, kind='stable')
        cumulative = np.cumsum(weights[order]) / weights.sum()
        return float(distance[order[min(int(np.searchsorted(cumulative, .9)), len(order)-1)]])
    rows = []
    for label, selected in [('self', ~nonself), ('1', hops == 1), ('2', hops == 2),
                            ('3', hops == 3), ('>3_or_unreachable', (hops > 3) | (hops < 0))]:
        n = int(selected.sum())
        rows.append(dict(bin=label, n_sources=n, total_energy=float(energy[selected].sum()),
                         total_mass=float(mass[selected].sum()) if valid_total else None,
                         mean_mass_per_source=float(mass[selected].mean()) if n and valid_total else None,
                         nonself_mass=float(nonself_mass[selected].sum()) if valid_nonself else None,
                         mean_nonself_mass_per_source=float(nonself_mass[selected].mean()) if n and valid_nonself else None))
    if distance_edges is None:
        distance_edges = np.linspace(0, max(float(distance.max()), 1e-12) * (1 + 1e-12), 16)
    distance_rows = []
    for lo, hi in zip(distance_edges[:-1], distance_edges[1:]):
        selected = nonself & (distance >= lo) & (distance < hi)
        n = int(selected.sum())
        distance_rows.append(dict(lower=float(lo), upper=float(hi), n_sources=n,
            total_mass=float(mass[selected].sum()) if valid_total else None,
            nonself_mass=float(nonself_mass[selected].sum()) if valid_nonself else None,
            mean_mass_per_source=float(mass[selected].mean()) if n and valid_total else None,
            mean_nonself_mass_per_source=float(nonself_mass[selected].mean()) if n and valid_nonself else None))
    return dict(self_mass=float(mass[anchor]) if valid_total else None, radius90=radius(mass),
                nonself_radius90=radius(nonself_mass), total_raw_energy=total,
                nonself_raw_energy=nonself_total,
                beyond3hop_mass=float(mass[(hops > 3) | (hops < 0)].sum()) if valid_total else None,
                beyond3hop_nonself_mass=float(nonself_mass[(hops > 3) | (hops < 0)].sum()) if valid_nonself else None,
                degenerate_zero_total=not valid_total, degenerate_zero_nonself=not valid_nonself,
                hop_rows=rows, distance_rows=distance_rows,
                normalized_mass=mass, normalized_nonself_mass=nonself_mass)


@contextmanager
def temporary_inducing_coordinates(svgp, coordinates):
    original = svgp.inducing_index_points
    before = array_digest(original)
    replacement = torch.as_tensor(coordinates, dtype=original.dtype, device=original.device)
    try:
        svgp.inducing_index_points = replacement
        yield
    finally:
        svgp.inducing_index_points = original
        assert array_digest(svgp.inducing_index_points) == before


def _drift(current, baseline, whitener):
    difference = np.asarray(current) - np.asarray(baseline)
    white = difference @ np.asarray(whitener)
    return dict(rms_raw=float(np.sqrt(np.mean(difference ** 2))),
                median_spot_L2_raw=float(np.median(np.linalg.norm(difference, axis=1))),
                rms_whitened=float(np.sqrt(np.mean(white ** 2))),
                median_spot_L2_whitened=float(np.median(np.linalg.norm(white, axis=1))),
                max_abs_raw=float(np.max(np.abs(difference))))


def coordinate_geometry_experiment(model, tp, X, Y, edge_index, raw_coordinates,
                                  scaler_mean, scaler_scale, graph_builder,
                                  whiteners, baseline=None, atol=1e-9, rtol=1e-10):
    """Joint geometry counterfactual, not frozen-inducing query perturbation.

    Apply the same raw-coordinate transformation to observed AND inducing
    locations, then map through the ORIGINAL fixed StandardScaler.  Powers of
    two preserve distances' ordering mathematically; actual native graph edge
    hashes must still match.  Nonmatching transforms are explicitly excluded.
    Reflection is a distance-preserving control.  Original tensors are restored
    in finally.  No scale/whitener refit and no learned parameter changes.
    """
    assert not model.training
    raw_coordinates = np.asarray(raw_coordinates, dtype=np.float64)
    center, scale = np.asarray(scaler_mean), np.asarray(scaler_scale)
    original_U = model.svgp_dict[str(tp)].inducing_index_points
    original_U_hash = array_digest(original_U)
    U = original_U.detach().cpu().numpy()
    assert np.allclose((raw_coordinates-center)/scale, X.detach().cpu().numpy(), atol=1e-12, rtol=1e-12)
    raw_U = U * scale + center
    reference_edge = array_digest(edge_index)
    assert array_digest(graph_builder(raw_coordinates)) == reference_edge
    anchors = np.arange(len(Y))
    with torch.no_grad():
        base_gp = gp_anchor_mean(model, tp, X, Y, anchors).cpu().numpy()
        base_gat = model.GATEncoder(Y, edge_index)[0].cpu().numpy()
    if baseline is not None:
        assert np.allclose(base_gp, np.asarray(baseline)[:, :10], atol=atol, rtol=rtol)
        assert np.array_equal(base_gat, np.asarray(baseline)[:, 10:])
    rows, arrays = [], {}
    transforms = [('identity', 1., False), ('scale_0.5', .5, False),
                  ('scale_2', 2., False), ('rigid_reflection', 1., True)]
    for name, multiplier, reflected in transforms:
        coefficient = -multiplier if reflected else multiplier
        raw_new = coefficient * raw_coordinates
        raw_new_U = coefficient * raw_U
        transformed_X = (raw_new-center)/scale
        transformed_U = (raw_new_U-center)/scale
        rebuilt = graph_builder(raw_new)
        current_edge_hash = array_digest(rebuilt)
        graph_exact = current_edge_hash == reference_edge
        row = dict(transform=name, raw_coordinate_multiplier=coefficient,
            baseline_graph_sha256=reference_edge, rebuilt_graph_sha256=current_edge_hash,
            rebuilt_graph_exact=graph_exact, fixed_original_scaler=True,
            joint_observed_and_inducing_transform=True, fixed_kernel_scale=True,
            learned_parameters_changed=False,
            status='evaluated' if graph_exact else 'excluded_graph_changed',
            interpretation='Geometry response, not evidence of predictive usefulness')
        if graph_exact:
            changed = torch.as_tensor(transformed_X, dtype=X.dtype, device=X.device)
            with temporary_inducing_coordinates(model.svgp_dict[str(tp)], transformed_U), torch.no_grad():
                gp = gp_anchor_mean(model, tp, changed, Y, anchors).cpu().numpy()
                gat = model.GATEncoder(Y, edge_index)[0].cpu().numpy()
            assert np.array_equal(gat, base_gat), 'Same graph/features changed deterministic GAT means'
            row.update(GP=_drift(gp, base_gp, whiteners['GP']),
                       GAT=_drift(gat, base_gat, whiteners['GAT']))
            if name in ('identity', 'rigid_reflection'):
                assert np.allclose(gp, base_gp, atol=atol, rtol=rtol), (name, row['GP'])
            arrays[name] = dict(GP=gp, GAT=gat, coordinates=transformed_X)
        rows.append(row)
    assert array_digest(model.svgp_dict[str(tp)].inducing_index_points) == original_U_hash
    return dict(rows=rows, arrays=arrays, baseline_GP=base_gp, baseline_GAT=base_gat,
                inducing_restored_exactly=True, baseline_graph_sha256=reference_edge,
                perturbations_prespecified=True, no_retraining=True)


def coordinate_api(model, tp, X, Y, edge_index, raw_coordinates, scaler_mean,
                   scaler_scale, edge_rebuilder, whiteners, baseline_latent=None,
                   atol=1e-9, rtol=1e-10):
    """Runner-facing alias; arrays are separate from serializable row metadata."""
    return coordinate_geometry_experiment(model,tp,X,Y,edge_index,raw_coordinates,
        scaler_mean,scaler_scale,edge_rebuilder,whiteners,baseline=baseline_latent,
        atol=atol,rtol=rtol)
