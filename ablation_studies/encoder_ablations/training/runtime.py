"""Model-only installation and observation of the unchanged native sampler.

No loader factory or prepare_dataloader function is replaced here. Verification
reads metadata and tensors already constructed by native source; it never draws
a batch, constructs a graph, calls a neural model, or consumes random numbers.
"""
from contextlib import contextmanager
from pathlib import Path
import hashlib
import json


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def _array_digest(value):
    import numpy as np
    if hasattr(value, 'detach'):
        value = value.detach().cpu().contiguous().numpy()
    array = np.ascontiguousarray(value)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode())
    digest.update(str(array.shape).encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def _neighbor_values(sampler):
    specification = sampler.num_neighbors
    if callable(getattr(specification, 'get_values', None)):
        values = specification.get_values()
    elif isinstance(specification, (list, tuple)):
        values = specification
    else:
        values = specification.values
    if isinstance(values, dict) or not isinstance(values, (list, tuple)):
        raise AssertionError('Expected homogeneous finite neighbor metadata')
    if any(isinstance(value, bool) or int(value) != value for value in values):
        raise AssertionError('Neighbor counts must be integers')
    return [int(value) for value in values]


def verify_native_loaders(dataloaders, config):
    """Return per-timepoint specs and attach them without iterating a loader.

    graph_k is the native graph's requested cutoff. Finite fanout is a cap on
    incoming sampling and does not prove complete closure around every seed.
    graph_sha256 hashes the actual source edge_index (dtype, shape, bytes).
    """
    from torch_geometric.loader import NeighborLoader
    cutoff, maximum = config['knn_cutoff'], config['max_neighbors']
    if (cutoff, maximum) not in ((6, 30), (12, 60)):
        raise AssertionError('Unexpected graph rule')
    if config['batch_size'] != 512:
        raise AssertionError('Native training batch size must remain512')
    expected_timepoints = list(config['timepoints'])
    if list(dataloaders['dataloaders']) != expected_timepoints:
        raise AssertionError('Loader timepoint order differs from the source configuration')
    results = {}
    for timepoint in expected_timepoints:
        loader = dataloaders['dataloaders'][timepoint]
        if type(loader) is not NeighborLoader:
            raise AssertionError('NeighborLoader was wrapped or replaced')
        graph = loader.data
        n_nodes = int(graph.num_nodes)
        if n_nodes != int(dataloaders['N_train'][timepoint]):
            raise AssertionError('Graph and source training spot counts differ')
        scale = round(1 / 1000 * n_nodes)
        fanout = max(30, cutoff * scale)
        graph_k = min(maximum, cutoff * scale)
        sampler = loader.node_sampler
        actual = _neighbor_values(sampler)
        if actual != [fanout, fanout]:
            raise AssertionError(('Native two-hop finite fanout changed', actual, fanout))
        subgraph = getattr(sampler.subgraph_type, 'value', sampler.subgraph_type)
        if subgraph != 'induced':
            raise AssertionError('Native induced subgraph setting changed')
        if loader.batch_size != config['batch_size']:
            raise AssertionError('Training batch size changed')
        if loader.num_workers != 0 or loader.drop_last:
            raise AssertionError('Native loader worker/drop-last defaults changed')
        if type(loader.sampler).__name__ != 'SequentialSampler':
            raise AssertionError('Native unshuffled seed-node sampler changed')
        if graph.data_index.shape[0] != n_nodes or graph.edge_index.shape[0] != 2:
            raise AssertionError('Unexpected native graph metadata shape')
        specification = {
            'sampling_hops': 2, 'num_neighbors': actual,
            'native_fanout': fanout,
            'fanout_rule': 'max(30, knn_cutoff * round(1/1000 * n_spots))',
            'subgraph_type': subgraph, 'batch_size': int(loader.batch_size),
            'num_workers': int(loader.num_workers), 'drop_last': bool(loader.drop_last),
            'seed_sampler': type(loader.sampler).__name__,
            'n_nodes': n_nodes, 'n_edges': int(graph.edge_index.shape[1]),
            'graph_k': graph_k, 'knn_cutoff': cutoff, 'max_neighbors': maximum,
            'graph_sha256': _array_digest(graph.edge_index),
            'data_index_sha256': _array_digest(graph.data_index),
            'no_loader_iteration': True,
            'context_qualification': 'Finite two-hop sampling union across batch seeds; induced edges may provide longer paths within that union.',
        }
        loader._native_sampling_spec = specification
        results[str(timepoint)] = specification
    return results


inspect_dataloaders = verify_native_loaders


@contextmanager
def installed(train_utils, trial, observer=None):
    """Install only the model factory and prove loader/prepare identities intact."""
    from torch_geometric.loader import NeighborLoader
    from model_adapter import install_training_adapter
    original_loader = train_utils.NeighborLoader
    original_prepare = train_utils.prepare_dataloader
    if original_loader is not NeighborLoader:
        raise AssertionError('Source NeighborLoader is already overridden')
    original_module = install_training_adapter(train_utils, trial, observer)
    try:
        if train_utils.NeighborLoader is not original_loader:
            raise AssertionError('Model installation changed NeighborLoader')
        yield
    finally:
        train_utils.SpaDOT = original_module
        if (train_utils.NeighborLoader is not original_loader or
                train_utils.prepare_dataloader is not original_prepare):
            raise AssertionError('Native loader or prepare function changed')
