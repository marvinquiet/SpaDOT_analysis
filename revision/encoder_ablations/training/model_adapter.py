"""Model-only routing for the native-sampler encoder ablation study.

The frozen native SpaDOT module is used directly for dual controls. Replacements
are separate modules and never modify that class's module globals or the loader.
"""
from types import SimpleNamespace

from standalone_adapter import (
    FlexibleGATEncoder, StandaloneGATVAE, parameter_count, state_digest,
    make_model_factory as standalone_factory,
)
from replacement_adapter import (
    DualBranchReplacement, mlp_parameter_count,
    make_model_factory as replacement_factory,
)


def make_model_factory(original_class, trial, observer=None):
    family = trial['family']
    if family == 'gat_only':
        factory = standalone_factory(original_class, trial, observer)
    elif family in ('dual', 'mlp_gat', 'gat_gat_matched'):
        factory = replacement_factory(original_class, trial, observer)
    else:
        raise ValueError('Unknown prespecified model family: ' + str(family))

    def construct(config, dataloaders):
        # The graph-wide arm changes these two configuration values through the
        # unchanged native prepare_dataloader; no model or sampler override does.
        if (config['knn_cutoff'], config['max_neighbors']) not in ((6, 30), (12, 60)):
            raise AssertionError('Unexpected spatial-neighborhood configuration')
        model = factory(config, dataloaders)
        metadata = model._replacement_initialization
        metadata['model_factory_changes_sampler'] = False
        metadata['expected_native_sampling_hops'] = 2
        metadata['sampling_fanout_rule'] = 'max(30, knn_cutoff * round(1/1000 * n_spots))'
        metadata['graph_knn_cutoff'] = config['knn_cutoff']
        metadata['graph_max_neighbors'] = config['max_neighbors']
        if isinstance(observer, dict) and observer.get('capture_batch_hashes'):
            from runtime import _array_digest
            import weakref
            import inspect
            batches = observer.setdefault('first_training_batches', {})
            original_forward = model.forward.__func__
            signature = inspect.signature(original_forward)
            weak_model = weakref.ref(model)
            def observed_forward(*values, **keywords):
                module = weak_model()
                if module is None:
                    raise RuntimeError('Forward observer outlived its model')
                if module.training:
                    bound = signature.bind(module, *values, **keywords)
                    x, y, edges, tp, seed_count = list(bound.arguments.values())[1:]
                    if str(tp) not in batches:
                        batches[str(tp)] = {
                            'seed_count': int(seed_count), 'sampled_nodes': int(len(y)),
                            'sampled_edges': int(edges.shape[1]),
                            'coordinates_sha256': _array_digest(x),
                            'expression_sha256': _array_digest(y),
                            'edge_index_sha256': _array_digest(edges),
                            'seed_coordinates_sha256': _array_digest(x[:seed_count]),
                            'no_extra_loader_iteration_or_random_draw': True,
                        }
                return original_forward(module, *values, **keywords)
            # Native source calls forward directly, bypassing Module pre-hooks.
            # A weak reference avoids retaining the trained GPU model in a cycle.
            observed_forward.__signature__ = signature.replace(parameters=list(signature.parameters.values())[1:])
            model.forward = observed_forward
        return model

    return construct


def install_training_adapter(train_utils, trial, observer=None):
    """Replace only the trainer's constructor namespace; return it for cleanup."""
    original_module = train_utils.SpaDOT
    original_class = original_module.SpaDOT
    if original_class.__init__.__globals__.get('SpaDOT') is not original_class:
        raise AssertionError('Native SpaDOT class namespace is not intact')
    train_utils.SpaDOT = SimpleNamespace(
        SpaDOT=make_model_factory(original_class, trial, observer))
    if original_module.SpaDOT is not original_class:
        raise AssertionError('Native SpaDOT module global was modified')
    return original_module
