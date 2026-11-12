"""Historical forward/inference delegation and explicit loss-slot compatibility."""
import importlib
import sys
from pathlib import Path
from common import ROOT
sys.path.insert(0, str(ROOT / 'historical'))


def original_class(method):
    module = importlib.import_module('model.' + method['model'])
    assert Path(module.__file__).resolve() == (ROOT / 'historical/model' / (method['model'] + '.py')).resolve()
    return getattr(module, method['model'])


def raw_forward(model, method, x, y, edge_index, tp, batch_size):
    cls = original_class(method); family = method['family']
    if family == 'vae':
        return cls.forward(model, x=x, y=y, tp=tp)
    if family == 'svgp':
        return cls.forward(model, x=x, y=y, tp=tp, batch_size=batch_size)
    return cls.forward(model, x=x, y=y, edge_index=edge_index, tp=tp, batch_size=batch_size)


def raw_inference(model, method, X, Y, edge_index, tp):
    cls = original_class(method)
    if method['family'] in ['vae', 'svgp']:
        return cls.all_latent_samples(model, X, Y, tp)
    return cls.all_latent_samples(model, X, Y, edge_index, tp)


def compatible_class(method):
    cls = original_class(method); family = method['family']

    class HistoricalLossInterface(cls):
        def __init__(self, config, data):
            super().__init__(config, data)
            self.gammas = {}; self.kmeans_center_dict = {}
            self.kmeans_cluster_dict = {}; self.kmeans_index_dict = {}

        def forward(self, x, y, edge_index, tp, batch_size):
            values = raw_forward(self, method, x, y, edge_index, tp, batch_size)
            recon = values[0]; zero = recon.new_zeros(())
            if family == 'vae':
                # Native objective subtracts beta1*slot1. -KL here restores +beta1*MLP_KL.
                return recon, -values[1], zero, zero, values[2][:batch_size]
            if family == 'svgp':
                return recon, values[1], zero, zero, values[2][:batch_size]
            if family == 'graph':
                return recon, zero, values[1], zero, values[2][:batch_size]
            if family == 'dual':
                return recon, values[1], values[2], zero, values[3]
            return values

        def all_latent_samples(self, X, Y, edge_index, tp):
            return raw_inference(self, method, X, Y, edge_index, tp)

    return HistoricalLossInterface


def historical_objective(values, method, cfg, beta, kmeans):
    r = cfg['lambda1'] * values[0]; family = method['family']
    if family == 'vae':
        return r + beta * values[1]
    if family == 'svgp':
        return r - beta * values[1]
    if family == 'graph':
        return r + cfg['beta2'] * values[1]
    result = r - beta * values[1] + cfg['beta2'] * values[2]
    if family == 'alignment':
        result = result + cfg['omiga1'] * values[3] + cfg['omiga2'] * kmeans
    return result
