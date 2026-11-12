import hashlib
import json
import os
from pathlib import Path

ROOT = Path(os.environ['SPADOT_LOSS_STUDY_ROOT']).resolve()
METHODS = [
    dict(id='vae', model='SpaDOT_VAE', family='vae', weights=[0., 0., 0.],
         kl='positive Gaussian KL multiplied by cyclic beta1'),
    dict(id='svgp_only', model='SpaDOT_SVGPVAE', family='svgp', weights=[0., 0., 0.],
         kl='negative historical SVGP term subtracted with cyclic beta1'),
    dict(id='graphsage_only', model='SpaDOT_GCNVAE', family='graph', weights=[0., 0., 0.],
         kl='positive Gaussian KL multiplied by constant beta2'),
    dict(id='gat_only', model='SpaDOT_GATVAE', family='graph', weights=[0., 0., 0.],
         kl='positive Gaussian KL multiplied by constant beta2'),
    dict(id='dual_without_penalties', model='SpaDOT_SVGPGATVAE', family='dual', weights=[0., 0., 0.],
         kl='negative SVGP term with cyclic beta1; Gaussian GAT KL with beta2'),
    dict(id='dual_alignment_only', model='SpaDOT_SVGPGATVAE_alignment', family='alignment', weights=[.1, 0., 0.],
         kl='negative SVGP term with cyclic beta1; Gaussian GAT KL with beta2'),
    dict(id='without_ot', model='SpaDOT_SVGPGATVAE_alignment', family='alignment', weights=[.1, .1, 0.],
         kl='negative SVGP term with cyclic beta1; Gaussian GAT KL with beta2'),
]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path); tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, default=str) + '\n'); tmp.replace(path)


def check_frozen(manifest):
    assert all(sha(ROOT / p) == h for p, h in manifest['frozen_hashes'].items())
    assert all(sha(p) == h for p, h in manifest['original_historical_hashes'].items())
    assert all(sha(v['input']) == v['input_sha256'] for v in manifest['datasets'].values())
    assert all(sha(p) == h for p, h in manifest['parent_frozen_hashes'].items())
