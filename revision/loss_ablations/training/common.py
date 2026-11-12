import hashlib
import json
import os
from pathlib import Path

ROOT = Path(os.environ['SPADOT_LOSS_STUDY_ROOT']).resolve()
COMMIT = '04253a2a5275c478a2fd34f097cee369cdcd477b'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(value, indent=2, default=str) + '\n')
    temp.replace(path)


def check_frozen(manifest):
    for name, expected in manifest['frozen_hashes'].items():
        assert sha(ROOT / name) == expected, name
    for name, expected in manifest['protected_original_hashes'].items():
        assert sha(name) == expected, name
    for info in manifest['datasets'].values():
        assert sha(info['input']) == info['input_sha256'], info['input']
