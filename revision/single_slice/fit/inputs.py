"""Input and provenance helpers for the independent native method fits."""

import hashlib
import json
from pathlib import Path


def digest(path):
    with Path(path).open("rb") as stream:
        checksum = hashlib.sha256()
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def load_case(manifest, dataset):
    cases = json.loads(Path(manifest).read_text())
    if dataset not in cases:
        raise ValueError(f"Dataset or scenario {dataset!r} is absent from {manifest}")
    return cases[dataset]


def save_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)
