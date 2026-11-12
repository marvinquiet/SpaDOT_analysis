"""Register raw H5AD inputs before fitting GraphST, SEDR, spaGCN, spaVAE,
spatiAlign, SpatialPCA or PRECAST independently at each stage.

Run with --raw-root DIR --output FILE and repeat --case NAME PATH DATATYPE.
PATH is relative to DIR, or absolute. DATATYPE is 10X or Stereo. Input files
must provide the timepoint, annotation, pixel_x and pixel_y observation fields
and matching spatial coordinates in obsm['spatial']. Provide raw counts in
layers['counts'] or X. Benchmark K is the number of annotated domains in each
stage. Annotation labels are otherwise used only for evaluation.
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np

from inputs import digest, save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--case", nargs=3, action="append", required=True,
                        metavar=("NAME", "PATH", "DATATYPE"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    records = {}
    for name, relative, datatype in args.case:
        if name in records or datatype not in {"10X", "Stereo"}:
            raise ValueError(f"Duplicate case or invalid data type: {name}, {datatype}")
        source = args.raw_root / relative
        data = ad.read_h5ad(source, backed="r")
        try:
            if not data.obs_names.is_unique or not data.var_names.is_unique:
                raise ValueError(f"Duplicate spot or gene identifiers in {source}")
            obs = data.obs
            if obs[["timepoint", "annotation"]].isna().any().any():
                raise ValueError(f"Missing stage or annotation in {source}")
            if not np.array_equal(data.obsm["spatial"], obs[["pixel_x", "pixel_y"]].to_numpy()):
                raise ValueError(f"Spatial coordinates disagree in {source}")
            counts = obs.groupby(obs.timepoint.astype(str), sort=True)["annotation"].nunique()
            records[name] = {"relative_input": relative, "sha256": digest(source),
                             "k": {key: int(value) for key, value in counts.items()},
                             "datatype": datatype}
        finally:
            data.file.close()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    save_json(args.output, records)


if __name__ == "__main__":
    main()
