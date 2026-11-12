from hashlib import sha256
from importlib.metadata import version
import ast
import json
from pathlib import Path
import argparse

import anndata
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits
from SpaDOT.utils import _analyze_utils

parser = argparse.ArgumentParser(description="Calculate ARI, PAS and CHAOS from each Mouse SpaDOT encoder branch.")
parser.add_argument("--latent", required=True, type=Path)
parser.add_argument("--output", required=True, type=Path)
parser.add_argument("--benchmark-scripts", default=Path("summarize_benchmark/scripts"), type=Path)
parser.add_argument("--verify-archive", action="store_true")
args = parser.parse_args()
out, benchmark, latent = args.output, args.benchmark_scripts, args.latent
out.mkdir(parents=True, exist_ok=False)
util = benchmark / "utils.py"
downstream = benchmark / "MouseEmbryo_Visium_Organogenesis/analyze_ATACvisium_noremoval_downstream.R"
readout = Path(_analyze_utils.__file__).resolve()

expected = {
    latent: "64c87116af1c32207692850641c7fd13de3b236ec9870cda134931063173e3df",
    util: "a5d31f792c73c7c141dacd3125f39ddaa9f8a527939c992424ee392a345ac717",
    readout: "6b4e43b0e5f88b3ddb8d064d6e510065a3648e2fe8f847fb54d384a5007578f2",
}
hashes = {str(path): sha256(path.read_bytes()).hexdigest() for path in (*expected, downstream)}
if args.verify_archive:
    for path, digest in expected.items():
        assert hashes[str(path)] == digest, (path, hashes[str(path)])
assert "compute_lisi" in downstream.read_text()
assert version("SpaDOT") == "1.0.7"

names = {"_compute_PAS", "_compute_CHAOS"}
functions = [node for node in ast.parse(util.read_text()).body if isinstance(node, ast.FunctionDef) and node.name in names]
assert {node.name for node in functions} == names
scope = {"np": np, "NearestNeighbors": NearestNeighbors, "StandardScaler": StandardScaler}
exec(compile(ast.Module(body=functions, type_ignores=[]), str(util), "exec"), scope)

a = anndata.read_h5ad(latent)
assert a.shape == (8255, 20) and a.obs_names.is_unique and np.isfinite(a.X).all()
np.testing.assert_array_equal(a.obs[["pixel_x", "pixel_y"]].to_numpy(), a.obsm["spatial"])
stages = {"t0": (2397, "E12.5"), "t1": (2312, "E13.5"), "t2": (3546, "E15.5")}
obs = a.obs[["timepoint", "pixel_x", "pixel_y"]].copy()
assert obs.timepoint.astype(str).value_counts().to_dict() == {t: n for t, (n, _) in stages.items()}
predictions, results = obs.copy(), []
with threadpool_limits(limits=1):
    for name, dims in (("GPVAE10", slice(0, 10)), ("GAT10", slice(10, 20)), ("combined20", slice(None))):
        read = _analyze_utils.KMeans_Clustering(a[:, dims].copy(), [5, 5, 5], seed=2026, n_init=10)
        assert read.obs_names.is_unique and set(read.obs_names) == set(obs.index)
        predictions[name] = read.obs.kmeans.reindex(obs.index)
        for tp, (count, stage) in stages.items():
            m = obs.timepoint == tp
            labels = predictions.loc[m, name].to_numpy().astype(str)
            truth = a.obs.loc[m, "annotation"].astype(str).str.replace("\n", " ").to_numpy()
            xy = obs.loc[m, ["pixel_x", "pixel_y"]].to_numpy(dtype=float)
            assert len(labels) == count and len(set(labels)) == 5
            results.append(dict(representation=name, timepoint=tp, stage=stage, n_spots=count,
                                PAS=scope["_compute_PAS"](xy, labels, k=10),
                                CHAOS=scope["_compute_CHAOS"](xy, labels),
                                ARI=adjusted_rand_score(truth, labels)))
predictions.to_csv(out / "labels.csv", index=False)
pd.DataFrame(results).to_csv(out / "pas_chaos_ari.csv", index=False)
(out / "provenance.json").write_text(json.dumps(dict(hashes=hashes, versions={k: version(k) for k in
    ("SpaDOT", "anndata", "scikit-learn", "numpy")}, dims={"GPVAE10":[0,10], "GAT10":[10,20]},
    K=[5,5,5], seed=2026, n_init=10, coordinate_order="saved latent spot order"), indent=2) + "\n")
