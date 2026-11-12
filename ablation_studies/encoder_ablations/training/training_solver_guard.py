"""Observe the release solver without changing costs, updates, RNG or returns."""
from contextlib import contextmanager
from pathlib import Path
import hashlib
import json
import sys
import time

SOLVER_SHA = "34d2c557045e7c2f0166e647684fa6b711cab26605ff80a6750ed43de9bec805"
EXPECTED = dict(epsilon=.05, epsilon0=1., lambda1=.1, lambda2=5., tau=1000.)


def digest(value):
    import numpy as np
    x = np.ascontiguousarray(value)
    h = hashlib.sha256()
    h.update(str(x.dtype).encode())
    h.update(str(x.shape).encode())
    h.update(x.tobytes())
    return h.hexdigest()


@contextmanager
def observe_training_solver(root, destination):
    """Wrap only a call to the already installed fixed solver, with no retries."""
    import numpy as np
    from SpaDOT.utils import _train_utils as U
    from SpaDOT.utils.OT_loss import ot_solvers as S
    root, destination = Path(root).resolve(), Path(destination).resolve()
    source = root / "source/SpaDOT/utils/OT_loss/ot_solvers.py"
    assert Path(S.__file__).resolve() == source
    assert hashlib.sha256(source.read_bytes()).hexdigest() == SOLVER_SHA
    assert U.compute_transport_map is S.compute_transport_map
    assert U.compute_transport_map.__globals__["optimal_transport_duality_gap"] is S.optimal_transport_duality_gap
    original = S.optimal_transport_duality_gap
    destination.mkdir(parents=True, exist_ok=False)
    calls = []
    header = dict(solver_sha256=SOLVER_SHA, expected_parameters=EXPECTED,
                  training_growth_iters=3, training_returned_growth_index=0,
                  no_algorithm_or_parameter_change=True, automatic_retry=False,
                  no_extra_neural_forward=True, no_RNG_calls=True, calls=calls)

    def write():
        target = destination / "execution.json"
        tmp = destination / "execution.json.tmp"
        tmp.write_text(json.dumps(header, indent=2, allow_nan=False) + "\n")
        tmp.replace(target)

    def observed(**params):
        for k, value in EXPECTED.items():
            assert params[k] == value, (k, params[k], value)
        assert sys.getprofile() is None, "Unexpected concurrent profiler"
        captured = {}
        def profile(frame, event, arg):
            if event == "return" and frame.f_code is original.__code__:
                for k in ["epsilon_i", "duality_gap", "pri", "dua"]:
                    if k in frame.f_locals:
                        captured[k] = float(frame.f_locals[k])
        record = dict(index=len(calls), started_unix=time.time(),
                      parameters={k: params[k] for k in EXPECTED},
                      tolerance=float(params["tolerance"]),
                      C_shape=list(np.shape(params["C"])),
                      C_sha256=digest(params["C"]), G_sha256=digest(params["G"]))
        calls.append(record)
        try:
            sys.setprofile(profile)
            try:
                result = original(**params)
            finally:
                sys.setprofile(None)
            assert set(captured) == {"epsilon_i", "duality_gap", "pri", "dua"}
            assert np.isfinite(list(captured.values())).all()
            assert captured["duality_gap"] <= params["tolerance"]
            assert np.isclose(captured["epsilon_i"], params["epsilon"], atol=1e-15, rtol=1e-14)
            assert np.isfinite(result).all() and (result >= 0).all()
            file = destination / ("call_%04d.npz" % record["index"])
            np.savez_compressed(file, C=params["C"], G=params["G"], returned=result)
            record.update(status="returned", convergence=captured,
                          returned_sha256=digest(result), vectors_file=file.name)
            return result
        except Exception as error:
            record.update(status="failed", exception_type=type(error).__name__, error=str(error))
            raise
        finally:
            record["finished_unix"] = time.time()
            write()

    S.optimal_transport_duality_gap = observed
    write()
    try:
        assert U.compute_transport_map.__globals__["optimal_transport_duality_gap"] is observed
        yield header
    finally:
        S.optimal_transport_duality_gap = original
        header["original_function_restored"] = True
        header["all_recorded_calls_returned"] = all(c["status"] == "returned" for c in calls)
        write()


def numerical_failure(status, root, family):
    """Narrow prospective classification; integrity/OOM/baseline errors stop."""
    if family == "dual" or status.get("stage") != "training":
        return False
    text = status.get("traceback", "")
    if str(Path(root) / "source/SpaDOT/utils/_train_utils.py") not in text:
        return False
    error = status.get("error", "")
    fixed_errors = {
        "Non-finite OT primal or dual objective", "OT did not converge before max_iter was reached",
        "Non-finite OT convergence diagnostic", "OT returned an invalid transport map"}
    if error in fixed_errors:
        return ("_update_OT_matrix" in text and "optimal_transport_duality_gap" in text
                and str(Path(root) / "source/SpaDOT/utils/OT_loss/ot_solvers.py") in text
                and status.get("exception_type") in {"FloatingPointError", "RuntimeError"})
    if error == "Non-finite recorded training loss":
        return status.get("exception_type") == "FloatingPointError" and "observed_update" in text
    return (status.get("exception_type") == "ValueError" and "_update_Kmeans" in text
            and "sklearn" in text
            and (error.startswith("Input X contains NaN.")
                 or error.startswith("Input X contains infinity or a value too large")))
