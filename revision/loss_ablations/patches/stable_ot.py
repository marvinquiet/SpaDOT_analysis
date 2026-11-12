"""Numerical OT fallback; same unbalanced objective and epsilon schedule.

Keep the released solver's exact result whenever it succeeds. On a detected
floating-point failure, solve the same problem with log-sum-exp updates. This
module changes no model, cost normalization, OT weights, or training schedule.
"""
import hashlib
import warnings
import numpy as np
from scipy.special import logsumexp, rel_entr, xlogy


def _objectives(C, p, q, log_a, log_b, epsilon, lambda1, lambda2):
    log_R = log_a[:, None] - C / epsilon + log_b[None, :]
    with np.errstate(under='ignore'):
        R = np.exp(log_R)
        K = np.exp(-C / epsilon)
    row = R.mean(axis=1)
    col = R.mean(axis=0)
    primal = (lambda1 * np.mean(rel_entr(row, p) - row + p)
              + lambda2 * np.mean(rel_entr(col, q) - col + q)
              + np.mean(epsilon * (xlogy(R, R) - R + K) + R * C))
    def conjugate(values, logs, weight):
        terms = np.zeros_like(values)
        positive = values > 0
        terms[positive] = -values[positive] * np.expm1(-epsilon * logs[positive] / weight)
        return weight * np.mean(terms)
    dual = (conjugate(p, log_a, lambda1) + conjugate(q, log_b, lambda2)
            - epsilon * np.mean(R - K))
    if not np.isfinite(primal) or not np.isfinite(dual) or not np.isfinite(R).all():
        raise FloatingPointError('Log-domain OT produced a non-finite objective or map')
    gap = abs(primal - dual) / max(abs(primal), np.finfo(np.float64).tiny)
    return R, float(primal), float(dual), float(gap)


def _relative_change(current, previous):
    shift = max(0., float(np.max(current)), float(np.max(previous)))
    x, y = np.exp(current - shift), np.exp(previous - shift)
    return np.linalg.norm(x - y) / (np.exp(-shift) + np.linalg.norm(x))


def solve_log(C, G, lambda1, lambda2, epsilon, batch_size, tolerance,
              tau, epsilon0, max_iter, profiling=False, **ignored):
    if ignored.get('use_C', False) or not ignored.get('use_Py', True):
        raise ValueError('Only the Python OT backend is available')
    C = np.asarray(C, dtype=np.float64)
    G = np.asarray(G, dtype=np.float64)
    if C.ndim != 2 or not all(C.shape) or G.shape != (C.shape[0],):
        raise ValueError('OT cost and growth shapes do not agree')
    if not np.isfinite(C).all() or not np.isfinite(G).all() or (G < 0).any():
        raise ValueError('OT inputs must be finite and growth must be nonnegative')
    if min(lambda1, lambda2, epsilon, epsilon0, tolerance, batch_size, max_iter) <= 0:
        raise ValueError('OT weights, epsilons, tolerance and iteration limits must be positive')
    if not G.any():
        return np.zeros_like(C), dict(backend='logsumexp', iterations=0,
            primal=0., dual=0., relative_duality_gap=0., converged=True)
    p, q = G.copy(), np.full(C.shape[1], np.mean(G))
    log_p = np.full(p.shape, -np.inf)
    log_p[p > 0] = np.log(p[p > 0])
    log_q = np.log(q)
    log_dx, log_dy = -np.log(C.shape[0]), -np.log(C.shape[1])
    log_a, log_b = np.zeros(C.shape[0]), np.zeros(C.shape[1])
    scale_factor = np.exp(-np.log(epsilon) / 5)
    epsilon_i = epsilon0 * scale_factor
    iterations = 0
    for scale in range(6):
        old_epsilon = epsilon_i
        epsilon_i /= scale_factor
        log_a = log_a * (old_epsilon / epsilon_i)
        log_b = log_b * (old_epsilon / epsilon_i)
        alpha1 = lambda1 / (lambda1 + epsilon_i)
        alpha2 = lambda2 / (lambda2 + epsilon_i)
        log_K = -C / epsilon_i
        threshold = tolerance if scale == 5 else 1e-6
        diagnostic = np.inf
        while diagnostic > threshold:
            for _ in range(int(batch_size) if scale == 5 else 5):
                old_a, old_b = log_a.copy(), log_b.copy()
                log_a = alpha1 * (log_p - logsumexp(log_K + log_b[None, :] + log_dy, axis=1))
                log_b = alpha2 * (log_q - logsumexp(log_K + log_a[:, None] + log_dx, axis=0))
                iterations += 1
                if iterations >= int(max_iter):
                    raise RuntimeError('Log-domain OT did not converge before max_iter was reached')
            if scale == 5:
                R, primal, dual, diagnostic = _objectives(
                    C, p, q, log_a, log_b, epsilon_i, lambda1, lambda2)
            else:
                diagnostic = max(_relative_change(log_a, old_a), _relative_change(log_b, old_b))
            if not np.isfinite(diagnostic):
                raise FloatingPointError('Non-finite log-domain OT convergence diagnostic')
        if profiling:
            print('Log-domain OT scale', scale, 'diagnostic', diagnostic, flush=True)
    assert (R >= 0).all() and np.isfinite(R).all()
    return R / C.shape[1], dict(backend='logsumexp', iterations=iterations,
        final_epsilon=float(epsilon_i), primal=primal, dual=dual,
        relative_duality_gap=float(diagnostic), converged=True)


def install(solver):
    """Install only in the new worker's imported module, never in frozen files."""
    native = solver.optimal_transport_duality_gap
    events = []
    def robust(C, G, lambda1, lambda2, epsilon, batch_size, tolerance,
               tau, epsilon0, max_iter, profiling=False, **ignored):
        args = dict(C=C, G=G, lambda1=lambda1, lambda2=lambda2, epsilon=epsilon,
            batch_size=batch_size, tolerance=tolerance, tau=tau, epsilon0=epsilon0,
            max_iter=max_iter, profiling=profiling, **ignored)
        event = dict(cost_sha256=hashlib.sha256(np.ascontiguousarray(C).tobytes()).hexdigest(),
                     growth_sha256=hashlib.sha256(np.ascontiguousarray(G).tobytes()).hexdigest())
        try:
            with warnings.catch_warnings(record=True) as captured, np.errstate(all='ignore'):
                result = native(**args)
            if not np.isfinite(result).all() or (result < 0).any():
                raise FloatingPointError('Released OT returned an invalid map')
            event.update(backend='released_numpy', exact_released_path=True,
                         captured_warnings=list(dict.fromkeys(str(w.message) for w in captured)))
        except FloatingPointError as error:
            result, info = solve_log(**args)
            event.update(info, exact_released_path=False, fallback_reason=str(error))
        events.append(event)
        return result
    solver.optimal_transport_duality_gap = robust
    return dict(native=native, events=events)
