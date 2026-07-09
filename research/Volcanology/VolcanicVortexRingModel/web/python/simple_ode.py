"""
simple_ode.py
==============
A small, dependency-free (numpy only) adaptive Dormand-Prince RK45 ODE
integrator with event detection, written to replace scipy.integrate.solve_ivp
for the in-browser (Pyodide) version of the volcanic vortex ring model.

Why not just use scipy in Pyodide? Pyodide *can* load scipy, but it's a
large (tens of MB) WebAssembly package and noticeably slows the first page
load. This module implements the same explicit Runge-Kutta-Dormand-Prince
5(4) scheme that scipy's own `RK45` method uses (Dormand & Prince, 1980),
with adaptive step-size control and simple linear-interpolation resampling
in place of scipy's true polynomial dense output -- accurate enough for
this model's height profiles, at a fraction of the dependency weight.

Validated against scipy.integrate.solve_ivp (LSODA) for the vortex ring
model's reference case; see test_simple_ode.py.

Johan Gilchrist
"""

import numpy as np

# Dormand-Prince RK45 (DOPRI5) Butcher tableau
_C2, _C3, _C4, _C5, _C6 = 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0

_A21 = 1 / 5
_A31, _A32 = 3 / 40, 9 / 40
_A41, _A42, _A43 = 44 / 45, -56 / 15, 32 / 9
_A51, _A52, _A53, _A54 = 19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729
_A61, _A62, _A63, _A64, _A65 = 9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656

# 5th-order solution weights (b2 = 0)
_B1, _B3, _B4, _B5, _B6 = 35 / 384, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84
# 4th-order (embedded) solution weights, for error estimation
_BS1, _BS3, _BS4, _BS5, _BS6, _BS7 = (
    5179 / 57600, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40,
)


class OdeResult:
    """Minimal result container, loosely mirroring scipy's OdeResult."""

    def __init__(self, t, y, status, message):
        self.t = np.asarray(t)               # shape (n_accepted_steps,)
        self.y = np.asarray(y).T              # shape (n_states, n_accepted_steps)
        self.status = status                  # 0 = reached t_bound, 1 = event, -1 = failed
        self.message = message

    @property
    def success(self):
        return self.status >= 0


def _error_norm(err, y_old, y_new, rtol, atol):
    scale = atol + rtol * np.maximum(np.abs(y_old), np.abs(y_new))
    return float(np.sqrt(np.mean((err / scale) ** 2)))


def _rk45_step(fun, t, y, h, k1):
    k2 = fun(t + _C2 * h, y + h * (_A21 * k1))
    k3 = fun(t + _C3 * h, y + h * (_A31 * k1 + _A32 * k2))
    k4 = fun(t + _C4 * h, y + h * (_A41 * k1 + _A42 * k2 + _A43 * k3))
    k5 = fun(t + _C5 * h, y + h * (_A51 * k1 + _A52 * k2 + _A53 * k3 + _A54 * k4))
    k6 = fun(t + _C6 * h, y + h * (_A61 * k1 + _A62 * k2 + _A63 * k3 + _A64 * k4 + _A65 * k5))
    y5 = y + h * (_B1 * k1 + _B3 * k3 + _B4 * k4 + _B5 * k5 + _B6 * k6)
    k7 = fun(t + h, y5)  # FSAL: reused as k1 of the next step if accepted
    y4 = y + h * (
        _BS1 * k1 + _BS3 * k3 + _BS4 * k4 + _BS5 * k5 + _BS6 * k6 + _BS7 * k7
    )
    err = y5 - y4
    return y5, err, k7


def solve_ode(fun, t_span, y0, max_step=np.inf, rtol=1e-6, atol=1e-9,
              event=None, first_step=None, max_steps=200000):
    """Integrate dy/dt = fun(t, y) from t_span[0] to t_span[1].

    Parameters
    ----------
    fun : callable(t, y) -> array_like
        Right-hand side. y is a 1D numpy array.
    t_span : (t0, t_bound)
    y0 : array_like
        Initial state.
    max_step : float
        Hard cap on the step size.
    rtol, atol : float
        Relative/absolute error tolerances (same convention as scipy).
    event : callable(t, y) -> float, optional
        If given, integration stops the first time this function's value
        crosses from positive to negative (same convention as a scipy
        solve_ivp event with direction=-1, terminal=True). The crossing is
        located to within `rtol` of the step size by bisection.
    first_step : float, optional
        Initial step size guess. Defaults to a small fraction of the span.
    max_steps : int
        Safety limit on the number of accepted steps.

    Returns
    -------
    OdeResult with .t (1D) and .y (n_states x n_points), evaluated only at
    the (adaptively chosen) accepted steps -- resample with
    numpy.interp / resample_solution() for a uniform output grid.
    """
    t0, t_bound = t_span
    y = np.array(y0, dtype=float)
    t = float(t0)
    ts = [t]
    ys = [y.copy()]

    k1 = np.asarray(fun(t, y), dtype=float)

    if first_step is None:
        h = min(max_step, (t_bound - t0) / 100.0) if np.isfinite(max_step) else (t_bound - t0) / 100.0
        h = max(h, 1e-8)
    else:
        h = first_step

    status = 0
    message = "Reached t_bound."
    n_steps = 0

    while t < t_bound and n_steps < max_steps:
        h = min(h, max_step, t_bound - t)
        if h <= 0:
            break

        y_new, err, k7 = _rk45_step(fun, t, y, h, k1)
        err_norm = _error_norm(err, y, y_new, rtol, atol)

        if err_norm <= 1.0 or h <= 1e-12:
            t_new = t + h

            if event is not None:
                v_old = event(t, y)
                v_new = event(t_new, y_new)
                if v_old > 0 and v_new <= 0:
                    t_event, y_event = _refine_event(fun, event, t, y, t_new, y_new)
                    ts.append(t_event)
                    ys.append(y_event)
                    status = 1
                    message = "A termination event occurred."
                    return OdeResult(ts, ys, status, message)

            t, y, k1 = t_new, y_new, k7
            ts.append(t)
            ys.append(y.copy())
            n_steps += 1

            # Step-size growth factor (standard PI-ish controller)
            factor = 5.0 if err_norm == 0 else min(5.0, max(0.2, 0.9 * err_norm ** (-0.2)))
            h = h * factor
        else:
            h = h * max(0.1, 0.9 * err_norm ** (-0.25))

    if n_steps >= max_steps:
        status = -1
        message = "Maximum number of steps exceeded."

    return OdeResult(ts, ys, status, message)


def _refine_event(fun, event, t_lo, y_lo, t_hi, y_hi, tol=1e-6, max_iter=50):
    """Bisect for the root of `event` between (t_lo,y_lo) [positive] and
    (t_hi,y_hi) [non-positive], re-integrating with a single RK45 step of
    the appropriate size each time (adequate since the bracket is already
    small -- one accepted adaptive step)."""
    v_lo = event(t_lo, y_lo)
    for _ in range(max_iter):
        t_mid = 0.5 * (t_lo + t_hi)
        h = t_mid - t_lo
        k1 = np.asarray(fun(t_lo, y_lo), dtype=float)
        y_mid, _err, _k7 = _rk45_step(fun, t_lo, y_lo, h, k1)
        v_mid = event(t_mid, y_mid)
        if abs(t_hi - t_lo) < tol:
            return t_mid, y_mid
        if (v_mid > 0) == (v_lo > 0):
            t_lo, y_lo, v_lo = t_mid, y_mid, v_mid
        else:
            t_hi, y_hi = t_mid, y_mid
    return t_hi, y_hi


def resample_solution(sol, t_uniform):
    """Linearly interpolate an OdeResult onto a uniform grid t_uniform.
    Returns an array of shape (n_states, len(t_uniform))."""
    n_states = sol.y.shape[0]
    out = np.empty((n_states, len(t_uniform)))
    for i in range(n_states):
        out[i, :] = np.interp(t_uniform, sol.t, sol.y[i, :])
    return out
