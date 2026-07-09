"""
vortex_ring_ode.py (web / Pyodide version)
============================================
Identical physics to the desktop package's vortex_ring_ode.py, just wired
to the dependency-free simple_ode.solve_ode() integrator instead of
scipy.integrate.solve_ivp (see simple_ode.py for why). The right-hand-side
function and entrainment closure are unchanged.

State vector y = [V, M, B, alpha]
    V     : ring volume                         [m^3]      ~ (4/3) pi r^3
    M     : ring volume x vertical velocity      [m^4/s]    ~ V*u
    B     : ring volume x reduced gravity        [m^4/s^2]  ~ V*g'
    alpha : entrainment coefficient (diagnostic state)

    dV/dz = 3 * alpha * V^(2/3)
    dM/dz = (2/3) * B * V / M
    dB/dz = -N^2(z) * V

Johan Gilchrist
"""

import numpy as np
from atmosphere import brunt_vaisala_N2


def make_rhs(z0_asl, atm_mode, atm_profile, T_offset, c_coef, alpha_max=1.0,
             entrainment_mode="variable", alpha_const=0.1):
    """Return a `rhs(z_above_vent, y)` closure. See the desktop package's
    vortex_ring_ode.py for full documentation of entrainment_mode/alpha_max."""

    def rhs(z_above_vent, y):
        V, M, B, _alpha_state = y
        r = (3.0 * V / (4.0 * np.pi)) ** (1.0 / 3.0)
        u = M / V

        if entrainment_mode == "constant":
            alpha = alpha_const
        else:
            alpha = abs(B) / (2.0 * np.pi * c_coef * (r * u) ** 2)
            alpha = min(alpha, alpha_max)

        N2 = brunt_vaisala_N2(z0_asl + z_above_vent, atm_mode, atm_profile, T_offset)

        dVdz = 3.0 * alpha * V ** (2.0 / 3.0)
        dMdz = (2.0 / 3.0) * B * V / M
        dBdz = -N2 * V
        dAlphadz = alpha - _alpha_state

        return np.array([dVdz, dMdz, dBdz, dAlphadz])

    return rhs


def make_apex_event(u_min=1e-3):
    """Return an `event(z_above_vent, y) -> float` for simple_ode.solve_ode,
    which stops integration the first time it crosses from positive to
    non-positive (same convention as the desktop package)."""

    def event(z_above_vent, y):
        V, M, _B, _alpha_state = y
        return M / V - u_min

    return event
