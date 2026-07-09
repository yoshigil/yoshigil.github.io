"""
vortex_ring_ode.py
===================
Height-domain ODE right-hand side for the 1D variable-entrainment
buoyant vortex ring / thermal rise model, generalized to a real
(compressible) atmosphere.

State vector y = [V, M, B, alpha]
    V     : ring volume                         [m^3]      ~ (4/3) pi r^3
    M     : ring volume x vertical velocity      [m^4/s]    ~ V*u
    B     : ring volume x reduced gravity        [m^4/s^2]  ~ V*g'
    alpha : entrainment coefficient (diagnostic state, updated algebraically
            each step so it can be output alongside V,M,B by the solver)

Governing equations (Morton, Taylor & Turner 1956; Turner 1962; Sparks et
al. 1997), integrated with height z (above the vent) as the independent
variable instead of time:

    dV/dz = 3 * alpha * V^(2/3)
    dM/dz = (2/3) * B * V / M
    dB/dz = -N^2(z) * V

The only physics that changes relative to the analog tank-experiment model
is:
  (1) B is defined via potential-temperature excess rather than a fixed-fluid
      density excess, so that the ring's own adiabatic expansion/cooling as
      it rises to lower ambient pressure is captured automatically; and
  (2) N^2(z) is the REAL, height-varying stratification of the atmosphere
      (standard atmosphere or a supplied sounding), evaluated at the true
      altitude z0_asl + z, rather than a single constant value.

The variable-entrainment closure (alpha as a function of the local buoyancy
and momentum flux) is unchanged from thermalfuncup_height_var_alpha.m.
"""

import numpy as np
from atmosphere import brunt_vaisala_N2


def make_rhs(z0_asl, atm_mode, atm_profile, T_offset, c_coef, alpha_max=1.0,
             entrainment_mode="variable", alpha_const=0.1):
    """Return a `rhs(z_above_vent, y)` closure bound to the given source
    altitude, atmosphere settings, and entrainment closure.

    entrainment_mode : {'variable','constant'}
        'variable' (default) uses the Ri_0-style closure below (same
        functional form as thermalfuncup_height_var_alpha.m).
        'constant' instead uses a single fixed alpha_const for the whole
        rise (analogous to thermalfuncup_height_const_alpha.m /
        the classic constant-entrainment plume/thermal model).
    alpha_const : float
        Fixed entrainment coefficient used when entrainment_mode='constant'.
        Ignored otherwise.
    alpha_max : float
        Hard cap on the variable entrainment coefficient (ignored in
        'constant' mode). A volcanic ring exiting into cold, high-altitude
        ambient air can have a MUCH larger relative buoyancy anomaly than
        the small-density-contrast analog tank experiments the alpha(Ri)
        closure was built for (fractional potential-temperature excess of
        order 1 rather than a few percent). Left uncapped, the closure's
        positive feedback (higher alpha -> more dilution of momentum ->
        even higher alpha) can produce a genuine finite-height blow-up of
        V. Capping alpha is both a standard numerical safeguard for
        variable-entrainment closures of this type (see e.g. Kaminski et
        al. 2005) and physically grounded here: fitted near-vent ("jet
        development") entrainment coefficients for Etna's Rings vent
        events reach up to ~0.95 (see
        INGV/Etna_Rings_Vent/Rings_alphae.m), so alpha_max = 1.0 is a
        reasonable, data-informed ceiling. Set to np.inf to disable.
    """

    def rhs(z_above_vent, y):
        V, M, B, _alpha_state = y
        r = (3.0 * V / (4.0 * np.pi)) ** (1.0 / 3.0)
        u = M / V

        if entrainment_mode == "constant":
            alpha = alpha_const
        else:
            # Variable entrainment coefficient, Ri_0-style closure
            # (Turner 1979, Ch. 6; same functional form as
            # thermalfuncup_height_var_alpha.m), capped at alpha_max (see above).
            alpha = abs(B) / (2.0 * np.pi * c_coef * (r * u) ** 2)
            alpha = min(alpha, alpha_max)

        N2 = brunt_vaisala_N2(z0_asl + z_above_vent, atm_mode, atm_profile, T_offset)

        dVdz = 3.0 * alpha * V ** (2.0 / 3.0)
        dMdz = (2.0 / 3.0) * B * V / M
        dBdz = -N2 * V
        dAlphadz = alpha - _alpha_state  # forces the 4th state to track alpha

        return [dVdz, dMdz, dBdz, dAlphadz]

    return rhs


def make_apex_event(u_min=1e-3):
    """Return a `solve_ivp` event that fires when the ring's vertical
    velocity u = M/V decays to `u_min` [m/s]. The height-domain momentum
    equation dM/dz = (2/3) B V / M is singular as u -> 0 (the ring's rise
    apex), so integration should be stopped at/near that point rather than
    letting the solver grind through ever-shrinking step sizes trying to
    push past it.
    """

    def event(z_above_vent, y):
        V, M, _B, _alpha_state = y
        return M / V - u_min

    event.terminal = True
    event.direction = -1
    return event
