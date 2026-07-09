"""
model_core.py
==============
Shared model-driver logic for the 1D volcanic vortex ring rise model.

Both the command-line script (run_volcanic_vortex_ring_model.py) and the
GUI (vortex_ring_gui.py) call `run_vortex_ring_model()` here, so the two
front-ends can never drift out of sync with each other.

Johan Gilchrist
"""

import numpy as np
from scipy.integrate import solve_ivp

from atmosphere import (
    atmosphere_TP,
    potential_temperature,
    brunt_vaisala_N2,
    G0,
    KAPPA,
    P0REF,
)
from vortex_ring_ode import make_rhs, make_apex_event


class ModelError(RuntimeError):
    """Raised for invalid model inputs or a failed integration."""


def run_vortex_ring_model(
    z0_asl,
    r0,
    u0,
    Tj0_C,
    atm_mode="standard",
    atm_profile=None,
    Ta0_target_C=15.0,
    entrainment_mode="variable",
    c_coef=1.0,
    alpha_max=1.0,
    alpha_const=0.1,
    zmax_rise=3000.0,
    n_pts=10000,
    alphaup=0.18,
    u_min=1e-3,
):
    """Integrate the 1D volcanic vortex ring rise model and return a dict
    of results and diagnostics.

    Parameters
    ----------
    z0_asl : float
        Vent altitude above sea level [m].
    r0, u0, Tj0_C : float
        Initial ring radius [m], exit velocity [m/s], exit temperature [degC].
    atm_mode : {'standard','profile'}
        'standard' uses the built-in ICAO/US Standard Atmosphere 1976,
        shifted so that its temperature AT THE VENT equals Ta0_target_C
        (see T_offset below). 'profile' uses a user-supplied sounding
        (atm_profile, from atmosphere.load_atm_profile) as-is.
    atm_profile : dict or None
        Required if atm_mode='profile' (see atmosphere.load_atm_profile).
    Ta0_target_C : float
        Desired/measured ambient temperature at the vent altitude [degC].
        Only used when atm_mode='standard'.
    entrainment_mode : {'variable','constant'}
        'variable' uses the Ri_0-style closure (alpha varies with height,
        the model's default). 'constant' uses a single fixed alpha_const
        for the whole rise.
    c_coef : float
        O(1) coefficient in the variable-alpha closure (ignored if
        entrainment_mode='constant').
    alpha_max : float
        Cap on the variable entrainment coefficient (ignored if
        entrainment_mode='constant'). See vortex_ring_ode.make_rhs.
    alpha_const : float
        Fixed entrainment coefficient (used only if entrainment_mode='constant').
    zmax_rise : float
        Maximum height above the vent to integrate to [m].
    n_pts : int
        Number of points in the returned height grid.
    alphaup : float
        Top-hat entrainment coefficient used only in the plume-scaling
        comparison law (not in the ODE itself).
    u_min : float
        Vertical velocity [m/s] at which integration is stopped (rise apex).

    Returns
    -------
    dict with height profiles (z_out, V, M, B, alpha_k_z, u, r, gprime,
    Ta_z, Pa_z, Tp_z) and scalar diagnostics (Ri0, alpha0, N2_vent,
    scaling-law estimates, z_LNB, alpha_avg, alpha_avg_LNB, apex_reached,
    T_offset, and the echoed input parameters).
    """
    if r0 <= 0 or u0 <= 0:
        raise ModelError("Vent radius and exit velocity must both be > 0.")
    if entrainment_mode not in ("variable", "constant"):
        raise ModelError("entrainment_mode must be 'variable' or 'constant'.")
    if entrainment_mode == "constant" and alpha_const <= 0:
        raise ModelError("alpha_const must be > 0 for constant-entrainment mode.")

    Tj0 = Tj0_C + 273.15

    # Resolve T_offset so the standard atmosphere gives Ta0_target_C at the
    # vent. Ignored (0) for a user-supplied profile -- supply a real
    # sounding with the correct near-vent temperature in that case instead.
    if atm_mode == "profile":
        if atm_profile is None or "z" not in atm_profile:
            raise ModelError(
                "atm_mode='profile' requires atm_profile data "
                "(see atmosphere.load_atm_profile)."
            )
        T_offset = 0.0
    else:
        Ta0_std, _ = atmosphere_TP(z0_asl, "standard", None, 0.0)
        T_offset = (Ta0_target_C + 273.15) - Ta0_std

    Ta0, Pa0 = atmosphere_TP(z0_asl, atm_mode, atm_profile, T_offset)
    theta_a0 = potential_temperature(Ta0, Pa0)
    theta_j0 = potential_temperature(Tj0, Pa0)  # ring assumed at local ambient P

    V0 = (4.0 / 3.0) * np.pi * r0 ** 3
    M0 = V0 * u0
    B0 = V0 * G0 * (theta_j0 - theta_a0) / theta_a0
    K0 = u0 * r0  # [m^2/s] initial circulation scale
    gprime0 = B0 / V0  # [m/s^2] initial reduced gravity (potential-temperature based)

    # Source (initial) Richardson number, using the SAME reduced-gravity
    # definition as the rest of the model (potential temperature -- see
    # atmosphere.brunt_vaisala_N2) and as the Ri_z height profile below, so
    # this is exactly Ri_z's starting value and is consistent with alpha0.
    Ri0 = abs(gprime0) * r0 / u0 ** 2

    if entrainment_mode == "constant":
        alpha0 = alpha_const
    else:
        alpha0 = min(abs(B0) / (2.0 * np.pi * c_coef * K0 ** 2), alpha_max)

    y0 = [V0, M0, B0, alpha0]

    zmin_rise = 2.0 * r0
    z_eval = np.linspace(zmin_rise, zmax_rise, n_pts)

    rhs = make_rhs(
        z0_asl, atm_mode, atm_profile, T_offset, c_coef,
        alpha_max=alpha_max, entrainment_mode=entrainment_mode,
        alpha_const=alpha_const,
    )
    apex_event = make_apex_event(u_min=u_min)

    sol = solve_ivp(
        rhs,
        (z_eval[0], z_eval[-1]),
        y0,
        method="LSODA",
        rtol=1e-8,
        atol=1e-10,
        max_step=5.0,
        events=apex_event,
        dense_output=True,
    )

    if sol.t.size < 2:
        raise ModelError(f"Integration failed immediately: {sol.message}")

    apex_reached = sol.status == 1

    z_out = np.linspace(z_eval[0], sol.t[-1], n_pts)
    V, M, B, alpha_k_z = sol.sol(z_out)
    u = M / V
    r = (3.0 * V / (4.0 * np.pi)) ** (1.0 / 3.0)
    gprime = B / V
    Ri_z = np.abs(gprime) * r / u ** 2  # local Richardson number height profile

    Ta_z, Pa_z = atmosphere_TP(z0_asl + z_out, atm_mode, atm_profile, T_offset)
    theta_a_z = potential_temperature(Ta_z, Pa_z)
    theta_p_z = theta_a_z + gprime * theta_a_z / G0
    Tp_z = theta_p_z * (Pa_z / P0REF) ** KAPPA

    # Classical scaling-law estimates (approximate: assume constant N^2,
    # evaluated at the vent; only leading-order in a real, height-varying
    # atmosphere).
    Cth, Cpl = 1.0, 2.3
    N2_vent = brunt_vaisala_N2(z0_asl, atm_mode, atm_profile, T_offset)
    alpha_for_scaling = alpha_const if entrainment_mode == "constant" else alpha0
    zmax_thermal_scaling = (
        Cth * alpha_for_scaling ** (-3.0 / 4.0) * B0 ** 0.25 * N2_vent ** (-0.25)
        - r0 / alpha_for_scaling
    )
    zmax_vortring_scaling = Cth * K0 ** 1.5 * B0 ** (-0.5) * N2_vent ** (-0.25)
    zmax_plume_scaling = Cpl * alphaup ** (-0.5) * B0 ** 0.25 * N2_vent ** (-0.375)

    ilnb = int(np.argmin(np.abs(gprime)))
    z_LNB = z_out[ilnb]
    alpha_avg = float(np.mean(alpha_k_z))
    alpha_avg_LNB = float(np.mean(alpha_k_z[: ilnb + 1]))

    return dict(
        # height profiles
        z_out=z_out, V=V, M=M, B=B, alpha_k_z=alpha_k_z, u=u, r=r, gprime=gprime,
        Ri_z=Ri_z, Ta_z=Ta_z, Pa_z=Pa_z, Tp_z=Tp_z,
        # echoed inputs
        z0_asl=z0_asl, r0=r0, u0=u0, Tj0_C=Tj0_C,
        atm_mode=atm_mode, entrainment_mode=entrainment_mode,
        c_coef=c_coef, alpha_max=alpha_max, alpha_const=alpha_const,
        T_offset=T_offset,
        # scalar diagnostics
        Ta0=Ta0, Pa0=Pa0, Ri0=Ri0, alpha0=alpha0, N2_vent=N2_vent,
        zmax_thermal_scaling=zmax_thermal_scaling,
        zmax_vortring_scaling=zmax_vortring_scaling,
        zmax_plume_scaling=zmax_plume_scaling,
        z_LNB=z_LNB, alpha_avg=alpha_avg, alpha_avg_LNB=alpha_avg_LNB,
        apex_reached=apex_reached,
    )


def format_report(res):
    """Format the diagnostic summary of a run_vortex_ring_model() result
    dict as a multi-line string (used by both the CLI script and the GUI)."""
    lines = []
    lines.append(f"Vent altitude:                        {res['z0_asl']:8.1f} m a.s.l.")
    lines.append(f"Ambient T, P at vent:                  {res['Ta0']-273.15:6.1f} degC, {res['Pa0']/100:7.1f} hPa")
    lines.append(f"Initial ring T, radius, velocity:      {res['Tj0_C']:6.1f} degC, {res['r0']:5.2f} m, {res['u0']:5.2f} m/s")
    lines.append(f"Entrainment mode:                      {res['entrainment_mode']}")
    lines.append(f"Initial entrainment coefficient:        {res['alpha0']: .4f}")
    lines.append(f"Source Richardson number Ri0:           {res['Ri0']: .4f}")
    lines.append(f"N^2 at vent:                            {res['N2_vent']: .3e} s^-2")
    lines.append("--- classical scaling estimates (approx, const-N^2) ---")
    lines.append(f"  thermal scaling:   {res['zmax_thermal_scaling']:8.1f} m")
    lines.append(f"  vortex-ring scal.: {res['zmax_vortring_scaling']:8.1f} m")
    lines.append(f"  plume scaling:     {res['zmax_plume_scaling']:8.1f} m")
    lines.append("--- integrated model result ---")
    apex_note = "" if res["apex_reached"] else "  [integration stopped before reaching the apex]"
    lines.append(f"Max height risen above vent:           {res['z_out'][-1]:8.1f} m{apex_note}")
    lines.append(f"Level of neutral buoyancy (rise):       {res['z_LNB']:8.1f} m")
    lines.append(f"Level of neutral buoyancy (a.s.l.):     {res['z0_asl']+res['z_LNB']:8.1f} m")
    lines.append(f"Height-averaged alpha_e (0-LNB):        {res['alpha_avg_LNB']:8.4f}")
    lines.append(f"Height-averaged alpha_e (full rise):    {res['alpha_avg']:8.4f}")
    return "\n".join(lines)
