"""
run_volcanic_vortex_ring_model.py
==================================
1D variable-entrainment vortex-ring / thermal rise model for NATURAL
VOLCANIC vortex rings rising through a real, stratified atmosphere.

This is the "volcanic" counterpart of the analog buoyant-vortex-ring
tank-experiment model (OneDThermalModel_PureFluid_height.m /
thermalfuncup_height_var_alpha.m). The governing entrainment/momentum/
buoyancy ODEs are unchanged; what changes is the physics of the ambient
medium and the source:

  - The working fluid is air (ideal gas): density depends on temperature
    AND on the local ambient pressure (which decreases with altitude), not
    on a fixed lab-fluid density.
  - Buoyancy is tracked via potential temperature so that the model
    automatically captures the adiabatic expansion/cooling the ring
    undergoes as it rises to lower ambient pressure (see atmosphere.py).
  - The ambient stratification N^2(z) is computed from a real atmospheric
    T(z), P(z) profile -- either the built-in ICAO/US Standard Atmosphere
    1976, or a user-supplied sounding (see load_atm_profile) -- evaluated
    at the TRUE altitude (vent altitude + height risen).
  - The vent (source) can be placed at any altitude above sea level.
  - No particle/ash loading is tracked (pure-gas ring), matching the
    "PureFluid" scope of the original model and the user's request that
    ring density be treated as that of air at the ring's temperature.
  - The entrainment coefficient can be the model's variable, height-
    dependent Ri_0-style closure (default) OR a single user-specified
    constant value (see entrainment_mode below).

This script is a thin driver around model_core.run_vortex_ring_model(),
which is shared with the GUI (vortex_ring_gui.py) so both stay consistent.

Source condition defaults below are for Etna's "Rings" vent (Bocca Nuova
crater, ~3380 m a.s.l.), based on thermal-camera measurements of vortex
ring events described in this project (see
INGV/Etna_Rings_Vent/EtnaRings_diam_alpha_vel_250902.m).

Johan Gilchrist
"""

import matplotlib
matplotlib.use("Agg")  # headless-safe; comment out if running interactively
import matplotlib.pyplot as plt

from model_core import run_vortex_ring_model, format_report
from atmosphere import load_atm_profile  # noqa: F401  (see example below)

# =============================================================================
# USER-DEFINED PARAMETERS
# =============================================================================

# --- Entrainment ---------------------------------------------------------
entrainment_mode = "variable"  # 'variable' (height-dependent, Ri_0-style
                                # closure) or 'constant' (single fixed alpha)
c_coef = 1.0        # O(1) coefficient in the variable-alpha closure
                     # (matches the original tank-model convention in
                     # thermalfuncup_height_var_alpha.m;
                     # EtnaRings_diam_alpha_vel_250902.m instead used
                     # c=0.26, Turner 1962 -- change here if desired).
                     # Ignored if entrainment_mode = 'constant'.
alpha_max = 1.0      # hard cap on the variable entrainment coefficient --
                     # see the detailed comment in vortex_ring_ode.make_rhs().
                     # Needed because volcanic source/ambient temperature
                     # contrasts are far larger than the tank experiments
                     # the alpha(Ri) closure was tuned on; without a cap the
                     # model can blow up to NaN within a few source radii of
                     # the vent. 1.0 is consistent with the largest near-vent
                     # entrainment coefficients fitted to real Rings-vent
                     # events (~0.95). Ignored if entrainment_mode = 'constant'.
alpha_const = 0.1    # fixed entrainment coefficient used only if
                     # entrainment_mode = 'constant'.

# --- Source (vent) conditions -- EDIT THESE FOR YOUR EVENT ------------------
z0_asl = 3380.0     # [m] vent altitude above sea level.
                     # Default: Etna's Bocca Nuova crater ("Rings" vent).
r0 = 5.0             # [m] vent/initial ring radius -- measured value for the
                     # event being modeled (observed range across 11 Rings
                     # events: ~0-19 m).
u0 = 20.0            # [m/s] initial (exit) vertical velocity -- measured
                     # value for the event being modeled
                     # (observed range: ~0.05-33 m/s).
Tj0_C = 450.0        # [degC] initial ring temperature at the vent -- measured
                     # value for Etna's "Rings" vent.

# --- Atmosphere --------------------------------------------------------------
atm_mode = "standard"  # 'standard' = built-in ICAO/US Standard Atmosphere 1976
                        # 'profile'  = user-supplied T(z), P(z) sounding
atm_profile = None
Ta0_target_C = 15.0     # [degC] measured/desired ambient temperature AT THE
                        # VENT ALTITUDE. The standard atmosphere alone is
                        # generally too cold near a volcanic summit (local
                        # heating, etc.), so model_core solves automatically
                        # for the uniform temperature shift that makes the
                        # standard atmosphere match this value at z0_asl.
                        # Ignored if atm_mode = 'profile' (supply a real
                        # sounding with the correct near-vent temperature
                        # instead in that case).

# Example: use a measured/radiosonde profile instead of the standard atmosphere
# atm_profile = load_atm_profile("etna_sounding.csv",
#                                 temp_units="C", press_units="hPa")
# atm_mode = "profile"

# --- Integration range --------------------------------------------------------
zmax_rise = 3000.0     # [m] max height ABOVE THE VENT to integrate to

# =============================================================================
# RUN
# =============================================================================
res = run_vortex_ring_model(
    z0_asl=z0_asl, r0=r0, u0=u0, Tj0_C=Tj0_C,
    atm_mode=atm_mode, atm_profile=atm_profile, Ta0_target_C=Ta0_target_C,
    entrainment_mode=entrainment_mode, c_coef=c_coef, alpha_max=alpha_max,
    alpha_const=alpha_const, zmax_rise=zmax_rise,
)

if not res["apex_reached"]:
    print(f"[info] integration reached zmax_rise = {zmax_rise:.0f} m above "
          f"the vent without the ring stopping (u never dropped to the "
          f"apex threshold).")
else:
    print(f"[info] ring reached its rise apex (u -> 0) at z = "
          f"{res['z_out'][-1]:.1f} m above the vent.")

print(format_report(res))

# =============================================================================
# PLOTS
# =============================================================================
z_out, V, M, B = res["z_out"], res["V"], res["M"], res["B"]
gprime, u, r, alpha_k_z = res["gprime"], res["u"], res["r"], res["alpha_k_z"]
Ri_z = res["Ri_z"]
Tp_z, Ta_z = res["Tp_z"], res["Ta_z"]

# Figure 1: fluxes (Volume, Momentum, Buoyancy)
fig, axs = plt.subplots(1, 3, figsize=(13, 5))
axs[0].plot(V, z_out, ".")
axs[0].set_title("Volume"); axs[0].set_xlabel("V [m^3]"); axs[0].set_ylabel("Height above vent [m]")
axs[1].plot(M, z_out, ".")
axs[1].set_title("Momentum"); axs[1].set_xlabel("M [m^4/s]")
axs[2].plot(B, z_out, ".")
axs[2].set_title("Buoyancy"); axs[2].set_xlabel("B [m^4/s^2]")
fig.tight_layout()
fig.savefig("volcanic_vortex_ring_fluxes.png", dpi=150)

# Figure 2: profiles (velocity, radius, reduced gravity, temperature)
fig2, axs2 = plt.subplots(1, 4, figsize=(16, 5))
axs2[0].plot(u, z_out, ".")
axs2[0].set_xlabel("Vertical velocity [m/s]"); axs2[0].set_ylabel("Height above vent [m]")
axs2[1].plot(r, z_out, ".")
axs2[1].set_xlabel("Radius [m]")
axs2[2].plot(gprime, z_out, ".")
axs2[2].set_xlabel("g' [m/s^2]")
axs2[3].plot(Tp_z - 273.15, z_out, ".", label="ring")
axs2[3].plot(Ta_z - 273.15, z_out, "--", label="ambient")
axs2[3].set_xlabel("Temperature [deg C]"); axs2[3].legend()
fig2.tight_layout()
fig2.savefig("volcanic_vortex_ring_profiles.png", dpi=150)

# Figure 3: entrainment and Richardson number
# Ri = |g'| r / u^2 diverges as u -> 0 near the rise apex, so (as in the
# original tank-experiment model) it's plotted on a log x-axis.
fig3, axs3 = plt.subplots(1, 2, figsize=(9, 5))
axs3[0].plot(alpha_k_z, z_out, ".")
axs3[0].set_xlabel(r"$\alpha_e$"); axs3[0].set_ylabel("Height above vent [m]"); axs3[0].set_xlim(0, 1)
axs3[1].semilogx(Ri_z, z_out, ".")
axs3[1].set_xlabel("Richardson number, Ri")
fig3.tight_layout()
fig3.savefig("volcanic_vortex_ring_entrainment_Ri.png", dpi=150)

print("\nSaved plots: volcanic_vortex_ring_fluxes.png, volcanic_vortex_ring_profiles.png, "
      "volcanic_vortex_ring_entrainment_Ri.png")
