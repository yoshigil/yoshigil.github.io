"""
atmosphere.py
=============
Atmospheric thermodynamics helpers for the volcanic vortex ring rise model.

Provides:
    - standard_atmosphere(z)   : ICAO / US Standard Atmosphere 1976, T(z), P(z)
    - load_atm_profile(path)   : read a user-supplied height/T/P sounding
    - atmosphere_TP(z, ...)    : dispatcher -> standard atmosphere or user profile
    - potential_temperature(T,P)
    - air_density(T, z, ...)   : ideal-gas density of a parcel at temperature T
                                  sitting at the local ambient pressure P(z)
                                  (i.e. assumes fast dynamic pressure equilibration
                                  between the ring and its surroundings -- the
                                  standard assumption used in 1D plume/thermal
                                  theory, e.g. Morton 1957; Woods 1988)
    - brunt_vaisala_N2(z, ...) : real, height-dependent stratification frequency
                                  N^2(z) = (g/theta_a) d(theta_a)/dz, computed
                                  from the supplied/standard atmosphere profile

Johan Gilchrist
"""

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants (dry air)
# ---------------------------------------------------------------------------
G0 = 9.80665          # [m/s^2] standard gravitational acceleration
RD = 287.053           # [J/(kg K)] specific gas constant, dry air
CP = 1004.7             # [J/(kg K)] specific heat of dry air, constant pressure
KAPPA = RD / CP          # ~0.2854, Poisson constant
P0REF = 100000.0         # [Pa] reference pressure for potential temperature (1000 hPa)

# ---------------------------------------------------------------------------
# ICAO / US Standard Atmosphere 1976 layer definitions (valid 0-84852 m)
# ---------------------------------------------------------------------------
_ZB = np.array([0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 84852.0])
_LB = np.array([-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002, 0.0])
_TB = np.zeros_like(_ZB)
_PB = np.zeros_like(_ZB)
_TB[0] = 288.15
_PB[0] = 101325.0
for _i in range(1, len(_ZB)):
    _dz = _ZB[_i] - _ZB[_i - 1]
    _Tb, _Lb, _Pb = _TB[_i - 1], _LB[_i - 1], _PB[_i - 1]
    _Tt = _Tb + _Lb * _dz
    if abs(_Lb) > 1e-12:
        _Pt = _Pb * (_Tb / _Tt) ** (G0 / (RD * _Lb))
    else:
        _Pt = _Pb * np.exp(-G0 * _dz / (RD * _Tb))
    _TB[_i] = _Tt
    _PB[_i] = _Pt


def standard_atmosphere(z):
    """ICAO / US Standard Atmosphere 1976 temperature and pressure.

    Parameters
    ----------
    z : scalar or array_like
        Geometric altitude above sea level [m]. Valid for 0-84852 m.

    Returns
    -------
    T : temperature [K] (scalar if z is scalar, else ndarray)
    P : pressure [Pa]    (scalar if z is scalar, else ndarray)
    """
    z_arr = np.atleast_1d(np.asarray(z, dtype=float))
    T = np.empty_like(z_arr)
    P = np.empty_like(z_arr)
    idx = np.clip(np.searchsorted(_ZB, z_arr, side="right") - 1, 0, len(_ZB) - 1)
    for k in range(z_arr.size):
        i = idx[k]
        zb, Tb, Lb, Pb = _ZB[i], _TB[i], _LB[i], _PB[i]
        dz = z_arr[k] - zb
        T[k] = Tb + Lb * dz
        if abs(Lb) > 1e-12:
            P[k] = Pb * (Tb / T[k]) ** (G0 / (RD * Lb))
        else:
            P[k] = Pb * np.exp(-G0 * dz / (RD * Tb))
    if np.isscalar(z) or np.ndim(z) == 0:
        return float(T[0]), float(P[0])
    return T, P


def load_atm_profile(path, height_col=None, temp_col=None, press_col=None,
                      temp_units="K", press_units="Pa", delimiter=","):
    """Load a user-supplied atmospheric sounding (e.g. radiosonde/reanalysis).

    The file should be a text/CSV table with (at minimum) columns for height,
    temperature and pressure. Column names are auto-detected (case-insensitive,
    matching 'height'/'alt'/'z', 'temp'/'t', 'press'/'p') unless explicitly
    given via height_col/temp_col/press_col.

    Parameters
    ----------
    path : str
        Path to the profile file.
    height_col, temp_col, press_col : str, optional
        Explicit column names to use.
    temp_units : {'K','C'}
        Units of the temperature column in the file.
    press_units : {'Pa','hPa'}
        Units of the pressure column in the file.
    delimiter : str
        Column delimiter (default ',').

    Returns
    -------
    profile : dict with keys 'z' [m], 'T' [K], 'P' [Pa], sorted by ascending z.
    """
    import csv

    with open(path, "r", newline="") as f:
        reader = csv.reader(f, delimiter=delimiter)
        rows = [r for r in reader if any(cell.strip() for cell in r)]
    header = [h.strip().lower() for h in rows[0]]
    data_rows = rows[1:]

    def find_col(explicit, descriptive_keywords, exact_shorthand):
        # Matches in two passes to avoid false positives from short,
        # single-letter shorthand (e.g. a bare 't' or 'p' column keyword
        # must NOT match as a substring of an unrelated column like
        # 'height_m', which happens to contain the letter 't').
        if explicit is not None:
            return header.index(explicit.strip().lower())
        for i, h in enumerate(header):
            if any(k in h for k in descriptive_keywords):
                return i
        for i, h in enumerate(header):
            if h in exact_shorthand:
                return i
        raise ValueError(
            f"Could not find a column matching {descriptive_keywords} "
            f"(or exactly named {sorted(exact_shorthand)}) in header {header}"
        )

    iz = find_col(height_col, ["height", "altitude", "alt"], {"z"})
    it = find_col(temp_col, ["temp"], {"t"})
    ip = find_col(press_col, ["press"], {"p"})

    z = np.array([float(r[iz]) for r in data_rows])
    T = np.array([float(r[it]) for r in data_rows])
    P = np.array([float(r[ip]) for r in data_rows])

    if temp_units.upper() == "C":
        T = T + 273.15
    if press_units.lower() == "hpa":
        P = P * 100.0

    order = np.argsort(z)
    return {"z": z[order], "T": T[order], "P": P[order]}


def atmosphere_TP(z_asl, mode="standard", profile=None, T_offset=0.0):
    """Return ambient temperature and pressure at altitude(s) z_asl [m a.s.l.].

    Parameters
    ----------
    z_asl : scalar or array_like
        Altitude above sea level [m].
    mode : {'standard','profile'}
        'standard' uses the built-in ICAO/US Standard Atmosphere 1976.
        'profile' interpolates the user-supplied `profile` dict
        (see load_atm_profile). Linear extrapolation is used outside range.
    profile : dict, required if mode='profile'
        {'z':..., 'T':..., 'P':...} in SI units (m, K, Pa).
    T_offset : float
        Optional uniform additive temperature offset [K] applied to the
        standard atmosphere only -- handy for shifting the standard profile
        to match a known local/near-vent surface temperature without needing
        a full custom sounding. Ignored when mode='profile'.

    Returns
    -------
    T [K], P [Pa] (scalars if z_asl is scalar, else ndarrays)
    """
    z_arr = np.atleast_1d(np.asarray(z_asl, dtype=float))
    if mode == "profile":
        if profile is None or "z" not in profile:
            raise ValueError(
                "mode='profile' requires a profile dict with 'z','T','P' "
                "(see load_atm_profile). Use mode='standard' otherwise."
            )
        T = np.interp(z_arr, profile["z"], profile["T"])
        P = np.interp(z_arr, profile["z"], profile["P"])
    else:
        T, P = standard_atmosphere(z_arr)
        T = np.atleast_1d(T) + T_offset
        P = np.atleast_1d(P)
    if np.isscalar(z_asl) or np.ndim(z_asl) == 0:
        return float(T[0]), float(P[0])
    return T, P


def potential_temperature(T, P):
    """Dry potential temperature theta = T * (P0ref/P)^kappa."""
    return T * (P0REF / P) ** KAPPA


def air_density(T, z_asl, mode="standard", profile=None, T_offset=0.0):
    """Ideal-gas density [kg/m^3] of a dry-air-like parcel at temperature T
    [K], assumed to be in pressure equilibrium with the local ambient
    atmosphere at altitude z_asl [m a.s.l.].

    rho = P_ambient(z_asl) / (Rd * T)

    This mirrors the "AirDensityFunc(T,z)" convention used elsewhere in this
    project, generalized to take the ambient pressure from a full T(z),P(z)
    atmospheric profile (standard or user-supplied) instead of a fixed value.
    """
    _, P = atmosphere_TP(z_asl, mode, profile, T_offset)
    rho = P / (RD * T)
    return rho, P


def brunt_vaisala_N2(z_asl, mode="standard", profile=None, T_offset=0.0, dz=1.0):
    """Real, height-dependent Brunt-Vaisala frequency squared of the ambient
    atmosphere, N^2(z) = (g/theta_a) d(theta_a)/dz, evaluated at altitude
    z_asl [m a.s.l.] using a centered finite difference of step `dz` [m].

    Using potential temperature (rather than a raw density difference)
    correctly generalizes the constant-N^2 Boussinesq stratification used in
    the tank-experiment model to a real, compressible atmosphere: it
    automatically accounts for the adiabatic expansion/cooling that a rising
    ring undergoes purely from moving to lower ambient pressure, separate
    from any mixing with entrained air.
    """
    Tm, Pm = atmosphere_TP(z_asl, mode, profile, T_offset)
    T1, P1 = atmosphere_TP(np.asarray(z_asl) - dz, mode, profile, T_offset)
    T2, P2 = atmosphere_TP(np.asarray(z_asl) + dz, mode, profile, T_offset)
    theta_m = potential_temperature(Tm, Pm)
    theta1 = potential_temperature(T1, P1)
    theta2 = potential_temperature(T2, P2)
    dtheta_dz = (theta2 - theta1) / (2.0 * dz)
    return G0 / theta_m * dtheta_dz
