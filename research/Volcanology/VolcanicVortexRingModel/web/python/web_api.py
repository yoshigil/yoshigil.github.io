"""
web_api.py
============
Thin JSON-in / JSON-out wrapper around model_core, meant to be called
directly from JavaScript (via a Pyodide PyProxy function call with a single
string argument) rather than exercising Pyodide's more elaborate JS<->Python
object-conversion machinery. Keeping the browser/Python boundary to plain
strings keeps the front-end code simple and avoids edge cases in converting
nested dicts/None/numpy types across the JS<->Python boundary.

Johan Gilchrist
"""

import json

from model_core import run_vortex_ring_model, format_report, ModelError
from atmosphere import load_atm_profile


def run_from_web(params_json):
    """params_json: a JSON string with keys matching
    model_core.run_vortex_ring_model()'s parameters (plus 'profile_path',
    'profile_temp_units', 'profile_press_units' when atm_mode='profile',
    where profile_path is a path already written into Pyodide's virtual
    filesystem by the caller, e.g. via pyodide.FS.writeFile()).

    Returns a JSON string: either the run_vortex_ring_model() result dict
    (with an added 'report_text' field), or {"error": "..."} on failure.
    """
    try:
        p = json.loads(params_json)

        atm_mode = p.get("atm_mode", "standard")
        atm_profile = None
        if atm_mode == "profile":
            profile_path = p.get("profile_path")
            if not profile_path:
                raise ModelError("No atmospheric profile file was provided.")
            atm_profile = load_atm_profile(
                profile_path,
                temp_units=p.get("profile_temp_units", "K"),
                press_units=p.get("profile_press_units", "Pa"),
            )

        res = run_vortex_ring_model(
            z0_asl=float(p["z0_asl"]),
            r0=float(p["r0"]),
            u0=float(p["u0"]),
            Tj0_C=float(p["Tj0_C"]),
            atm_mode=atm_mode,
            atm_profile=atm_profile,
            Ta0_target_C=float(p.get("Ta0_target_C", 15.0)),
            entrainment_mode=p.get("entrainment_mode", "variable"),
            c_coef=float(p.get("c_coef", 1.0)),
            alpha_max=float(p.get("alpha_max", 1.0)),
            alpha_const=float(p.get("alpha_const", 0.1)),
            zmax_rise=float(p.get("zmax_rise", 3000.0)),
        )
        res["report_text"] = format_report(res)
        return json.dumps(res)
    except Exception as exc:  # noqa: BLE001 -- deliberately broad: any
        # failure here should surface as a readable message in the browser
        # UI rather than an unhandled JS-side exception.
        return json.dumps({"error": str(exc)})
