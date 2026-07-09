VOLCANIC VORTEX RING RISE MODEL
================================
1D variable-entrainment vortex-ring / thermal rise model for natural
volcanic vortex rings rising through a real, stratified atmosphere.
Available as both a MATLAB version and a Python version (the Python
version also has an optional desktop GUI).

Author: Johan Gilchrist


FOLDER CONTENTS
---------------
matlab/
    run_volcanic_vortex_ring_model.m   Main script -- edit parameters at
                                        the top and run.
    vortex_ring_ode_volcanic.m         ODE right-hand side (height-domain).
    apex_event.m                       Stops integration at the rise apex.
    atmosphere_TP.m                    Ambient T(z), P(z) dispatcher.
    standard_atmosphere.m              Built-in ICAO/US Standard Atmosphere 1976.
    load_atm_profile.m                 Loads a user-supplied sounding file.
    potential_temperature.m            Potential-temperature helper.
    air_density.m                      Ideal-gas density helper.
    brunt_vaisala_N2.m                 Height-dependent stratification N^2(z).

python/
    run_volcanic_vortex_ring_model.py  Command-line script -- edit
                                        parameters at the top and run.
    vortex_ring_gui.py                 Desktop GUI (Tkinter) -- set
                                        parameters in a form instead of
                                        editing the script.
    model_core.py                      Shared model-driver logic used by
                                        both of the above.
    vortex_ring_ode.py                 ODE right-hand side (height-domain).
    atmosphere.py                      Atmosphere helpers (standard
                                        atmosphere, profile loading,
                                        potential temperature, N^2).
    requirements.txt                   Python package dependencies.

web/
    index.html                         Browser version of the model + GUI.
                                        Open this file (or host the whole
                                        web/ folder) to run it -- no
                                        installation, no server required.
    style.css                          Page styling.
    app.js                             Front-end logic (loads Pyodide,
                                        runs the model, renders plots).
    python/                            Same model physics as python/ above,
                                        adapted to run in-browser via
                                        Pyodide (WebAssembly) instead of
                                        scipy.

README.txt                             This file.


WHAT YOU NEED
-------------
- MATLAB: no toolboxes required beyond base MATLAB (uses only ode23s,
  readtable, inputParser -- all standard). R2016b or newer (for local
  functions in load_atm_profile.m).
- Python: version 3.9 or newer, plus the packages listed in
  python/requirements.txt (numpy, scipy, matplotlib). The GUI additionally
  needs Tkinter, which normally comes bundled with Python itself (see
  troubleshooting below if it's missing or fails to start).


===============================================================================
MATLAB -- SETUP AND USE (Mac and Windows, same steps)
===============================================================================
1. Open MATLAB.
2. In the MATLAB "Current Folder" browser (or via `cd`), navigate to the
   matlab/ folder from this package.
3. Open run_volcanic_vortex_ring_model.m.
4. Edit the "USER-DEFINED PARAMETERS" section near the top (vent altitude,
   radius, exit velocity, exit temperature, ambient temperature,
   entrainment settings, etc.) -- values are commented with their meaning
   and typical ranges.
5. Click Run (or press F5). Three figure windows will appear (fluxes;
   velocity/radius/g'/temperature profiles; entrainment and Richardson
   number profiles), and diagnostics print to the Command Window.


===============================================================================
PYTHON -- SETUP (Mac)
===============================================================================
These steps create an isolated Python environment ("virtual environment")
so this package's dependencies don't interfere with anything else on your
system.

1. Open Terminal.

2. Check you have Python 3.9+ available:
       python3 --version

3. (Recommended) Confirm that Python's Tkinter works BEFORE creating the
   environment, since a virtual environment inherits whatever Tk support
   the Python you use to create it has:
       python3 -m tkinter
   A small test window should pop up. If it errors out (e.g. "no display
   name and no $DISPLAY environment variable"), see TROUBLESHOOTING below
   -- try a different python3 (e.g. the one from python.org, or Homebrew's
   `python-tk`) before continuing.

4. Create and activate a virtual environment (from inside the python/
   folder of this package):
       cd path/to/VolcanicVortexRingModel/python
       python3 -m venv vvr-env
       source vvr-env/bin/activate

5. Install the required packages:
       pip install -r requirements.txt

6. Run the GUI:
       python vortex_ring_gui.py

   ...or the command-line script:
       python run_volcanic_vortex_ring_model.py

7. Next time, you only need to reactivate the environment (step 4's
   `source vvr-env/bin/activate`) before running the scripts again --
   no need to recreate it or reinstall packages.


===============================================================================
PYTHON -- SETUP (Windows)
===============================================================================
1. Install Python from https://www.python.org/downloads/ (the official
   installer). During install, leave "tcl/tk and IDLE" checked (it is
   checked by default) so Tkinter is included, and check "Add python.exe
   to PATH".

2. Open Command Prompt (or PowerShell).

3. Create and activate a virtual environment (from inside the python\
   folder of this package):
       cd path\to\VolcanicVortexRingModel\python
       python -m venv vvr-env

   Command Prompt:
       vvr-env\Scripts\activate.bat
   PowerShell:
       vvr-env\Scripts\Activate.ps1
   (If PowerShell blocks the activation script with an "execution policy"
   error, either use Command Prompt instead, or run PowerShell as
   Administrator once and execute:
       Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
   then retry.)

4. Install the required packages:
       pip install -r requirements.txt

5. Run the GUI:
       python vortex_ring_gui.py

   ...or the command-line script:
       python run_volcanic_vortex_ring_model.py

6. Next time, you only need to reactivate the environment (step 3's
   activate command) before running the scripts again.

Note: Tkinter is native on Windows (it doesn't rely on X11 or any
external display server the way it can on Mac/Linux), so Windows users are
generally less likely to hit GUI-startup problems than Mac users are.


===============================================================================
WEB VERSION -- HOSTING ON YOUR WEBSITE (johangilchrist.com)
===============================================================================
The web/ folder is a fully self-contained, static version of the model and
GUI that runs entirely in the visitor's browser -- there is no server-side
code and nothing is sent anywhere (all computation happens locally on the
visitor's machine via Pyodide, which runs Python compiled to WebAssembly).

To deploy it:
1. Upload the entire web/ folder as-is to your static host, keeping the
   python/ subfolder alongside index.html (the page fetches those .py
   files at load time, so the relative folder structure must be
   preserved).
2. Link to web/index.html from your site's navigation, or embed it in an
   existing page with an iframe, e.g.:
       <iframe src="/VolcanicVortexRingModel/web/index.html"
               style="width:100%; height:1400px; border:none;"></iframe>
3. That's it -- no build step, no database, no backend. Any static file
   host (the same one already serving your HTML site) works.

Notes:
- First load takes a few seconds while the browser downloads and
  initializes Pyodide + numpy (from a CDN); results appear instantly on
  every run after that.
- Requires a reasonably modern browser (Chrome, Firefox, Safari, Edge --
  anything from the last ~3 years). No mobile-specific testing has been
  done, though the layout is responsive down to narrow screens.
- The web version has the same physics and produces the same results as
  the desktop Python version -- it was validated against it to 4 decimal
  places before packaging.


===============================================================================
TROUBLESHOOTING
===============================================================================

"_tkinter.TclError: no display name and no $DISPLAY environment variable"
(Mac only)
    This means the Python you're running was built against an X11 Tk
    instead of native macOS Quartz/Cocoa Tk (this is common with some
    MacPorts Python/Tk installs, but not usually python.org or Homebrew
    installs) -- it's trying to open an X11 server that isn't running.
    Fixes, in order of ease:
      1. Use a different Python to create your virtual environment --
         e.g. Python from python.org, or Homebrew's python
         (`brew install python-tk` to make sure Tk is included).
      2. Install and launch XQuartz (https://www.xquartz.org), then retry.
      3. If you specifically need MacPorts Python, reinstall Tk with the
         Quartz variant: `sudo port install tk +quartz`

ImportError / dlopen error mentioning "libgfortran" when importing scipy
    This means the Python environment's scipy/numpy install is broken or
    mismatched (common in Anaconda "base" environments that have drifted
    over time or mixed installation channels). This is NOT a problem with
    these scripts. Fix: create a fresh virtual environment as described
    above (via `venv`, not conda) and `pip install -r requirements.txt`
    into it -- a clean, single-source install avoids the mismatch. If you
    are using conda instead of venv and don't want to abandon it, create a
    brand new conda environment rather than trying to repair an existing
    one:
        conda create -n vvr python=3.12 numpy scipy matplotlib tk -c conda-forge -y
        conda activate vvr

"Solving environment" hangs forever in conda
    This is a known slow point with conda's default solver, especially in
    large/old "base" environments or when mixing the `defaults` and
    `conda-forge` channels. Prefer `venv` (see above) for this package, or
    if using conda, always create a small, fresh environment from a single
    channel rather than modifying an existing large one.

GUI runs but the model integration seems to fail or produce NaN/blow up
    This can happen with the "variable" (height-dependent) entrainment
    mode if the source conditions imply an extremely large buoyancy
    contrast relative to the ambient air (e.g. a very hot, very slow
    ring). The scripts already cap the entrainment coefficient
    (alpha_max, default 1.0) as a safeguard -- if you still see problems,
    try increasing the exit velocity, lowering alpha_max, or switching to
    "constant" entrainment mode with a modest alpha (e.g. 0.1-0.2).


===============================================================================
QUICK PARAMETER REFERENCE
===============================================================================
- Vent altitude (m above sea level), ring radius (m), exit velocity (m/s)
  and exit temperature (degC) are all set either at the top of the
  scripts (edit directly) or via the GUI's input fields.
- Ambient conditions: use the built-in standard atmosphere (optionally
  shifted to match a known near-vent temperature) or load a real sounding
  file (CSV with height/temperature/pressure columns) via the GUI's
  "Load profile file..." option or load_atm_profile()/load_atm_profile.m.
- Entrainment: "variable" (height-dependent, the model's default closure)
  or "constant" (a single value you specify) -- both available from the
  GUI and as script parameters.
