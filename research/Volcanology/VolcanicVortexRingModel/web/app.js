/*
 * app.js
 * ========
 * Front-end logic for the browser version of the volcanic vortex ring
 * rise model. Loads Pyodide (Python compiled to WebAssembly), fetches the
 * model's Python source files, and calls into web_api.run_from_web() --
 * the exact same physics as the desktop MATLAB/Python package -- entirely
 * client-side. Plots are rendered with Plotly.js.
 *
 * Johan Gilchrist
 */

const PY_FILES = ["atmosphere.py", "simple_ode.py", "vortex_ring_ode.py", "model_core.py", "web_api.py"];
const UPLOADED_PROFILE_PATH = "/home/pyodide/uploaded_profile.csv";

let pyodide = null;
let runFromWeb = null;
let uploadedProfileReady = false;
let lastResult = null;

const $ = (id) => document.getElementById(id);

function radioValue(name) {
  const el = document.querySelector(`input[name="${name}"]:checked`);
  return el ? el.value : null;
}

function setStatus(msg, isError) {
  const el = $("run-status");
  el.textContent = msg || "";
  el.style.color = isError ? "#b3261e" : "";
}

// ---------------------------------------------------------------------------
// Pyodide bootstrap
// ---------------------------------------------------------------------------
async function init() {
  try {
    pyodide = await loadPyodide();
    await pyodide.loadPackage(["numpy"]);

    for (const filename of PY_FILES) {
      const resp = await fetch(`python/${filename}`);
      if (!resp.ok) {
        throw new Error(`Failed to fetch python/${filename} (HTTP ${resp.status})`);
      }
      const src = await resp.text();
      pyodide.FS.writeFile(`/home/pyodide/${filename}`, src);
    }

    pyodide.runPython(`
import sys
if "/home/pyodide" not in sys.path:
    sys.path.insert(0, "/home/pyodide")
from web_api import run_from_web
`);
    runFromWeb = pyodide.globals.get("run_from_web");

    $("loading-banner").classList.add("hidden");
    $("app").classList.remove("hidden");
    $("run-btn").disabled = false;
    $("run-btn").textContent = "Run model";
  } catch (err) {
    $("loading-banner").textContent =
      "Failed to load the Python environment: " + err.message +
      ". Try reloading the page; if the problem persists, your browser or " +
      "network connection may be blocking the Pyodide/CDN scripts.";
    console.error(err);
  }
}

// ---------------------------------------------------------------------------
// Form show/hide logic
// ---------------------------------------------------------------------------
function updateAtmModeUI() {
  const mode = radioValue("atm_mode");
  $("standard-atm-fields").style.display = mode === "standard" ? "" : "none";
  $("profile-atm-fields").style.display = mode === "profile" ? "" : "none";
}

function updateEntrainModeUI() {
  const mode = radioValue("entrainment_mode");
  $("variable-entrain-fields").style.display = mode === "variable" ? "" : "none";
  $("constant-entrain-fields").style.display = mode === "constant" ? "" : "none";
}

document.querySelectorAll('input[name="atm_mode"]').forEach((r) =>
  r.addEventListener("change", updateAtmModeUI)
);
document.querySelectorAll('input[name="entrainment_mode"]').forEach((r) =>
  r.addEventListener("change", updateEntrainModeUI)
);

// ---------------------------------------------------------------------------
// Profile file upload
// ---------------------------------------------------------------------------
$("profile-file").addEventListener("change", async (evt) => {
  const file = evt.target.files[0];
  uploadedProfileReady = false;
  if (!file) {
    $("profile-file-status").textContent = "No file selected.";
    return;
  }
  try {
    const text = await file.text();
    if (pyodide) {
      pyodide.FS.writeFile(UPLOADED_PROFILE_PATH, text);
    }
    uploadedProfileReady = true;
    $("profile-file-status").textContent = `Loaded "${file.name}" (${file.size} bytes).`;
  } catch (err) {
    $("profile-file-status").textContent = "Failed to read file: " + err.message;
  }
});

// ---------------------------------------------------------------------------
// Run the model
// ---------------------------------------------------------------------------
function collectParams() {
  const num = (id) => {
    const v = parseFloat($(id).value);
    return Number.isFinite(v) ? v : NaN;
  };

  const params = {
    z0_asl: num("z0_asl"),
    r0: num("r0"),
    u0: num("u0"),
    Tj0_C: num("Tj0_C"),
    atm_mode: radioValue("atm_mode"),
    Ta0_target_C: num("Ta0_target_C"),
    entrainment_mode: radioValue("entrainment_mode"),
    c_coef: num("c_coef"),
    alpha_max: num("alpha_max"),
    alpha_const: num("alpha_const"),
    zmax_rise: num("zmax_rise"),
  };

  if (params.atm_mode === "profile") {
    params.profile_path = UPLOADED_PROFILE_PATH;
    params.profile_temp_units = $("profile_temp_units").value;
    params.profile_press_units = $("profile_press_units").value;
  }

  return params;
}

function validateParams(params) {
  const requiredNumeric = ["z0_asl", "r0", "u0", "Tj0_C", "zmax_rise"];
  for (const key of requiredNumeric) {
    if (Number.isNaN(params[key])) {
      return `Please enter a valid number for "${key}".`;
    }
  }
  if (params.r0 <= 0) return "Vent radius must be greater than 0.";
  if (params.u0 <= 0) return "Initial velocity must be greater than 0.";
  if (params.zmax_rise <= 2 * params.r0) {
    return "Max height above vent must be greater than 2 x vent radius.";
  }
  if (params.atm_mode === "standard" && Number.isNaN(params.Ta0_target_C)) {
    return 'Please enter a valid "Ambient temperature at vent altitude".';
  }
  if (params.atm_mode === "profile" && !uploadedProfileReady) {
    return "Please choose an atmospheric profile file, or switch to Standard atmosphere.";
  }
  if (params.entrainment_mode === "variable") {
    if (Number.isNaN(params.c_coef) || params.c_coef <= 0) return "c coefficient must be a positive number.";
    if (Number.isNaN(params.alpha_max) || params.alpha_max <= 0) return "alpha cap must be a positive number.";
  } else {
    if (Number.isNaN(params.alpha_const) || params.alpha_const <= 0) return "Constant alpha must be a positive number.";
  }
  return null;
}

$("model-form").addEventListener("submit", async (evt) => {
  evt.preventDefault();
  if (!runFromWeb) {
    setStatus("Python environment is still loading -- please wait.", true);
    return;
  }

  const params = collectParams();
  const validationError = validateParams(params);
  if (validationError) {
    setStatus(validationError, true);
    return;
  }

  const runBtn = $("run-btn");
  runBtn.disabled = true;
  runBtn.textContent = "Running…";
  setStatus("");

  try {
    const resultJsonStr = runFromWeb(JSON.stringify(params));
    const result = JSON.parse(resultJsonStr);

    if (result.error) {
      setStatus("Model error: " + result.error, true);
      return;
    }

    lastResult = result;
    renderReport(result);
    renderPlots(result);
    $("download-csv-btn").disabled = false;
    setStatus("Done.");
  } catch (err) {
    setStatus("Unexpected error: " + err.message, true);
    console.error(err);
  } finally {
    runBtn.disabled = false;
    runBtn.textContent = "Run model";
  }
});

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------
function renderReport(res) {
  $("report-text").textContent = res.report_text;
}

function makeLayout(xTitle, extra) {
  return Object.assign(
    {
      margin: { t: 10, r: 15, b: 45, l: 55 },
      xaxis: { title: xTitle },
      yaxis: { title: "Height above vent [m]" },
      font: { size: 11 },
    },
    extra || {}
  );
}

function scatterTrace(x, y, extra) {
  return Object.assign(
    { x, y, mode: "markers", type: "scatter", marker: { size: 3 } },
    extra || {}
  );
}

const PLOTLY_CONFIG = { responsive: true, displaylogo: false };

function renderPlots(res) {
  const z = res.z_out;

  Plotly.newPlot("plot-V", [scatterTrace(res.V, z)], makeLayout("V [m^3]"), PLOTLY_CONFIG);
  Plotly.newPlot("plot-M", [scatterTrace(res.M, z)], makeLayout("M [m^4/s]"), PLOTLY_CONFIG);
  Plotly.newPlot("plot-B", [scatterTrace(res.B, z)], makeLayout("B [m^4/s^2]"), PLOTLY_CONFIG);

  Plotly.newPlot("plot-u", [scatterTrace(res.u, z)], makeLayout("Vertical velocity [m/s]"), PLOTLY_CONFIG);
  Plotly.newPlot("plot-r", [scatterTrace(res.r, z)], makeLayout("Radius [m]"), PLOTLY_CONFIG);
  Plotly.newPlot("plot-gprime", [scatterTrace(res.gprime, z)], makeLayout("g' [m/s^2]"), PLOTLY_CONFIG);

  const Tring = res.Tp_z.map((t) => t - 273.15);
  const Tamb = res.Ta_z.map((t) => t - 273.15);
  Plotly.newPlot(
    "plot-T",
    [
      scatterTrace(Tring, z, { name: "ring" }),
      { x: Tamb, y: z, mode: "lines", name: "ambient", line: { dash: "dash" } },
    ],
    makeLayout("Temperature [deg C]", { showlegend: true, legend: { x: 0.55, y: 1 } }),
    PLOTLY_CONFIG
  );

  Plotly.newPlot(
    "plot-alpha",
    [scatterTrace(res.alpha_k_z, z)],
    makeLayout("alpha_e", { xaxis: { title: "alpha_e", range: [0, 1] } }),
    PLOTLY_CONFIG
  );
  Plotly.newPlot(
    "plot-Ri",
    [scatterTrace(res.Ri_z, z)],
    makeLayout("Richardson number, Ri", { xaxis: { title: "Richardson number, Ri", type: "log" } }),
    PLOTLY_CONFIG
  );
}

// ---------------------------------------------------------------------------
// CSV export
// ---------------------------------------------------------------------------
$("download-csv-btn").addEventListener("click", () => {
  if (!lastResult) return;
  const res = lastResult;
  const header = [
    "z_above_vent_m", "V_m3", "M_m4s", "B_m4s2", "alpha_e", "Ri",
    "u_ms", "r_m", "gprime_ms2", "Tring_K", "Tambient_K", "Pambient_Pa",
  ];
  const rows = [header.join(",")];
  for (let i = 0; i < res.z_out.length; i++) {
    rows.push([
      res.z_out[i], res.V[i], res.M[i], res.B[i], res.alpha_k_z[i], res.Ri_z[i],
      res.u[i], res.r[i], res.gprime[i], res.Tp_z[i], res.Ta_z[i], res.Pa_z[i],
    ].join(","));
  }
  const blob = new Blob([rows.join("\n")], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "vortex_ring_results.csv";
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
});

// ---------------------------------------------------------------------------
// Boot
// ---------------------------------------------------------------------------
updateAtmModeUI();
updateEntrainModeUI();
init();
