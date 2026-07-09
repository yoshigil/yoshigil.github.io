"""
vortex_ring_gui.py
====================
Desktop GUI (Tkinter) for the 1D volcanic vortex ring rise model.

Lets the user set the source/vent conditions (radius, velocity,
temperature, vent altitude), choose between the built-in standard
atmosphere or a loaded atmospheric-profile file, and choose between the
model's variable (height-dependent) entrainment closure or a single
user-specified constant entrainment coefficient. Results are shown as a
text diagnostic summary plus two tabs of embedded plots, and can be
exported to CSV/PNG.

Run with:   python3 vortex_ring_gui.py

All of the physics/integration lives in model_core.py (shared with the
command-line script run_volcanic_vortex_ring_model.py), so this file is
just the front-end.

Johan Gilchrist
"""

import traceback
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from model_core import run_vortex_ring_model, format_report, ModelError
from atmosphere import load_atm_profile


class VortexRingGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Volcanic Vortex Ring Rise Model")
        self.geometry("1200x750")

        self.atm_profile = None       # loaded profile dict, once a file is picked
        self.atm_profile_path = None  # display string
        self.last_result = None       # most recent run_vortex_ring_model() output

        self._build_widgets()
        self._on_atm_mode_change()
        self._on_entrainment_mode_change()

    # ------------------------------------------------------------------
    # Widget construction
    # ------------------------------------------------------------------
    def _build_widgets(self):
        root = ttk.Frame(self, padding=8)
        root.pack(side=tk.LEFT, fill=tk.Y)

        plot_frame = ttk.Frame(self, padding=4)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self._build_input_panel(root)
        self._build_plot_panel(plot_frame)

    def _build_input_panel(self, parent):
        # --- Source (vent) conditions ---------------------------------
        src = ttk.LabelFrame(parent, text="Source (vent) conditions", padding=8)
        src.pack(fill=tk.X, pady=4)

        self.var_z0 = tk.StringVar(value="3380")
        self.var_r0 = tk.StringVar(value="5")
        self.var_u0 = tk.StringVar(value="20")
        self.var_Tj0 = tk.StringVar(value="450")

        self._labeled_entry(src, "Vent altitude a.s.l. [m]:", self.var_z0, 0)
        self._labeled_entry(src, "Vent/initial ring radius [m]:", self.var_r0, 1)
        self._labeled_entry(src, "Initial (exit) velocity [m/s]:", self.var_u0, 2)
        self._labeled_entry(src, "Initial ring temperature [degC]:", self.var_Tj0, 3)

        # --- Atmosphere --------------------------------------------------
        atm = ttk.LabelFrame(parent, text="Atmosphere", padding=8)
        atm.pack(fill=tk.X, pady=4)

        self.var_atm_mode = tk.StringVar(value="standard")
        ttk.Radiobutton(atm, text="Standard atmosphere (ICAO/US 1976)",
                         variable=self.var_atm_mode, value="standard",
                         command=self._on_atm_mode_change).grid(
            row=0, column=0, columnspan=3, sticky="w")
        ttk.Radiobutton(atm, text="Load profile file...",
                         variable=self.var_atm_mode, value="profile",
                         command=self._on_atm_mode_change).grid(
            row=1, column=0, columnspan=3, sticky="w")

        self.var_Ta0 = tk.StringVar(value="15")
        self.lbl_Ta0 = ttk.Label(atm, text="Ambient T at vent altitude [degC]:")
        self.lbl_Ta0.grid(row=2, column=0, sticky="w", pady=2)
        self.entry_Ta0 = ttk.Entry(atm, textvariable=self.var_Ta0, width=12)
        self.entry_Ta0.grid(row=2, column=1, sticky="w", pady=2)

        self.btn_browse = ttk.Button(atm, text="Browse...", command=self._browse_profile)
        self.btn_browse.grid(row=3, column=0, sticky="w", pady=2)
        self.lbl_profile_path = ttk.Label(atm, text="(no file selected)", foreground="gray")
        self.lbl_profile_path.grid(row=3, column=1, columnspan=2, sticky="w", pady=2)

        self.lbl_tunits = ttk.Label(atm, text="Profile temp. units:")
        self.lbl_tunits.grid(row=4, column=0, sticky="w")
        self.var_tunits = tk.StringVar(value="K")
        self.combo_tunits = ttk.Combobox(atm, textvariable=self.var_tunits,
                                          values=["K", "C"], width=5, state="readonly")
        self.combo_tunits.grid(row=4, column=1, sticky="w")

        self.lbl_punits = ttk.Label(atm, text="Profile press. units:")
        self.lbl_punits.grid(row=5, column=0, sticky="w")
        self.var_punits = tk.StringVar(value="Pa")
        self.combo_punits = ttk.Combobox(atm, textvariable=self.var_punits,
                                          values=["Pa", "hPa"], width=5, state="readonly")
        self.combo_punits.grid(row=5, column=1, sticky="w")

        # --- Entrainment ---------------------------------------------------
        ent = ttk.LabelFrame(parent, text="Entrainment", padding=8)
        ent.pack(fill=tk.X, pady=4)

        self.var_ent_mode = tk.StringVar(value="variable")
        ttk.Radiobutton(ent, text="Variable (height-dependent, Ri-based closure)",
                         variable=self.var_ent_mode, value="variable",
                         command=self._on_entrainment_mode_change).grid(
            row=0, column=0, columnspan=3, sticky="w")
        ttk.Radiobutton(ent, text="Constant (single user-specified value)",
                         variable=self.var_ent_mode, value="constant",
                         command=self._on_entrainment_mode_change).grid(
            row=1, column=0, columnspan=3, sticky="w")

        self.var_c = tk.StringVar(value="1.0")
        self.var_alpha_max = tk.StringVar(value="1.0")
        self.var_alpha_const = tk.StringVar(value="0.1")

        self.lbl_c = ttk.Label(ent, text="c coefficient (Turner 1962):")
        self.lbl_c.grid(row=2, column=0, sticky="w", pady=2)
        self.entry_c = ttk.Entry(ent, textvariable=self.var_c, width=10)
        self.entry_c.grid(row=2, column=1, sticky="w", pady=2)

        self.lbl_amax = ttk.Label(ent, text="alpha cap (alpha_max):")
        self.lbl_amax.grid(row=3, column=0, sticky="w", pady=2)
        self.entry_amax = ttk.Entry(ent, textvariable=self.var_alpha_max, width=10)
        self.entry_amax.grid(row=3, column=1, sticky="w", pady=2)

        self.lbl_aconst = ttk.Label(ent, text="constant alpha:")
        self.lbl_aconst.grid(row=4, column=0, sticky="w", pady=2)
        self.entry_aconst = ttk.Entry(ent, textvariable=self.var_alpha_const, width=10)
        self.entry_aconst.grid(row=4, column=1, sticky="w", pady=2)

        # --- Integration range ----------------------------------------
        integ = ttk.LabelFrame(parent, text="Integration", padding=8)
        integ.pack(fill=tk.X, pady=4)
        self.var_zmax = tk.StringVar(value="3000")
        self._labeled_entry(integ, "Max height above vent to try [m]:", self.var_zmax, 0)

        # --- Actions -------------------------------------------------------
        actions = ttk.Frame(parent, padding=4)
        actions.pack(fill=tk.X, pady=8)
        self.btn_run = ttk.Button(actions, text="Run model", command=self._on_run)
        self.btn_run.pack(side=tk.LEFT, padx=2)
        self.btn_save_csv = ttk.Button(actions, text="Save results (CSV)...",
                                        command=self._on_save_csv, state="disabled")
        self.btn_save_csv.pack(side=tk.LEFT, padx=2)
        self.btn_save_png = ttk.Button(actions, text="Save plots (PNG)...",
                                        command=self._on_save_png, state="disabled")
        self.btn_save_png.pack(side=tk.LEFT, padx=2)

        self.status_var = tk.StringVar(value="Ready.")
        ttk.Label(parent, textvariable=self.status_var, foreground="blue").pack(
            fill=tk.X, pady=(0, 4))

        # --- Diagnostics text -----------------------------------------
        diag = ttk.LabelFrame(parent, text="Diagnostics", padding=4)
        diag.pack(fill=tk.BOTH, expand=True, pady=4)
        self.txt_diag = scrolledtext.ScrolledText(diag, width=48, height=20,
                                                    font=("Courier New", 9))
        self.txt_diag.pack(fill=tk.BOTH, expand=True)

    def _labeled_entry(self, parent, label, var, row):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=2)
        ttk.Entry(parent, textvariable=var, width=12).grid(
            row=row, column=1, sticky="w", pady=2)

    def _build_plot_panel(self, parent):
        self.notebook = ttk.Notebook(parent)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # Tab 1: fluxes (Volume, Momentum, Buoyancy)
        self.fig_fluxes = plt.Figure(figsize=(10, 5))
        self.axs_fluxes = self.fig_fluxes.subplots(1, 3)
        tab1 = ttk.Frame(self.notebook)
        self.notebook.add(tab1, text="Fluxes")
        self.canvas_fluxes = FigureCanvasTkAgg(self.fig_fluxes, master=tab1)
        self.canvas_fluxes.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas_fluxes, tab1).pack(fill=tk.X)

        # Tab 2: profiles (velocity, radius, reduced gravity, temperature)
        self.fig_profiles = plt.Figure(figsize=(10, 5))
        self.axs_profiles = self.fig_profiles.subplots(1, 4)
        tab2 = ttk.Frame(self.notebook)
        self.notebook.add(tab2, text="Profiles")
        self.canvas_profiles = FigureCanvasTkAgg(self.fig_profiles, master=tab2)
        self.canvas_profiles.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas_profiles, tab2).pack(fill=tk.X)

        # Tab 3: entrainment coefficient and Richardson number
        self.fig_entrain = plt.Figure(figsize=(6, 5))
        self.axs_entrain = self.fig_entrain.subplots(1, 2)
        tab3 = ttk.Frame(self.notebook)
        self.notebook.add(tab3, text="Entrainment / Ri")
        self.canvas_entrain = FigureCanvasTkAgg(self.fig_entrain, master=tab3)
        self.canvas_entrain.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas_entrain, tab3).pack(fill=tk.X)

    # ------------------------------------------------------------------
    # Widget state callbacks
    # ------------------------------------------------------------------
    def _on_atm_mode_change(self):
        is_standard = self.var_atm_mode.get() == "standard"
        state_std = "normal" if is_standard else "disabled"
        state_prof = "disabled" if is_standard else "normal"
        self.entry_Ta0.configure(state=state_std)
        self.lbl_Ta0.configure(foreground="black" if is_standard else "gray")
        self.btn_browse.configure(state=state_prof)
        self.combo_tunits.configure(state=("readonly" if not is_standard else "disabled"))
        self.combo_punits.configure(state=("readonly" if not is_standard else "disabled"))

    def _on_entrainment_mode_change(self):
        is_variable = self.var_ent_mode.get() == "variable"
        state_var = "normal" if is_variable else "disabled"
        state_const = "disabled" if is_variable else "normal"
        self.entry_c.configure(state=state_var)
        self.entry_amax.configure(state=state_var)
        self.entry_aconst.configure(state=state_const)

    def _browse_profile(self):
        path = filedialog.askopenfilename(
            title="Select atmospheric profile file",
            filetypes=[("CSV files", "*.csv"), ("Text files", "*.txt"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            profile = load_atm_profile(
                path,
                temp_units=self.var_tunits.get(),
                press_units=self.var_punits.get(),
            )
        except Exception as exc:
            messagebox.showerror("Failed to load profile", str(exc))
            return
        self.atm_profile = profile
        self.atm_profile_path = path
        display = path if len(path) < 45 else "..." + path[-42:]
        self.lbl_profile_path.configure(text=display, foreground="black")
        self.status_var.set(f"Loaded profile with {len(profile['z'])} levels.")

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------
    def _parse_float(self, var, name):
        try:
            return float(var.get())
        except ValueError:
            raise ModelError(f"'{name}' must be a number (got '{var.get()}').")

    def _on_run(self):
        try:
            z0_asl = self._parse_float(self.var_z0, "Vent altitude")
            r0 = self._parse_float(self.var_r0, "Vent radius")
            u0 = self._parse_float(self.var_u0, "Initial velocity")
            Tj0_C = self._parse_float(self.var_Tj0, "Initial ring temperature")
            zmax_rise = self._parse_float(self.var_zmax, "Max height")

            atm_mode = self.var_atm_mode.get()
            atm_profile = None
            Ta0_target_C = 15.0
            if atm_mode == "standard":
                Ta0_target_C = self._parse_float(self.var_Ta0, "Ambient temperature")
            else:
                if self.atm_profile is None:
                    raise ModelError("Select an atmospheric profile file first "
                                      "(or switch to 'Standard atmosphere').")
                atm_profile = self.atm_profile

            entrainment_mode = self.var_ent_mode.get()
            c_coef = self._parse_float(self.var_c, "c coefficient") if entrainment_mode == "variable" else 1.0
            alpha_max = self._parse_float(self.var_alpha_max, "alpha cap") if entrainment_mode == "variable" else 1.0
            alpha_const = self._parse_float(self.var_alpha_const, "constant alpha") if entrainment_mode == "constant" else 0.1

            self.status_var.set("Running...")
            self.btn_run.configure(state="disabled")
            self.update_idletasks()

            res = run_vortex_ring_model(
                z0_asl=z0_asl, r0=r0, u0=u0, Tj0_C=Tj0_C,
                atm_mode=atm_mode, atm_profile=atm_profile, Ta0_target_C=Ta0_target_C,
                entrainment_mode=entrainment_mode, c_coef=c_coef,
                alpha_max=alpha_max, alpha_const=alpha_const,
                zmax_rise=zmax_rise,
            )
            self.last_result = res
            self._show_report(res)
            self._update_plots(res)
            self.btn_save_csv.configure(state="normal")
            self.btn_save_png.configure(state="normal")
            self.status_var.set("Done.")
        except ModelError as exc:
            self.status_var.set("Input error.")
            messagebox.showerror("Invalid input", str(exc))
        except Exception as exc:
            self.status_var.set("Error.")
            messagebox.showerror("Model run failed", f"{exc}\n\n{traceback.format_exc()}")
        finally:
            self.btn_run.configure(state="normal")

    def _show_report(self, res):
        self.txt_diag.delete("1.0", tk.END)
        self.txt_diag.insert(tk.END, format_report(res))

    def _update_plots(self, res):
        z_out, V, M, B = res["z_out"], res["V"], res["M"], res["B"]
        gprime, u, r, alpha_k_z = res["gprime"], res["u"], res["r"], res["alpha_k_z"]
        Ri_z = res["Ri_z"]
        Tp_z, Ta_z = res["Tp_z"], res["Ta_z"]

        # Tab 1: fluxes (Volume, Momentum, Buoyancy)
        for axi in self.axs_fluxes:
            axi.clear()
        ax = self.axs_fluxes
        ax[0].plot(V, z_out, "."); ax[0].set_title("Volume"); ax[0].set_xlabel("V [m^3]")
        ax[0].set_ylabel("Height above vent [m]")
        ax[1].plot(M, z_out, "."); ax[1].set_title("Momentum"); ax[1].set_xlabel("M [m^4/s]")
        ax[2].plot(B, z_out, "."); ax[2].set_title("Buoyancy"); ax[2].set_xlabel("B [m^4/s^2]")
        self.fig_fluxes.tight_layout()
        self.canvas_fluxes.draw_idle()

        # Tab 2: profiles (velocity, radius, reduced gravity, temperature)
        for axi in self.axs_profiles:
            axi.clear()
        ax2 = self.axs_profiles
        ax2[0].plot(u, z_out, "."); ax2[0].set_xlabel("Vertical velocity [m/s]")
        ax2[0].set_ylabel("Height above vent [m]")
        ax2[1].plot(r, z_out, "."); ax2[1].set_xlabel("Radius [m]")
        ax2[2].plot(gprime, z_out, "."); ax2[2].set_xlabel("g' [m/s^2]")
        ax2[3].plot(Tp_z - 273.15, z_out, ".", label="ring")
        ax2[3].plot(Ta_z - 273.15, z_out, "--", label="ambient")
        ax2[3].set_xlabel("Temperature [deg C]"); ax2[3].legend()
        self.fig_profiles.tight_layout()
        self.canvas_profiles.draw_idle()

        # Tab 3: entrainment coefficient and Richardson number. Ri = |g'| r
        # / u^2 diverges as u -> 0 near the rise apex, so (as in the
        # original tank-experiment model) it's plotted on a log x-axis.
        for axi in self.axs_entrain:
            axi.clear()
        ax3 = self.axs_entrain
        ax3[0].plot(alpha_k_z, z_out, ".")
        ax3[0].set_xlabel("alpha_e"); ax3[0].set_ylabel("Height above vent [m]")
        ax3[0].set_xlim(0, 1)
        ax3[1].semilogx(Ri_z, z_out, ".")
        ax3[1].set_xlabel("Richardson number, Ri")
        self.fig_entrain.tight_layout()
        self.canvas_entrain.draw_idle()

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------
    def _on_save_csv(self):
        if self.last_result is None:
            return
        path = filedialog.asksaveasfilename(
            title="Save results as CSV", defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
        )
        if not path:
            return
        res = self.last_result
        try:
            import csv
            with open(path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["z_above_vent_m", "V_m3", "M_m4s", "B_m4s2",
                                  "alpha_e", "Ri", "u_ms", "r_m", "gprime_ms2",
                                  "Tring_K", "Tambient_K", "Pambient_Pa"])
                for i in range(len(res["z_out"])):
                    writer.writerow([
                        res["z_out"][i], res["V"][i], res["M"][i], res["B"][i],
                        res["alpha_k_z"][i], res["Ri_z"][i], res["u"][i], res["r"][i],
                        res["gprime"][i], res["Tp_z"][i], res["Ta_z"][i],
                        res["Pa_z"][i],
                    ])
            self.status_var.set(f"Saved results to {path}")
        except Exception as exc:
            messagebox.showerror("Save failed", str(exc))

    def _on_save_png(self):
        if self.last_result is None:
            return
        path = filedialog.asksaveasfilename(
            title="Save plots as PNG (base name, all tabs will be saved)",
            defaultextension=".png", filetypes=[("PNG files", "*.png")],
        )
        if not path:
            return
        base = path[:-4] if path.lower().endswith(".png") else path
        try:
            self.fig_fluxes.savefig(base + "_fluxes.png", dpi=150)
            self.fig_profiles.savefig(base + "_profiles.png", dpi=150)
            self.fig_entrain.savefig(base + "_entrainment_Ri.png", dpi=150)
            self.status_var.set(
                f"Saved {base}_fluxes.png, {base}_profiles.png, "
                f"{base}_entrainment_Ri.png"
            )
        except Exception as exc:
            messagebox.showerror("Save failed", str(exc))


if __name__ == "__main__":
    try:
        app = VortexRingGUI()
    except tk.TclError as exc:
        import sys
        print("Could not start the GUI: Tkinter has no display to draw on.")
        print(f"  ({exc})")
        print()
        print("This is an environment/Tk-configuration issue, not a problem with")
        print("this script. It usually happens when the Python you're running was")
        print("built against an X11 Tk instead of native macOS Quartz (common with")
        print("some MacPorts Python/Tk installs) -- it tries to open an X server")
        print("that isn't running. To fix it, try one of:")
        print("  1. Run the script with a different Python that has native Tk")
        print("     support -- e.g. a conda environment's own interpreter")
        print("     ('python vortex_ring_gui.py' instead of an explicit")
        print("     /opt/local/bin/python3.x path), or the python.org macOS installer.")
        print("  2. Install and launch XQuartz (https://www.xquartz.org), then retry.")
        print("  3. If using MacPorts, reinstall Tk with the Quartz variant:")
        print("       sudo port install tk +quartz")
        sys.exit(1)
    app.mainloop()
