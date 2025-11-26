import tkinter as tk
from tkinter import ttk
import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

# ------------------- RCA kinetic model -------------------
def rca_model(P, t, params):
    E, T, dNTP, Pr, Mg = (
        params["E"],
        params["T"],
        params["dNTP"],
        params["Pr"],
        params["Mg"]
    )

    vmax = 50.0                     # bp/s per polymerase
    Km_template = 1e-9
    Km_primer = 1e-9
    Km_Mg = 1e-3

    primer_factor = Pr / (Km_primer + Pr)
    Mg_factor = Mg / (Km_Mg + Mg)
    template_factor = T / (Km_template + T)

    rate = vmax * E * primer_factor * Mg_factor * template_factor
    return rate

def simulate(params, tmax=3600):
    t = np.linspace(0, tmax, 400)
    P0 = 0
    P = odeint(rca_model, P0, t, args=(params,))
    return t, P.flatten()

def final_yield(param_name, base_params, sweep_range):
    yields = []
    for val in sweep_range:
        params = base_params.copy()
        params[param_name] = val
        _, P = simulate(params)
        yields.append(P[-1])
    return yields

# ------------------- GUI -------------------
root = tk.Tk()
root.title("RCA Digital Twin – Live Dashboard")

# Main frame with 2 columns
main = tk.Frame(root)
main.pack(fill="both", expand=True)

# Left: controls
controls = tk.Frame(main)
controls.grid(row=0, column=0, sticky="ns")

# Right: plots
plot_frame = tk.Frame(main)
plot_frame.grid(row=0, column=1, sticky="nsew")

# Expand plot area
main.columnconfigure(1, weight=1)
main.rowconfigure(0, weight=1)

# ------------------- Sliders -------------------
def make_slider(label, from_, to, resolution, row):
    tk.Label(controls, text=label).grid(row=row, column=0)
    s = tk.Scale(controls, from_=from_, to=to, resolution=resolution,
                 orient="horizontal", length=250, command=lambda x:update_plots())
    s.grid(row=row, column=1)
    return s

poly = make_slider("Polymerase (nM)", 0.01, 100, 0.01, 0)
templ = make_slider("Template (nM)", 0.001, 10, 0.001, 1)
primer = make_slider("Primer (nM)", 0.001, 10, 0.001, 2)
dntp = make_slider("dNTP (mM)", 0.01, 10, 0.01, 3)
mg = make_slider("Mg2+ (mM)", 0.01, 10, 0.01, 4)

time_slider = make_slider("Time (min)", 5, 240, 1, 5)

# Default values
poly.set(1)
templ.set(1)
primer.set(1)
dntp.set(1)
mg.set(2)
time_slider.set(60)

# ------------------- Matplotlib Figure -------------------
fig, axs = plt.subplots(2, 3, figsize=(12, 6))
axs = axs.flatten()

canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.get_tk_widget().pack(fill="both", expand=True)

# ------------------- Update function -------------------
def update_plots():
    params = {
        "E": poly.get()*1e-9,
        "T": templ.get()*1e-9,
        "Pr": primer.get()*1e-9,
        "dNTP": dntp.get()*1e-3,
        "Mg": mg.get()*1e-3
    }

    t, P = simulate(params, tmax=time_slider.get()*60)

    sweep = np.linspace(0.1, 10, 50)

    axs[0].clear()
    axs[0].plot(t/60, P)
    axs[0].set_title("Time Series Product")
    axs[0].set_xlabel("Time (min)")
    axs[0].set_ylabel("Product")

    axs[1].clear()
    axs[1].plot(sweep, final_yield("T", params, sweep*1e-9))
    axs[1].set_title("Final Yield vs Template")

    axs[2].clear()
    axs[2].plot(sweep, final_yield("E", params, sweep*1e-9))
    axs[2].set_title("Final Yield vs Polymerase")

    axs[3].clear()
    axs[3].plot(sweep, final_yield("Pr", params, sweep*1e-9))
    axs[3].set_title("Final Yield vs Primer")

    axs[4].clear()
    axs[4].plot(sweep, final_yield("dNTP", params, sweep*1e-3))
    axs[4].set_title("Final Yield vs dNTP")

    axs[5].clear()
    axs[5].plot(sweep, final_yield("Mg", params, sweep*1e-3))
    axs[5].set_title("Final Yield vs Mg2+")

    fig.tight_layout()
    canvas.draw()

update_plots()
root.mainloop()
