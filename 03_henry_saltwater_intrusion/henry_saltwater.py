"""
Project 03 - Seawater intrusion: the Henry problem with MODFLOW 6
==================================================================

A vertical cross-section (2 m long, 1 m deep) of a coastal aquifer.
Fresh groundwater enters from the land side (left) and flows to the sea
(right). Seawater is heavier (1025 kg/m3 vs 1000 kg/m3), so it slides in
underneath the fresh water and forms a *saltwater wedge*.

This needs two-way coupling: the salt concentration (GWT) changes the
water density, and the density changes the flow (GWF). MODFLOW 6 does this
with the BUY (buoyancy) package.

We run two scenarios:
  1. Henry's classic case (freshwater inflow 5.7 m3/d per m of coast)
  2. Half the freshwater inflow - e.g. less recharge or pumping inland
and compare how far the seawater intrudes.

Packages used: GWF (DIS, IC, NPF, BUY, CHD, WEL, OC) + GWT (DIS, IC, ADV, DSP,
MST, SSM, OC) + GWF-GWT exchange

Run:  python henry_saltwater.py    (takes about 1-3 minutes)
"""
import shutil
from pathlib import Path

import flopy
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LinearSegmentedColormap

MF6_PATH = r"D:\python\mf6.exe"   # change if mf6 is elsewhere


def find_mf6():
    for candidate in ("mf6", MF6_PATH):
        if shutil.which(candidate) or Path(candidate).is_file():
            return candidate
    raise FileNotFoundError("mf6 not found - set MF6_PATH at the top of the script")


HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

# ----------------------------------------------------------------------
# 1. Henry (1964) parameters, as in the MODFLOW 6 example "ex-gwt-henry"
# ----------------------------------------------------------------------
Lx, Lz = 2.0, 1.0                 # length and depth (m)
nlay, ncol = 40, 80               # cells: 2.5 cm x 2.5 cm
delr, delz = Lx / ncol, Lz / nlay
K = 864.0                         # hydraulic conductivity (m/d)
porosity = 0.35
diffc = 0.57024                   # diffusion coefficient (m2/d)
c_sea, rho_fresh, rho_sea = 35.0, 1000.0, 1024.5
drhodc = (rho_sea - rho_fresh) / c_sea          # 0.7 kg/m3 per g/L
t_end, nstp = 0.5, 500            # 0.5 days, 500 time steps
save_every = 10                   # save concentration every 10 steps (for the animation)


def build_and_run(name, q_fresh):
    ws = HERE / "model" / name
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=str(ws), exe_name=find_mf6())
    flopy.mf6.ModflowTdis(sim, time_units="days", nper=1, perioddata=[(t_end, nstp, 1.0)])

    # ---------------- flow model ----------------
    gwf = flopy.mf6.ModflowGwf(sim, modelname="flow", save_flows=True)
    ims_f = flopy.mf6.ModflowIms(sim, complexity="MODERATE", linear_acceleration="BICGSTAB",
                                 outer_dvclose=1e-8, inner_dvclose=1e-9, filename="flow.ims")
    sim.register_ims_package(ims_f, [gwf.name])
    botm = np.linspace(Lz - delz, 0.0, nlay)
    flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=1, ncol=ncol, delr=delr, delc=1.0,
                            top=Lz, botm=botm)
    flopy.mf6.ModflowGwfic(gwf, strt=Lz)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=K, save_specific_discharge=True)

    # BUY: density = 1000 + 0.7 * concentration, taken from the transport model "trans"
    flopy.mf6.ModflowGwfbuy(gwf, denseref=rho_fresh,
                            packagedata=[(0, drhodc, 0.0, "trans", "CONCENTRATION")])

    # Land side: fresh water enters (WEL, spread over all layers, concentration 0)
    wel = [[(k, 0, 0), q_fresh / nlay, 0.0] for k in range(nlay)]
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data=wel, auxiliary="CONCENTRATION", pname="WEL-1")
    # Sea side: sea level, with seawater concentration and density
    chd = [[(k, 0, ncol - 1), Lz, c_sea, rho_sea] for k in range(nlay)]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd,
                            auxiliary=["CONCENTRATION", "DENSITY"], pname="CHD-1")
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="flow.cbc", head_filerecord="flow.hds",
                           saverecord=[("BUDGET", "LAST"), ("HEAD", "LAST")])

    # ---------------- transport model ----------------
    gwt = flopy.mf6.ModflowGwt(sim, modelname="trans")
    ims_t = flopy.mf6.ModflowIms(sim, complexity="MODERATE", linear_acceleration="BICGSTAB",
                                 outer_dvclose=1e-6, inner_dvclose=1e-7, filename="trans.ims")
    sim.register_ims_package(ims_t, [gwt.name])
    flopy.mf6.ModflowGwtdis(gwt, nlay=nlay, nrow=1, ncol=ncol, delr=delr, delc=1.0,
                            top=Lz, botm=botm)
    flopy.mf6.ModflowGwtic(gwt, strt=c_sea)            # start: aquifer full of seawater
    flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    flopy.mf6.ModflowGwtdsp(gwt, diffc=diffc, alh=0.0, ath1=0.0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
    flopy.mf6.ModflowGwtssm(gwt, sources=[("WEL-1", "AUX", "CONCENTRATION"),
                                          ("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtoc(gwt, concentration_filerecord="trans.ucn",
                           saverecord=[("CONCENTRATION", "FREQUENCY", save_every)])

    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="flow", exgmnameb="trans")
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    assert ok, f"MODFLOW failed - see {ws / 'mfsim.lst'}"
    return gwf, gwt


scenarios = {"Henry case (Q = 5.70 m³/d)": 5.7024,
             "Half the fresh inflow (Q = 2.85 m³/d)": 5.7024 / 2}
results = {}
for label, q in scenarios.items():
    print(f"Running: {label} ...")
    name = "henry" if q > 5 else "henry_half"
    gwf, gwt = build_and_run(name, q)
    ucn = gwt.output.concentration()
    results[label] = dict(gwf=gwf, times=np.array(ucn.get_times()),
                          conc=ucn.get_alldata()[:, :, 0, :])      # (time, layer, column)

x = (np.arange(ncol) + 0.5) * delr
z = Lz - (np.arange(nlay) + 0.5) * delz


def toe_position(c):
    """How far inland (m, from the sea boundary) the 50 % seawater line reaches at the bottom."""
    bottom = c[-1, :] / c_sea
    inland = np.where(bottom >= 0.5)[0]
    return Lx - x[inland.min()] if inland.size else 0.0


for label, r in results.items():
    print(f"{label:40s}: 50 % isochlor toe {toe_position(r['conc'][-1]):.2f} m inland")

# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
INK, INK2 = "#0b0b0b", "#52514e"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]
cmap = LinearSegmentedColormap.from_list("salt", ["#fcfcfb"] + SEQ)
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": "white"})
levels = np.linspace(0, c_sea, 15)


def draw(ax, c, gwf=None, title=""):
    cf = ax.contourf(x, z, c, levels=levels, cmap=cmap)
    cs = ax.contour(x, z, c / c_sea, levels=[0.1, 0.5, 0.9], colors=INK, linewidths=[0.8, 1.6, 0.8])
    ax.clabel(cs, fmt=lambda v: f"{v:.0%}", fontsize=8)
    if gwf is not None:
        spd = gwf.output.budget().get_data(text="DATA-SPDIS")[-1]
        qx, _, qz = flopy.utils.postprocessing.get_specific_discharge(spd, gwf)
        s = 4
        ax.quiver(x[2::s], z[2::s], qx[2::s, 0, 2::s], qz[2::s, 0, 2::s],
                  color=INK2, scale=90, width=0.002)
    ax.set_aspect("equal")
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Lz)
    ax.set_title(title, loc="left", color=INK, fontsize=10)
    box = dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85)
    ax.text(0.015, 0.04, "LAND side\nfresh water in", transform=ax.transAxes, fontsize=8,
            color=INK, va="bottom", bbox=box)
    ax.text(0.985, 0.04, "SEA side", transform=ax.transAxes, fontsize=8, color="white",
            ha="right", va="bottom")
    return cf


# (a) Final state, both scenarios
fig, axes = plt.subplots(2, 1, figsize=(8, 7.2), sharex=True)
for ax, (label, r) in zip(axes, results.items()):
    cf = draw(ax, r["conc"][-1], r["gwf"],
              f"{label}  -  toe {toe_position(r['conc'][-1]):.2f} m inland")
    ax.set_ylabel("Elevation (m)")
axes[-1].set_xlabel("Distance from land boundary (m)")
fig.colorbar(cf, ax=axes, label="Salt concentration (g/L)", shrink=0.8)
fig.savefig(FIG / "saltwater_wedge.png", dpi=150, bbox_inches="tight")

# (b) Position of the 50 % line along the bottom over time
fig, ax = plt.subplots(figsize=(7, 3.8))
for (label, r), color, mk in zip(results.items(), ["#2a78d6", "#eb6834"], ["o", "s"]):
    toe = [toe_position(c) for c in r["conc"]]
    ax.plot(r["times"], toe, "-", marker=mk, markevery=5, color=color, lw=2, ms=5, label=label)
ax.set_xlabel("Time (days)")
ax.set_ylabel("50 % toe inland from sea (m)")
ax.set_title("How fast the saltwater wedge settles", loc="left", color=INK)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(True, axis="y", color="#e4e3df", lw=0.8)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(FIG / "toe_position.png", dpi=150, bbox_inches="tight")

# (c) Animation of the classic case
r = results[list(results)[0]]
fig, ax = plt.subplots(figsize=(8, 4.2))


def frame(i):
    ax.clear()
    draw(ax, r["conc"][i], title=f"Freshwater pushing seawater out  -  t = {r['times'][i]:.2f} d")
    ax.set_xlabel("Distance from land boundary (m)")
    ax.set_ylabel("Elevation (m)")


anim = FuncAnimation(fig, frame, frames=len(r["times"]), interval=120)
anim.save(FIG / "henry_animation.gif", writer=PillowWriter(fps=8), dpi=90)
print(f"Figures saved in {FIG}")
plt.close(fig)
plt.show()

