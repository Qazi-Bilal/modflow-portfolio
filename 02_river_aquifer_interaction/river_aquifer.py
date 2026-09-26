"""
Project 02 - River-aquifer interaction and streamflow depletion
================================================================

A river crosses an unconfined aquifer. Groundwater flows from the uplands
(west) toward the river, which is fed by it (a "gaining" river). River stage
and rain recharge change with the seasons. In year 2 a water-supply well is
switched on 500 m from the river.

Question a water manager would ask: *how much of the water pumped by the
well is actually taken from the river?*  (= streamflow depletion)

We answer it by running two scenarios - with and without the well - and
comparing the river-aquifer exchange (RIV package budget).

Packages used: DIS, IC, NPF, STO, CHD, RIV, RCH, WEL, OC (GWF only)

Run:  python river_aquifer.py
"""
import shutil
from pathlib import Path

import flopy
import matplotlib.pyplot as plt
import numpy as np

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
# 1. Parameters (metres, days)
# ----------------------------------------------------------------------
nrow, ncol, cell = 40, 60, 100.0          # 4 km x 6 km, 100 m cells
top, bot = 50.0, 0.0
K, Sy, Ss = 15.0, 0.15, 1e-5              # unconfined sandy gravel
h_west = 42.0                             # upland boundary head
riv_col = 45                              # river runs north-south here (x = 4.55 km)
riv_cond = 1.0 * cell * 20.0 / 1.0        # K_bed * length * width / bed thickness  (m2/d)
riv_bot = 26.0
well_cell = (0, 20, 40)                   # 500 m west of the river
Q_well = 4000.0                           # m3/d  (~ a town of 20 000 people)
n_years, month = 3, 365.25 / 12
well_on_month = 12                        # well starts at the beginning of year 2

months = np.arange(n_years * 12)
# seasonal signals (Northern-hemisphere-like): wet winter / spring, dry summer
stage = 30.0 + 1.0 * np.cos(2 * np.pi * (months - 3) / 12)            # peak in April
recharge = 0.0006 * (1 + np.cos(2 * np.pi * (months - 1) / 12))      # peak in Feb, ~0 in Aug


def build_and_run(name, with_well):
    ws = HERE / "model" / name
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=str(ws), exe_name=find_mf6())

    # Stress period 0 = steady state (starting condition), then monthly periods
    perioddata = [(1.0, 1, 1.0)] + [(month, 6, 1.0)] * len(months)
    flopy.mf6.ModflowTdis(sim, time_units="days", nper=len(perioddata), perioddata=perioddata)
    flopy.mf6.ModflowIms(sim, complexity="MODERATE", outer_dvclose=1e-4)

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True,
                               newtonoptions="UNDER_RELAXATION")
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=nrow, ncol=ncol, delr=cell, delc=cell,
                            top=top, botm=bot)
    flopy.mf6.ModflowGwfic(gwf, strt=36.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=K, save_specific_discharge=True)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, sy=Sy, ss=Ss,
                            steady_state={0: True}, transient={1: True})

    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, i, 0), h_west] for i in range(nrow)])

    # RIV: [cell, stage, conductance, river bottom]  - one entry per month
    riv = {0: [[(0, i, riv_col), stage.mean(), riv_cond, riv_bot] for i in range(nrow)]}
    rch = {0: recharge.mean()}
    for m in months:
        riv[m + 1] = [[(0, i, riv_col), stage[m], riv_cond, riv_bot] for i in range(nrow)]
        rch[m + 1] = recharge[m]
    flopy.mf6.ModflowGwfriv(gwf, stress_period_data=riv, pname="RIV")
    flopy.mf6.ModflowGwfrcha(gwf, recharge=rch)

    if with_well:
        wel = {0: [], well_on_month + 1: [[well_cell, -Q_well]]}
        flopy.mf6.ModflowGwfwel(gwf, stress_period_data=wel)

    flopy.mf6.ModflowGwfoc(gwf, head_filerecord=f"{name}.hds", budget_filerecord=f"{name}.cbc",
                           saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")])
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    assert ok, f"MODFLOW failed - see {ws / 'mfsim.lst'}"
    return gwf


print("Running scenario without well ...")
gwf_base = build_and_run("no_well", with_well=False)
print("Running scenario with well ...")
gwf_well = build_and_run("with_well", with_well=True)


def river_exchange(gwf):
    """Net river -> aquifer flow per stress period (m3/d). Negative = river gains."""
    cbc = gwf.output.budget()
    return np.array([rec["q"].sum() for rec in cbc.get_data(text="RIV")])[1:]


riv_base = river_exchange(gwf_base)
riv_well = river_exchange(gwf_well)
depletion = (riv_well - riv_base) / Q_well           # fraction of pumping taken from river
t_months = months + 1

yr3 = months >= 24
print(f"River gain without well : {-riv_base.mean():7.0f} m3/d (average)")
print(f"River gain with well    : {-riv_well.mean():7.0f} m3/d (average)")
print(f"Streamflow depletion in year 3: {100*depletion[yr3].mean():.0f} % of the pumping rate")

# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": "white"})

# (a) River gain over time, two scenarios
fig, ax = plt.subplots(figsize=(8, 4.2))
ax.axvspan(well_on_month + 0.5, t_months[-1] + 0.5, color="#f0efec", zorder=0)
ax.text(well_on_month + 1, 0.97, "well pumping", transform=ax.get_xaxis_transform(),
        color=INK2, fontsize=9, va="top")
ax.plot(t_months, -riv_base, "-o", color=BLUE, lw=2, ms=4, label="No well")
ax.plot(t_months, -riv_well, "-s", color=ORANGE, lw=2, ms=4, label="With well")
ax.annotate("No well", (t_months[-1], -riv_base[-1]), xytext=(6, 0), textcoords="offset points",
            va="center", color=INK, fontsize=9)
ax.annotate("With well", (t_months[-1], -riv_well[-1]), xytext=(6, 0), textcoords="offset points",
            va="center", color=INK, fontsize=9)
ax.axhline(0, color=INK2, lw=0.8)
ax.set_xlim(0.5, t_months[-1] + 0.5)
ax.set_xticks(np.arange(1, t_months[-1] + 1, 3))
ax.set_xlabel("Month")
ax.set_ylabel("Groundwater flow into river (m³/d)")
ax.set_title("River gain from the aquifer", loc="left", color=INK)
ax.grid(True, axis="y", color=GRID, lw=0.8)
ax.legend(frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(FIG / "river_gain.png", dpi=150, bbox_inches="tight")

# (b) Streamflow depletion fraction
fig, ax = plt.subplots(figsize=(8, 3.6))
on = t_months > well_on_month
ax.plot(t_months[on] - well_on_month, 100 * depletion[on], "-o", color=BLUE, lw=2, ms=4)
ax.set_ylim(0, 100)
ax.set_xlabel("Months since the well was switched on")
ax.set_ylabel("Share of pumping\ntaken from the river (%)")
ax.set_title("Streamflow depletion", loc="left", color=INK)
ax.grid(True, axis="y", color=GRID, lw=0.8)
fig.tight_layout()
fig.savefig(FIG / "streamflow_depletion.png", dpi=150, bbox_inches="tight")

# (c) Map: water table and flow direction at the end of summer, with well
kper_aug = 24 + 8                     # August of year 3 (stress period index)
head = gwf_well.output.head().get_data(kstpkper=(5, kper_aug))[0]
spd = gwf_well.output.budget().get_data(text="DATA-SPDIS", kstpkper=(5, kper_aug))[0]
qx, qy, _ = flopy.utils.postprocessing.get_specific_discharge(spd, gwf_well)
x = (np.arange(ncol) + 0.5) * cell / 1000
y = (nrow - np.arange(nrow) - 0.5) * cell / 1000
fig, ax = plt.subplots(figsize=(8, 5))
cf = ax.contourf(x, y, head, levels=np.arange(28, 42.01, 0.25), cmap="Blues_r", extend="both")
cs = ax.contour(x, y, head, levels=np.arange(28, 42.01, 1.0), colors=INK2, linewidths=0.6)
ax.clabel(cs, fmt="%.1f", fontsize=7)
s = 3
mag = np.hypot(qx[0], qy[0]) + 1e-12          # arrows show direction only (same length)
ax.quiver(x[1::s], y[1::s], (qx[0] / mag)[1::s, 1::s], (qy[0] / mag)[1::s, 1::s],
          color=INK2, scale=40, width=0.0025, headwidth=4)
ax.axvline(x[riv_col], color=BLUE, lw=4, alpha=0.8)
ax.text(x[riv_col] + 0.08, y[1], "River", color=INK, fontsize=9, va="top")
ax.plot(x[well_cell[2]], y[well_cell[1]], "v", color=ORANGE, ms=11, mec="white")
ax.annotate("Well", (x[well_cell[2]], y[well_cell[1]]), xytext=(-14, 12),
            textcoords="offset points", ha="right",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85), color=INK, fontsize=9)
ax.text(0.1, y[1], "Fixed head\n(uplands)", color=INK, fontsize=8, va="top")
fig.colorbar(cf, ax=ax, label="Water table (m)")
ax.set_aspect("equal")
ax.set_xlabel("x (km)")
ax.set_ylabel("y (km)")
ax.set_title("Water table and flow direction, August of year 3", loc="left", color=INK)
fig.tight_layout()
fig.savefig(FIG / "water_table_map.png", dpi=150, bbox_inches="tight")

print(f"Figures saved in {FIG}")
plt.show()
