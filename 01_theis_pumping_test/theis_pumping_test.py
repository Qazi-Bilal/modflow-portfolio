"""
Project 01 - Pumping test: MODFLOW 6 vs. the Theis solution
============================================================

A single well pumps a confined aquifer at a constant rate for one day.
We simulate the growing "cone of depression" with MODFLOW 6 and check it
against the classic Theis (1935) analytical solution. Then we do what a
hydrogeologist does with real field data: fit Theis to the simulated
drawdown and recover the aquifer's transmissivity (T) and storativity (S).

Packages used: DIS, IC, NPF, STO, CHD, WEL, OC (GWF only)

Run:  python theis_pumping_test.py
"""
import shutil
from pathlib import Path

import flopy
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import exp1


# ----------------------------------------------------------------------
# Where is mf6? Change MF6_PATH if the script cannot find it.
# ----------------------------------------------------------------------
MF6_PATH = r"D:\python\mf6.exe"


def find_mf6():
    for candidate in ("mf6", MF6_PATH):
        if shutil.which(candidate) or Path(candidate).is_file():
            return candidate
    raise FileNotFoundError("mf6 not found - set MF6_PATH at the top of the script")


HERE = Path(__file__).resolve().parent
WS = HERE / "model"            # MODFLOW files go here (not uploaded to GitHub)
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

# ----------------------------------------------------------------------
# 1. Aquifer and test parameters (units: metres, days)
# ----------------------------------------------------------------------
K = 10.0          # hydraulic conductivity (m/d)
b = 10.0          # aquifer thickness (m)
T = K * b         # transmissivity (m2/d)  = 100
S = 1.0e-4        # storativity (-)
Ss = S / b        # specific storage (1/m) - what MODFLOW needs
Q = 500.0         # pumping rate (m3/d)
h0 = 50.0         # initial head (m)
t_end = 1.0       # test duration (d)
obs_r = [10.0, 50.0, 200.0]   # observation wells: distance from pumping well (m)

# ----------------------------------------------------------------------
# 2. Grid: fine near the well, coarse far away (a "telescoping" grid).
#    The far boundary must be so far away that the cone never reaches it
#    in one day - otherwise the model is no longer "infinite" like Theis.
# ----------------------------------------------------------------------
dx0, growth, half_width = 2.0, 1.15, 20_000.0
half = [dx0]
while sum(half) < half_width:
    half.append(half[-1] * growth)
widths = np.array(half[::-1] + [dx0] + half)      # symmetric, well in the middle
n = len(widths)
ic = n // 2                                       # index of the well cell
centers = np.cumsum(widths) - widths / 2
centers -= centers[ic]                            # x = 0 at the well
print(f"Grid: {n} x {n} cells, cell size {widths.min():.0f} m to {widths.max():.0f} m")

# ----------------------------------------------------------------------
# 3. Build the MODFLOW 6 model
# ----------------------------------------------------------------------
sim = flopy.mf6.MFSimulation(sim_name="theis", sim_ws=str(WS), exe_name=find_mf6())
flopy.mf6.ModflowTdis(sim, time_units="days", nper=1,
                      perioddata=[(t_end, 60, 1.12)])   # 60 steps, growing in size
flopy.mf6.ModflowIms(sim, complexity="SIMPLE", outer_dvclose=1e-6, inner_dvclose=1e-8)

gwf = flopy.mf6.ModflowGwf(sim, modelname="theis")
flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=n, ncol=n, delr=widths, delc=widths,
                        top=h0 + 10, botm=h0 + 10 - b)
flopy.mf6.ModflowGwfic(gwf, strt=h0)
flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=K)                  # confined
flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=Ss, transient={0: True})

# Far boundary: fixed head (never reached by the cone in 1 day)
edge = [(0, i, j) for i in range(n) for j in range(n) if i in (0, n - 1) or j in (0, n - 1)]
flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[c, h0] for c in edge])
flopy.mf6.ModflowGwfwel(gwf, stress_period_data=[[(0, ic, ic), -Q]])
flopy.mf6.ModflowGwfoc(gwf, head_filerecord="theis.hds", saverecord=[("HEAD", "ALL")])

sim.write_simulation(silent=True)
ok, _ = sim.run_simulation(silent=True)
assert ok, "MODFLOW failed - see model/mfsim.lst"

# ----------------------------------------------------------------------
# 4. Read results: drawdown along the row through the well
# ----------------------------------------------------------------------
hds = gwf.output.head()
times = np.array(hds.get_times())
heads = hds.get_alldata()[:, 0, ic, :]               # (time, column)
right = slice(ic, n)                                 # well -> east edge
dd_model = {r: np.array([np.interp(r, centers[right], h0 - h[right])
                         for h in heads]) for r in obs_r}


def theis(r, t, T, S):
    """Theis drawdown s = Q/(4 pi T) * W(u),  u = r^2 S / (4 T t)."""
    u = r**2 * S / (4 * T * t)
    return Q / (4 * np.pi * T) * exp1(u)


# ----------------------------------------------------------------------
# 5. "Pumping test analysis": fit Theis to the 50 m observation well
# ----------------------------------------------------------------------
r_fit = 50.0
mask = dd_model[r_fit] > 0.01
(T_fit, S_fit), _ = curve_fit(lambda t, T, S: theis(r_fit, t, T, S),
                              times[mask], dd_model[r_fit][mask],
                              p0=[50.0, 1e-3], bounds=([1, 1e-7], [1e4, 1e-1]))
print(f"Fitted  T = {T_fit:6.1f} m2/d   (true {T:.1f})   error {100*(T_fit/T-1):+.1f} %")
print(f"Fitted  S = {S_fit:.2e}      (true {S:.1e})  error {100*(S_fit/S-1):+.1f} %")
for r in obs_r:
    err = np.max(np.abs(dd_model[r] - theis(r, times, T, S)))
    print(f"r = {r:5.0f} m : max |model - Theis| = {err*100:.1f} cm")

# ----------------------------------------------------------------------
# 6. Figures
# ----------------------------------------------------------------------
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
MARKERS = ["o", "s", "^"]
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": "white"})

# (a) Drawdown vs time - the classic pumping-test plot
fig, ax = plt.subplots(figsize=(7, 4.5))
tt = np.logspace(-4, 0, 200)
for r, c, m in zip(obs_r, COLORS, MARKERS):
    ax.plot(tt, theis(r, tt, T, S), color=c, lw=2)
    ax.plot(times[::3], dd_model[r][::3], m, color=c, ms=6, mfc="white", mew=1.5)
    ax.annotate(f"r = {r:.0f} m", xy=(1.0, theis(r, 1.0, T, S)), xytext=(6, 0),
                textcoords="offset points", va="center", color=INK, fontsize=9)
ax.plot([], [], "-", color=INK2, lw=2, label="Theis (analytical)")
ax.plot([], [], "o", color=INK2, mfc="white", label="MODFLOW 6")
ax.set_xscale("log")
ax.set_xlim(1e-4, 1.0)
ax.set_xlabel("Time since pumping started (days)")
ax.set_ylabel("Drawdown (m)")
ax.set_title("Drawdown at three observation wells", loc="left", color=INK)
ax.grid(True, which="major", color=GRID, lw=0.8)
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(FIG / "drawdown_vs_time.png", dpi=150, bbox_inches="tight")

# (b) Map of the cone of depression after 1 day
fig, ax = plt.subplots(figsize=(6, 5))
dd_map = h0 - hds.get_data(totim=times[-1])[0]
X, Y = np.meshgrid(centers, -centers)
zoom = np.abs(centers) <= 320
Z = dd_map[np.ix_(zoom, zoom)]
cf = ax.contourf(X[np.ix_(zoom, zoom)], Y[np.ix_(zoom, zoom)], Z,
                 levels=np.arange(0, 4.01, 0.25), cmap="Blues", extend="max")
cs = ax.contour(X[np.ix_(zoom, zoom)], Y[np.ix_(zoom, zoom)], Z,
                levels=[1, 1.5, 2, 3], colors=INK2, linewidths=0.8)
ax.clabel(cs, fmt="%.1f m", fontsize=8)
ax.plot(0, 0, marker="v", color=INK, ms=9, ls="none", label="Pumping well")
label_pos = {10.0: ((-2, -16), "right"), 50.0: ((2, -16), "left"), 200.0: ((0, 10), "center")}
for r, c, m in zip(obs_r, COLORS, MARKERS):
    ax.plot(r, 0, m, color=c, ms=8, mec="white", mew=1.2)
    off, ha = label_pos[r]
    ax.annotate(f"obs. {r:.0f} m", (r, 0), xytext=off, textcoords="offset points",
                ha=ha, fontsize=8, color=INK,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
fig.colorbar(cf, ax=ax, label="Drawdown (m)", shrink=0.85)
ax.set_aspect("equal")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_title("Cone of depression after 1 day", loc="left", color=INK)
ax.legend(frameon=False, loc="lower left", fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "cone_of_depression.png", dpi=150, bbox_inches="tight")

# (c) Parameter recovery
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(times[mask], dd_model[r_fit][mask], "o", color=COLORS[1], mfc="white", mew=1.5,
        ms=6, label="'Measured' (MODFLOW 6, r = 50 m)")
ax.plot(tt, theis(r_fit, tt, T_fit, S_fit), color=INK, lw=1.5,
        label=f"Theis fit: T = {T_fit:.0f} m²/d, S = {S_fit:.1e}")
ax.set_xscale("log")
ax.set_xlabel("Time since pumping started (days)")
ax.set_ylabel("Drawdown (m)")
ax.set_title("Recovering T and S from the drawdown curve", loc="left", color=INK)
ax.text(0.98, 0.05, f"True values: T = {T:.0f} m²/d, S = {S:.0e}", transform=ax.transAxes,
        ha="right", color=INK2, fontsize=9)
ax.grid(True, color=GRID, lw=0.8)
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(FIG / "parameter_fit.png", dpi=150, bbox_inches="tight")

print(f"Figures saved in {FIG}")
plt.show()
