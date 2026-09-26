# Groundwater modelling with MODFLOW 6 and Python

Groundwater flow and solute-transport models built with **MODFLOW 6** (USGS) and **FloPy**.
Each project is a complete, runnable Python script with its own README that explains the problem,
the model design, the results and what I learned.

**Author:** Qazi Bilal ([@Qazi-Bilal](https://github.com/Qazi-Bilal))

## Projects

| # | Project | What it shows | Key MODFLOW 6 features |
|---|---|---|---|
| 01 | [Pumping test vs. Theis](01_theis_pumping_test/) | Verifying a transient model against an analytical solution, then recovering *T* and *S* by curve fitting | Transient `STO`, `WEL`, telescoping grid |
| 02 | [River–aquifer interaction](02_river_aquifer_interaction/) | Seasonal gaining river; how much of a well's pumping is taken from the river | `RIV`, `RCH`, unconfined + Newton, scenario comparison |
| 03 | [Seawater intrusion (Henry)](03_henry_saltwater_intrusion/) | A saltwater wedge in a coastal aquifer, and the effect of reduced fresh-water inflow | Coupled `GWF`–`GWT`, `BUY` density package |
| – | Tank experiment vs. Domenico *(learning exercise, root folder)* | 2D transport in a lab tank compared with the Domenico analytical solution; velocity and grain-size experiments | `GWT` with `DSP`, `SSM`, `CHD` concentrations |

<p>
  <img src="01_theis_pumping_test/figures/drawdown_vs_time.png" width="32%">
  <img src="02_river_aquifer_interaction/figures/water_table_map.png" width="32%">
  <img src="03_henry_saltwater_intrusion/figures/saltwater_wedge.png" width="32%">
</p>

## Getting started

1. Install Python 3.10+ and the packages:
   ```bash
   pip install -r requirements.txt
   ```
2. Get the MODFLOW 6 executable, either with `get-modflow :flopy` or from the
   [USGS executables page](https://github.com/MODFLOW-ORG/executables).
3. Run any project:
   ```bash
   cd 01_theis_pumping_test
   python theis_pumping_test.py
   ```
   Each script looks for `mf6` on the PATH. If it isn't found there, edit `MF6_PATH` at the top of the script.

The model input and output files are written to a `model/` folder inside each project. They are not
stored in this repository, because each script recreates them.

## Skills demonstrated

- Building structured-grid MODFLOW 6 models entirely in Python with FloPy
- Flow (GWF), transport (GWT), and density-coupled flow and transport
- Verification against analytical solutions (Theis, Domenico) and published benchmarks (Henry)
- Scenario analysis for water-management questions
- Post-processing heads, budgets and concentrations; publication-style figures and animations

## References

- Langevin, C.D. et al. (2017). *Documentation for the MODFLOW 6 Groundwater Flow Model.* USGS TM 6-A55.
- Langevin, C.D. et al. (2022). *Documentation for the MODFLOW 6 Groundwater Transport Model.* USGS TM 6-A61.
- Bakker, M. et al. (2016). Scripting MODFLOW model development using Python and FloPy. *Groundwater*, 54(5).
- Theis, C.V. (1935). The relation between the lowering of the piezometric surface and the rate and
  duration of discharge of a well using ground-water storage. *Trans. AGU*, 16.
- Henry, H.R. (1964). Effects of dispersion on salt encroachment in coastal aquifers. USGS WSP 1613-C.
