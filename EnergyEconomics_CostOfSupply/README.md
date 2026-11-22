# Replication package — Energy Economics paper

This folder of the repository [SelfGeneration_PV_BESS_BTM_Sizing](https://github.com/pmdeoliveiradejesus/SelfGeneration_PV_BESS_BTM_Sizing) contains the data and code of the paper below. The parent repository holds the sizing models of De Oliveira-De Jesus et al. (2026), including the MILP version solved with HiGHS; `model/sizing_model.py` here is that MILP version with all the cases analyzed in the paper selectable from the command line.

**Paper:** P. M. De Oliveira-De Jesus, *The cost of supplying a prosumer: levelized costs and contract price floors for optimally sized and operated behind-the-meter solar and storage*, submitted to *Energy Economics*.

**Contact:** Paulo M. De Oliveira-De Jesus, Department of Electrical and Electronic Engineering, Universidad de los Andes, Bogotá, Colombia.

## One command

```bash
pip install numpy pandas matplotlib
python reproduce.py
```

This takes a few seconds and needs no optimizer or solver licence. It:

- reads the model inputs in `data/inputs/` and the optimization results in `data/model_outputs/`;
- applies the closed-form accounting of Section 3 of the paper (Eqs. 1–10, Table 2);
- writes Tables 5–8 and the data of Figure 2 to `tables/` (CSV and LaTeX);
- writes Figures 1–3 to `figures/` (PDF);
- writes `key_numbers.txt`, which compares 189 numbers reported in the paper (tables, figures, and text) with the recomputed values.

The script exits with code 0 only if all of them match to the precision reported.

## What reproduces what

| Paper item | Output of `reproduce.py` | Source data |
|---|---|---|
| Table 4 (assumptions), $CRF$, $CRF_e$, $E_{b0}$, $CP_0$ | header of `key_numbers.txt` | parameters in `reproduce.py`, section 1; `data/inputs/` |
| Table 5 (designs, levelized costs, Proposition 1) | `tables/table5_designs.*` | `results_summary_optimal.csv`, `_nocontract.csv`, `_ember.csv` |
| Ember island figures (92.5 % served, 112.2 USD/MWh served) | `key_numbers.txt` | `results_summary_ember_island.csv` |
| Zero and negative prices (532 h, 247 h, 757 MWh, 623 MWh, 77 MWh) | `key_numbers.txt` | `data/inputs/lambda_spain_localtime.inc`, `hourly/solution_optimal.csv` |
| Table 6 (contract price bounds, Proposition 2), Eqs. (7) and (11) | `tables/table6_bounds.*` | as Table 5 |
| Table 7 (30 % investment subsidy) | `tables/table7_subsidy.*` | `results_summary_optimal.csv`, `_sub30.csv` |
| Figure 2 (subsidy sweep, 0–50 %) | `figures/fig_sweep.pdf`, `tables/figure2_subsidy_sweep.*` | `results_summary_sub05.csv` … `_sub50.csv` |
| Table 8, Proposition 3 thresholds (16.0 USD/MWh, 156 and 43 USD/t) | `tables/table8_carbon.*` | as Table 7, plus emission factors in `reproduce.py`, section 1 |
| Figure 1 (decomposition of $LCOE_s$) | `figures/fig_decomposition.pdf` | as Tables 5 and 7 |
| Figure 3 (social cost vs carbon price) | `figures/fig_carbon.pdf` | as Table 8 |

## Data

### `data/inputs/` — hourly inputs of the sizing model (8 760 h, local time, 2024 calendar)

| File | Content | Source |
|---|---|---|
| `lambda_spain_localtime.inc` | Day-ahead spot price, Spain, 2024 [EUR/kWh] | OMIE, https://www.omie.es/ |
| `psi.inc` | Network access charges and regulated charges paid on withdrawals [EUR/kWh] | Spanish 6.1TD tariff, as in De Oliveira-De Jesus et al. (2026) |
| `periodo.inc` | Capacity tariff block (P1–P6) of each hour | Spanish 6.1TD calendar |
| `PpvuMadridSarah20052023_localtime.inc` | PV output per unit of capacity, representative year for Madrid built from satellite irradiance for 2005–2023 | PVGIS-SARAH, as in De Oliveira-De Jesus et al. (2026) |
| `PluDataCenter.inc` | Load factor of the data center (constant, equal to 1) | — |

### `data/model_outputs/` — results of each optimization run

- **`results_summary_<case>.csv`:** one file per run, holding the annual quantities in the columns *Variable, Value, Unit*. Examples are installed capacities, CAPEX, O&M, purchases $E_b$, capacity charges $CP$, sales $E_s$, the grid-only payments $E_{b0}$ and $CP_0$, NPV, IRR, and LCOEs.
- **The cases are:**
  - `optimal`: the NPV-maximizing design, energy and capacity optimized;
  - `nocontract`: the NPV-maximizing design with energy only optimized, contracted capacity fixed at 1 MW in the six blocks;
  - `ember`: the Ember configuration with grid connection;
  - `ember_island`: the Ember configuration as an island;
  - `subXX`: re-optimization with an investment subsidy of XX %.
- **`hourly/solution_optimal.csv`:** the hourly dispatch of the NPV-maximizing design, in kW. Its columns are:
  - `Pb`: purchases;
  - `Ppv`: PV used;
  - `Pd`: battery discharge;
  - `Ppvmax`: PV available;
  - `Ps`: sales;
  - `Pc`: battery charge;
  - `SOC`: state of charge in kWh;
  - `SpotPrice`: spot price in EUR/kWh.

All files are plain text and machine-readable.

## Regenerating the model outputs (optional)

`model/sizing_model.py` is the sizing and dispatch model of Appendix A: a Big-M MILP with 8 760 hourly periods. It writes the files in `data/model_outputs/`. Run it from the `model/` folder:

```bash
pip install gurobipy highspy numpy numpy-financial
cd model
python sizing_model.py --case optimal
python sizing_model.py --case nocontract
python sizing_model.py --case ember
python sizing_model.py --case ember_island
for s in 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50; do python sizing_model.py --case subsidy --sub $s; done
```

**Solvers.**
- `gurobipy` is used only to build the model and write it in MPS format. The free, size-limited `gurobipy` licence is enough for that step.
- The model is then solved with the open-source solver **HiGHS** (`highspy` 1.16.0, MIT licence). HiGHS solves an LP relaxation first, then the MILP from a start built from the relaxation.

**Run times and gaps.**
- Each run takes 1–2 minutes on a laptop.
- The runs at 40, 45, and 50 % subsidy stopped at the 900 s time limit, with relative MIP gaps of 0.47 %, 0.40 %, and 0.28 %.
- All other runs are proven optimal.

**Non-unique dispatch.** The 2024 Spanish spot price is exactly zero in 532 hours. In those hours, selling or curtailing PV gives the same NPV. Other solvers may therefore return a different split between `Ws` and curtailment, with the same NPV and the same $LCOE_s$.

## Computational cost of the subsidy sweep

`docs/computational_report_subsidy_sweep.pdf` documents the model size, the computing environment, the solve time and MIP gap of each run of Figure 2, and the costs of reproduction. The solver logs of every run are in `data/model_outputs/solver_logs/`.

## Software used

The package was run with:
- Python 3.13;
- numpy 2.5.3, pandas 3.0.5, and matplotlib 3.11.2;
- for the model only: highspy 1.16.0 and gurobipy 14.0.

## Licence

Code: MIT. Data: CC BY 4.0. The original sources of the input data (OMIE, PVGIS) retain their own terms.

## References

De Oliveira-De Jesus, P.M., Salas D., A.J., Yusta, J.M., 2026. Model for optimal sizing of photovoltaic and battery energy storage systems considering an efficient scheme for contracted power capacity. J. Energy Storage 180, 124014.
