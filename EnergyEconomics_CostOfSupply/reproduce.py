#!/usr/bin/env python3
"""
Replication script for
  P. M. De Oliveira-De Jesus, "The cost of supplying a prosumer: levelized costs
  and contract price floors for optimally sized and operated behind-the-meter
  solar and storage", submitted to Energy Economics.

One command reproduces every table, figure and number reported in Sections 4.2-4.5:

    python reproduce.py

The script does NOT run the optimizer. It reads
  data/inputs/          hourly inputs of the sizing model (prices, tariffs, PV profile, load)
  data/model_outputs/   annual results of each optimization run (one CSV per run) and the
                        hourly dispatch of the NPV-maximizing design
and applies the closed-form accounting of Section 3 (Eqs. 1-10, Table 2).
The optimization model that produced data/model_outputs is in model/ (see README.md).

Outputs are written to tables/, figures/ and key_numbers.txt. Every number quoted in
the paper is compared with the recomputed value; the script ends with a PASS/FAIL count.
Requirements: Python >= 3.9, numpy, pandas, matplotlib.
"""
import csv, os, re, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DIN = os.path.join(HERE, "data", "inputs")
DOUT = os.path.join(HERE, "data", "model_outputs")
TAB = os.path.join(HERE, "tables"); FIG = os.path.join(HERE, "figures")
os.makedirs(TAB, exist_ok=True); os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Parameters (Table 4; De Oliveira-De Jesus et al., 2026, J. Energy Storage 180, 124014)
# ---------------------------------------------------------------------------
P = dict(
    i_r=0.077,            # discount rate
    e=0.025,              # escalation of annual cash flows
    n=20,                 # project life [yr]
    c_pv=388.0,           # PV modules [USD/kWp]
    c_bess=185.0,         # battery [USD/kWh]
    c_inv=48.0,           # inverters [USD/kW]
    soft=1.2,             # soft-cost factor (CAPEX = soft x direct costs)
    om_pv=12.5,           # O&M PV [USD/kWp-yr]
    om_bess=5.9,          # O&M battery [USD/kWh-yr]
    usd_per_eur=1.1,      # exchange rate (2024 average)
    cap_tariff=[28.79, 15.08, 6.56, 5.17, 1.93, 0.92],  # Spanish 6.1TD capacity tariffs [EUR/kW-yr], blocks P1..P6
    load_kw=1000.0,       # constant load of the data center [kW]
    contracted_kw=1000.0, # contracted capacity without the project, all blocks [kW]
)
i, e, n = P["i_r"], P["e"], P["n"]
CRF = i * (1 + i) ** n / ((1 + i) ** n - 1)                    # Eq. (1)
CRFe = (i - e) / (1 - ((1 + e) / (1 + i)) ** n)                 # Eq. (1)
K = CRF / CRFe
WL = 8760.0                                                     # MWh/yr

# Carbon parameters (Section 4.5)
EMIS_2024_T = 27.0e6                       # tCO2 from Spanish power generation, 2024 (REE, 2025)
RES_GWH, RES_SHARE = 148_999.0, 0.568      # renewable generation 2024 and its share (REE, 2025)
EF_AVG = EMIS_2024_T / (RES_GWH / RES_SHARE * 1e3)   # tCO2/MWh, average of the system
EF_GAS = 0.37                              # tCO2/MWh, combined-cycle gas plant (REE, 2021)
PC_LIST = [50, 100, 200]                   # external carbon prices [USD/t] of Table 8

# ---------------------------------------------------------------------------
# 2. Helpers
# ---------------------------------------------------------------------------
def read_inc(name):
    """Read a GAMS-style .inc file with lines 't<hour> <value>'."""
    vals = {}
    with open(os.path.join(DIN, name), encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s.startswith("t"):
                p = re.split(r"\s+", s)
                vals[int(p[0][1:])] = float(p[1])
    return np.array([vals.get(t, 0.0) for t in range(1, 8761)])

def summary(name):
    """Annual results of one optimization run (written by the model script)."""
    out = {}
    with open(os.path.join(DOUT, name)) as f:
        for r in csv.DictReader(f):
            try: out[r["Variable"]] = float(r["Value"])
            except ValueError: out[r["Variable"]] = r["Value"]
    return out

def irr(cash):
    """Internal rate of return by bisection."""
    f = lambda r: sum(c / (1 + r) ** k for k, c in enumerate(cash))
    lo, hi = -0.99, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(lo) * f(mid) <= 0: hi = mid
        else: lo = mid
    return (lo + hi) / 2

def levelize(d, sub_rate=0.0):
    """Eqs. (2)-(6): levelized components per MWh of load [USD/MWh] and NPV [USD]."""
    I0, OM = d["CAPEX I0"], d["O&M of assets"]
    Eb, CP, Es = d["Energy purchases Eb"], d["Capacity charges CP"], d["Energy sales Es"]
    OPEX0 = d["Purchases without project Eb0"] + d["Capacity charges without project CP0"]
    SUB = sub_rate * I0
    r = dict(
        PV=d["PV capacity"] / 1e3, C=d["BESS capacity C"] / 1e3, I0=I0, SUB=SUB,
        Wb=d["Energy bought Wb"] / 1e3, Ws=d["Energy sold Ws"] / 1e3,
        G=K * (Eb + CP) / WL, S=K * Es / WL, Z=SUB * CRF / WL,  # Eq. (4)
        L0=K * OPEX0 / WL,                                      # Eq. (6)
    )
    r["LCOE"] = (I0 + OM / CRFe) * CRF / WL                     # USD/MWh
    r["LS"] = r["LCOE"] + r["G"] - r["S"] - r["Z"]              # Eq. (5)
    r["LSsoc"] = r["LCOE"] + r["G"] - r["S"]                    # Remark 1, Z = 0
    r["V"] = r["L0"] - r["LS"]
    CF = OPEX0 + Es - (OM + Eb + CP)                            # annual cash flow (first year)
    r["NPV"] = CF / CRFe - (I0 - SUB)                           # Eq. (2)
    r["NPVsoc"] = CF / CRFe - I0
    r["IRR"] = irr([-(I0 - SUB)] + [CF] * n) * 100
    r["IRRe"] = irr([-(I0 - SUB)] + [CF * (1 + e) ** (y - 1) for y in range(1, n + 1)]) * 100
    r["dW"] = WL - r["Wb"] + r["Ws"]                            # displaced grid energy [MWh/yr]
    r["chk"] = r["V"] - r["NPV"] * CRF / WL                     # Proposition 1 (must be 0)
    r["_d"] = d
    return r

def bounds(r):
    """Table 2: floors and ceilings of the three arrangements [USD/MWh of load]."""
    return {"A": (r["LS"], r["L0"]),
            "B": (r["LCOE"] - r["S"] - r["Z"], r["L0"] - r["G"]),
            "C": (r["LCOE"] - r["Z"], r["L0"] - r["G"] + r["S"])}

def X(r, pc, ef):
    """Eq. (8): levelized external benefit of displaced grid generation [USD/MWh]."""
    return K * pc * ef * r["dW"] / WL

CHECKS = []
def check(label, value, reported, tol):
    ok = abs(value - reported) <= tol
    CHECKS.append((label, value, reported, ok))

def write_table(name, df, caption):
    df.to_csv(os.path.join(TAB, name + ".csv"))
    with open(os.path.join(TAB, name + ".tex"), "w") as f:
        f.write(f"% {caption}\n" + df.to_latex(float_format=lambda x: f"{x:,.2f}"))

# ---------------------------------------------------------------------------
# 3. Grid-only cost from the hourly inputs (cross-check of E_b0 and CP_0)
# ---------------------------------------------------------------------------
lam = read_inc("lambda_spain_localtime.inc")   # day-ahead price [EUR/kWh], 2024 (OMIE)
psi = read_inc("psi.inc")                      # regulated charges on withdrawals [EUR/kWh]
plu = read_inc("PluDataCenter.inc")            # load factor (1 for the data center)
Eb0 = float(np.sum((lam + psi) * P["load_kw"] * plu)) * P["usd_per_eur"]
CP0 = sum(P["cap_tariff"]) * P["contracted_kw"] * P["usd_per_eur"]
check("E_b0 from inputs [USD/yr]", Eb0, 908223, 1)
check("CP_0 from inputs [USD/yr]", CP0, 64295, 1)
check("Hours with zero price", int(np.sum(lam == 0)), 532, 0)
check("Hours with negative price", int(np.sum(lam < 0)), 247, 0)

# ---------------------------------------------------------------------------
# 4. Designs (Table 5)
# ---------------------------------------------------------------------------
opt = levelize(summary("results_summary_optimal.csv"))
noc = levelize(summary("results_summary_nocontract.csv"))
emb = levelize(summary("results_summary_ember.csv"))
isl_d = summary("results_summary_ember_island.csv")
isl = levelize(isl_d)
# unserved load from the energy balance of the island run: W_u = W_L + W_c - w_pv - W_d
Wu = isl_d["Load energy Wl"] + isl_d["BESS charged Wc"] - isl_d["PV used wpv"] - isl_d["BESS discharged Wd"]
served = 1 - Wu / isl_d["Load energy Wl"]

for nm, r in [("NPV-maximizing", opt), ("Fixed capacity", noc), ("Ember", emb)]:
    assert abs(r["chk"]) < 1e-6, f"Proposition 1 fails for {nm}"
    check(f"{nm}: asset LCOE vs model output", r["LCOE"], r["_d"]["Asset LCOE 1 (Wl)"], 1e-6)
    check(f"{nm}: LCOE_s vs model output", r["LS"], r["_d"]["Supply LCOE (Wl)"], 1e-6)
    check(f"{nm}: NPV vs model output [USD]", r["NPV"], r["_d"]["NPV"], 1.0)

cols = {"NPV-maximizing": opt, "Fixed capacity contract": noc, "Ember configuration": emb}
t5 = pd.DataFrame({k: {
    "PV [MWp]": v["PV"], "Battery [MWh]": v["C"], "Purchases [MWh/yr]": v["Wb"],
    "Exports [MWh/yr]": v["Ws"], "CAPEX [MUSD]": v["I0"] / 1e6, "NPV [MUSD]": v["NPV"] / 1e6,
    "IRR constant CF [%]": v["IRR"], "IRR escalated CF [%]": v["IRRe"],
    "Asset LCOE": v["LCOE"], "G": v["G"], "S": v["S"], "LCOE_s": v["LS"],
    "LCOE_0": v["L0"], "V": v["V"]} for k, v in cols.items()})
write_table("table5_designs", t5, "Table 5. Designs and levelized costs (USD/MWh of load unless stated)")

rep5 = {"opt": (6.1, 9.3, 2432, 2011, 5.23, 3.62, 12.7, 15.1, 77.5, 31.9, 17.0, 92.4, 133.6, 41.2),
        "noc": (5.7, 8.1, 2771, 1896, 4.77, 2.91, 11.9, 14.3, 70.7, 45.4, 15.6, 100.5, 133.6, 33.1),
        "emb": (6.0, 17.0, 665, 0, 6.95, 1.88, 8.4, 10.8, 103.1, 9.0, 0.0, 112.2, 133.6, 21.4)}
keys = ["PV", "C", "Wb", "Ws", "I0", "NPV", "IRR", "IRRe", "LCOE", "G", "S", "LS", "L0", "V"]
scale = {"I0": 1e-6, "NPV": 1e-6}
for (lab, r) in [("opt", opt), ("noc", noc), ("emb", emb)]:
    for k, rv in zip(keys, rep5[lab]):
        v = r[k] * scale.get(k, 1)
        tol = 0.6 if k in ("Wb", "Ws") else (0.006 if k in ("I0", "NPV") else 0.051)
        check(f"Table 5 {lab} {k}", v, rv, tol)

check("Ember island: share of load served [%]", served * 100, 92.5, 0.05)
check("Ember island: asset LCOE per MWh served", isl["LCOE"] / served, 112.2, 0.05)
check("Ember: asset LCOE per MWh of load", emb["LCOE"], 103.1, 0.05)
check("LCOE_s / LCOE_0, NPV-maximizing [%]", opt["LS"] / opt["L0"] * 100, 69, 0.5)
check("Share of load bought from the grid [%]", opt["Wb"] / WL * 100, 28, 0.5)
check("Ratio V(opt)/V(Ember)", opt["V"] / emb["V"], 1.92, 0.01)
check("Supply cost below grid, Ember [%]", (1 - emb["LS"] / emb["L0"]) * 100, 16, 0.5)
check("Floor A: optimal minus Ember", emb["LS"] - opt["LS"], 19.8, 0.05)
check("Capacity charges, fixed contract [kUSD/yr]", noc["_d"]["Capacity charges CP"] / 1e3, 64, 0.5)
check("Capacity charges, optimal [kUSD/yr]", opt["_d"]["Capacity charges CP"] / 1e3, 1, 0.05)
alpha_114 = (114 - opt["LS"]) / opt["V"]
check("alpha such that pi_p = 114 (Eq. 11)", alpha_114, 0.52, 0.005)
check("LCOE_s minus asset LCOE, NPV-maximizing", opt["LS"] - opt["LCOE"], 14.9, 0.05)
check("LCOE_s: fixed capacity minus optimal", noc["LS"] - opt["LS"], 8.1, 0.05)
check("Asset LCOE: Ember minus optimal", emb["LCOE"] - opt["LCOE"], 25.6, 0.05)
check("G: optimal minus Ember", opt["G"] - emb["G"], 22.9, 0.05)
check("G: fixed capacity minus optimal", noc["G"] - opt["G"], 13.5, 0.05)
check("LCOE + G (threshold for negative LCOE_s)", opt["LCOE"] + opt["G"], 109.4, 0.05)

# Zero and negative price hours (Section 4.2), hourly dispatch of the NPV-maximizing design
h = pd.read_csv(os.path.join(DOUT, "hourly", "solution_optimal.csv"))
zn = h.SpotPrice <= 0
check("Battery charge at zero/negative prices [MWh]", h.Pc[zn].sum() / 1e3, 757, 0.5)
check("PV curtailed at zero/negative prices [MWh]", (h.Ppvmax - h.Ppv)[zn].sum() / 1e3, 623, 0.5)
check("Purchases at zero prices [MWh]", h.Pb[h.SpotPrice == 0].sum() / 1e3, 77, 0.6)
check("Exports at zero/negative prices [MWh]", h.Ps[zn].sum() / 1e3, 0, 0.01)

# ---------------------------------------------------------------------------
# 5. Contract price bounds (Table 6)
# ---------------------------------------------------------------------------
rows = {}
for nm, r in [("NPV-maximizing", opt), ("Fixed capacity contract", noc), ("Ember configuration", emb)]:
    b = bounds(r)
    for k in "ABC":
        assert abs(b[k][1] - b[k][0] - r["V"]) < 1e-9          # Proposition 2 (i)
    rows[(nm, "Floor")] = {k: b[k][0] for k in "ABC"}
    rows[(nm, "Ceiling")] = {k: b[k][1] for k in "ABC"}
    rows[(nm, "Price, alpha = 0.5")] = {k: b[k][0] + 0.5 * r["V"] for k in "ABC"}
t6 = pd.DataFrame(rows).T.rename(columns={"A": "A. Full-supply PPA", "B": "B. On-site PPA", "C": "C. Asset lease"})
write_table("table6_bounds", t6, "Table 6. Reservation prices [USD/MWh of load], no subsidies")
for (lab, kind, vals) in [("NPV-maximizing", "Floor", (92.4, 60.5, 77.5)), ("NPV-maximizing", "Ceiling", (133.6, 101.7, 118.7)),
                          ("NPV-maximizing", "Price, alpha = 0.5", (113.0, 81.1, 98.1)),
                          ("Fixed capacity contract", "Floor", (100.5, 55.1, 70.7)), ("Fixed capacity contract", "Ceiling", (133.6, 88.2, 103.8)),
                          ("Ember configuration", "Floor", (112.2, 103.1, 103.1)), ("Ember configuration", "Ceiling", (133.6, 124.6, 124.6))]:
    for c, rv in zip(t6.columns, vals):
        check(f"Table 6 {lab} {kind} {c[:2]}", t6.loc[(lab, kind), c], rv, 0.051)
check("Developer NPV with alpha = 0.5 [MUSD]", 0.5 * opt["NPV"] / 1e6, 1.81, 0.006)

# ---------------------------------------------------------------------------
# 6. Investment subsidy (Table 7) and sweep (Figure 2)
# ---------------------------------------------------------------------------
sub_same = levelize(summary("results_summary_optimal.csv"), 0.30)
sub_opt = levelize(summary("results_summary_sub30.csv"), 0.30)
check("Subsidized NPV vs model objective [USD]", sub_opt["NPV"], sub_opt["_d"]["NPV"], 1.0)
t7 = pd.DataFrame({k: {
    "PV [MWp]": v["PV"], "Battery [MWh]": v["C"], "Purchases [MWh/yr]": v["Wb"], "Exports [MWh/yr]": v["Ws"],
    "CAPEX [MUSD]": v["I0"] / 1e6, "Subsidy [MUSD]": v["SUB"] / 1e6, "Asset LCOE": v["LCOE"],
    "G": v["G"], "S": v["S"], "Z": v["Z"], "Private LCOE_s": v["LS"], "Social LCOE_s (Z=0)": v["LSsoc"],
    "Private NPV [MUSD]": v["NPV"] / 1e6, "Social NPV [MUSD]": v["NPVsoc"] / 1e6,
    "Floor A": bounds(v)["A"][0], "Ceiling A": bounds(v)["A"][1], "Floor B": bounds(v)["B"][0],
    "Ceiling B": bounds(v)["B"][1], "Floor C": bounds(v)["C"][0], "Ceiling C": bounds(v)["C"][1]}
    for k, v in {"No subsidy": opt, "Subsidy, same design": sub_same, "Subsidy, re-optimized": sub_opt}.items()})
write_table("table7_subsidy", t7, "Table 7. Effect of a 30 % investment subsidy")
rep7 = {"Subsidy, same design": {"Z": 17.8, "Private LCOE_s": 74.6, "Private NPV [MUSD]": 5.19, "Floor B": 42.7, "Floor C": 59.7,
                                  "Subsidy [MUSD]": 1.57},
        "Subsidy, re-optimized": {"PV [MWp]": 8.5, "Battery [MWh]": 17.1, "Purchases [MWh/yr]": 962, "Exports [MWh/yr]": 3356,
                                   "CAPEX [MUSD]": 8.20, "Subsidy [MUSD]": 2.46, "Asset LCOE": 121.6, "G": 12.4, "S": 35.5,
                                   "Z": 28.0, "Private LCOE_s": 70.6, "Social LCOE_s (Z=0)": 98.6, "Private NPV [MUSD]": 5.54,
                                   "Social NPV [MUSD]": 3.08, "Floor B": 58.2, "Ceiling B": 121.2, "Floor C": 93.7, "Ceiling C": 156.6}}
for col, d in rep7.items():
    for k, rv in d.items():
        tol = 0.6 if "MWh/yr" in k else (0.006 if "MUSD" in k else 0.051)
        check(f"Table 7 {col} {k}", t7.loc[k, col], rv, tol)
DWL = sub_opt["NPVsoc"] - opt["NPVsoc"]
DWL = -DWL
check("Deadweight loss [MUSD], Eq. (9)", DWL / 1e6, 0.54, 0.006)
check("Deadweight loss / subsidy [%]", DWL / sub_opt["SUB"] * 100, 22, 0.5)
check("Private gain [MUSD]", (sub_opt["NPV"] - opt["NPV"]) / 1e6, 1.91, 0.006)
check("DWL via Eq. (9) equals NPV difference", (sub_opt["LSsoc"] - opt["LSsoc"]) * WL / CRF, DWL, 1e-3)

# Subsidy sweep
sweep = [(0, opt)]
for rate in range(5, 55, 5):
    fn = f"results_summary_sub{rate:02d}.csv"
    if os.path.exists(os.path.join(DOUT, fn)):
        sweep.append((rate, levelize(summary(fn), rate / 100)))
sw = []
for rate, r in sweep:
    dwl = opt["NPVsoc"] - r["NPVsoc"]
    ddw = r["dW"] - opt["dW"]
    sw.append(dict(rate=rate, PV=r["PV"], C=r["C"], NPVpriv=r["NPV"] / 1e6, NPVsoc=r["NPVsoc"] / 1e6,
                   subsidy=r["SUB"] / 1e6, DWL=dwl / 1e6, DWL_share=dwl / r["SUB"] if r["SUB"] else np.nan,
                   pstar_avg=dwl * CRFe / (ddw * EF_AVG) if rate else np.nan,
                   pstar_gas=dwl * CRFe / (ddw * EF_GAS) if rate else np.nan))
sw = pd.DataFrame(sw).set_index("rate")
write_table("figure2_subsidy_sweep", sw, "Data of Figure 2: subsidy-rate sweep")
check("Sweep 10 %: DWL [MUSD]", sw.loc[10, "DWL"], 0.04, 0.005)
check("Sweep 10 %: DWL share [%]", sw.loc[10, "DWL_share"] * 100, 6, 0.5)
check("Sweep 50 %: DWL [MUSD]", sw.loc[50, "DWL"], 2.0, 0.05)
check("Sweep 50 %: DWL share [%]", sw.loc[50, "DWL_share"] * 100, 34, 0.5)
check("Sweep 50 %: battery [MWh]", sw.loc[50, "C"], 25.8, 0.05)
check("Sweep 50 %: PV [MWp]", sw.loc[50, "PV"], 11.5, 0.05)
check("Sweep 5 %: p* average EF [USD/t]", sw.loc[5, "pstar_avg"], 36, 0.5)
check("Sweep 50 %: p* average EF [USD/t]", sw.loc[50, "pstar_avg"], 268, 0.5)
check("Sweep 5 %: p* gas [USD/t]", sw.loc[5, "pstar_gas"], 10, 0.5)
check("Sweep 50 %: p* gas [USD/t]", sw.loc[50, "pstar_gas"], 75, 0.5)
# subsidy rate at which p* (average EF) equals 140 USD/t, by linear interpolation
xs, ys = sw.index[1:].values, sw["pstar_avg"].values[1:]
rate140 = float(np.interp(140, ys, xs))
check("Subsidy rate with p* = 140 USD/t [%]", rate140, 27, 0.5)

# ---------------------------------------------------------------------------
# 7. Carbon externality (Table 8, Proposition 3, Figure 3)
# ---------------------------------------------------------------------------
check("Average emission factor [t/MWh]", EF_AVG, 0.103, 0.0005)
dW = sub_opt["dW"] - opt["dW"]
check("Displaced energy, NPV-maximizing [MWh/yr]", opt["dW"], 8339, 0.5)
check("Displaced energy, subsidized [MWh/yr]", sub_opt["dW"], 11153, 0.5)
check("Additional displaced energy [MWh/yr]", dW, 2814, 0.5)
thr = DWL * CRFe / dW                                        # Eq. (10), USD/MWh
check("Threshold p_c*EF [USD/MWh], Eq. (10)", thr, 16.0, 0.05)
check("Threshold, average EF [USD/t]", thr / EF_AVG, 156, 0.5)
check("Threshold, gas EF [USD/t]", thr / EF_GAS, 43, 0.5)
t8 = {}
for lab, ef, pcs in [("No carbon", 0, [0]), ("Average EF", EF_AVG, PC_LIST), ("Gas EF", EF_GAS, PC_LIST)]:
    for pc in pcs:
        row = {}
        for nm, r in [("NPV-maximizing", opt), ("Subsidized, re-optimized", sub_opt), ("Ember", emb)]:
            Ls = r["LSsoc"] - X(r, pc, ef)
            row[nm + " LCOE_s"] = Ls
            row[nm + " NPV [MUSD]"] = (r["L0"] - Ls) * WL / CRF / 1e6
        t8[(lab, pc)] = row
t8 = pd.DataFrame(t8).T
write_table("table8_carbon", t8, "Table 8. Social cost of supply with a carbon externality")
rep8 = {("No carbon", 0): (92.4, 3.62, 98.6, 3.08, 112.2, 1.88),
        ("Average EF", 50): (86.5, 4.14, 90.7, 3.77, 106.4, 2.39), ("Average EF", 100): (80.6, 4.66, 82.8, 4.47, 100.7, 2.89),
        ("Average EF", 200): (68.8, 5.70, 67.0, 5.85, 89.3, 3.90), ("Gas EF", 50): (71.2, 5.49, 70.2, 5.57, 91.6, 3.69),
        ("Gas EF", 100): (50.0, 7.35, 41.9, 8.07, 71.0, 5.50), ("Gas EF", 200): (7.6, 11.08, -14.8, 13.05, 29.9, 9.12)}
for k, vals in rep8.items():
    for c, rv in zip(t8.columns, vals):
        check(f"Table 8 {k} {c}", t8.loc[k, c], rv, 0.006 if "[MUSD]" in c else 0.051)

# ---------------------------------------------------------------------------
# 8. Figures
# ---------------------------------------------------------------------------
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 0.6, "legend.frameon": False})
C_COST, C_UP, C_DN, C_SUB, C_RES = "#4a4a4a", "#c0504d", "#4f81bd", "#9bbb59", "#1f1f1f"

# Figure 1: decomposition of the cost of supply
fig, axs2 = plt.subplots(2, 2, figsize=(3.5, 3.4), sharey=True); axs = axs2.ravel()
L0 = opt["L0"]
for ax, (name, r) in zip(axs, [("(a) NPV-maximizing", opt), ("(b) Fixed capacity contract", noc),
                               ("(c) Ember configuration", emb), ("(d) Subsidized, re-optimized", sub_opt)]):
    LA, G, S, Z = r["LCOE"], r["G"], r["S"], r["Z"]; LS = r["LS"]
    steps = [("$LCOE$", 0, LA, C_COST), ("$+G$", LA, LA + G, C_UP), ("$-S$", LA + G, LA + G - S, C_DN)]
    if Z > 0: steps.append(("$-Z$", LA + G - S, LS, C_SUB))
    steps.append(("$LCOE_s$", 0, LS, C_RES)); m = len(steps)
    for k, (lab, a, b, col) in enumerate(steps):
        ax.bar(k, b - a, bottom=a, width=0.7, color=col, edgecolor="white", linewidth=0.6)
        if k < m - 1: ax.plot([k + 0.33, k + 0.67], [b, b], color="0.55", lw=0.5)
        if k in (0, m - 1): ax.text(k, b + 1.5, f"{b:.1f}", ha="center", va="bottom", fontsize=5)
        elif abs(b - a) < 0.05: ax.text(k, b + 1.5, "0", ha="center", va="bottom", fontsize=5)
        elif b > a: ax.text(k, b + 1.5, f"+{b - a:.1f}", ha="center", va="bottom", fontsize=5)
        else: ax.text(k, b - 1.5, f"−{a - b:.1f}", ha="center", va="top", fontsize=5)
    ax.axhline(L0, ls="--", lw=0.7, color="0.35")
    xa = m - 1 + 0.5
    ax.annotate("", xy=(xa, LS), xytext=(xa, L0), arrowprops=dict(arrowstyle="<->", lw=0.6, color="0.3", shrinkA=0, shrinkB=0))
    ax.text(xa + 0.12, (LS + L0) / 2, f"$V$\n{L0 - LS:.1f}", fontsize=5, va="center", ha="left")
    if Z > 0:
        ax.plot([m - 1 - 0.33, m - 1 + 0.33], [r["LSsoc"]] * 2, color=C_UP, lw=1.4)
        ax.text(m - 1, r["LSsoc"] + 1.5, f"social\n{r['LSsoc']:.1f}", fontsize=4.6, ha="center", va="bottom", color=C_UP)
    ax.set_xticks(range(m)); ax.set_xticklabels([s[0] for s in steps], fontsize=5.5)
    ax.set_title(name, fontsize=6.2); ax.set_xlim(-0.6, m - 1 + 1.25)
axs[0].set_ylabel("USD/MWh of load", fontsize=6); axs[2].set_ylabel("USD/MWh of load", fontsize=6)
for a_ in axs: a_.tick_params(labelsize=5.5, length=2); axs[0].set_ylim(0, 150)
axs[0].text(-0.5, L0 + 1.5, "$LCOE_0$", fontsize=5)
fig.tight_layout(w_pad=0.4, h_pad=0.6, pad=0.3); fig.savefig(os.path.join(FIG, "fig_decomposition.pdf")); plt.close(fig)

# Figure 2: subsidy sweep
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.3)); x = sw.index.values
a = axs[0]
a.plot(x, sw.PV, "o-", color=C_DN, ms=2.5, lw=1.0, label="PV [MWp]")
a.plot(x, sw.C, "s--", color=C_UP, ms=2.5, lw=1.0, label="Battery [MWh]")
a.axhline(17, color="0.6", lw=0.6, ls=":"); a.text(1, 17.5, "Ember battery", fontsize=5, color="0.4")
a.set_title("(a) Design", fontsize=6.5); a.legend(fontsize=5.5, loc="upper left")
a = axs[1]
a.plot(x, sw.NPVpriv, "o-", color=C_DN, ms=2.5, lw=1.0, label="Private NPV")
a.plot(x, sw.NPVsoc, "s--", color=C_COST, ms=2.5, lw=1.0, label="Social NPV")
a.plot(x, sw.subsidy, "^-", color=C_SUB, ms=2.5, lw=0.9, label="Subsidy outlay")
a.plot(x, sw.DWL, "d-", color=C_UP, ms=2.5, lw=1.0, label="Deadweight loss")
a.set_ylim(-0.3, 9.8); a.set_title("(b) Value [million USD]", fontsize=6.5); a.legend(fontsize=5.5, loc="upper left")
a = axs[2]; xs = x[1:]
a.plot(xs, sw.pstar_avg.values[1:], "o-", color=C_DN, ms=2.5, lw=1.0, label="Average EF, 0.103 t/MWh")
a.plot(xs, sw.pstar_gas.values[1:], "s--", color=C_UP, ms=2.5, lw=1.0, label="Gas-fired, 0.37 t/MWh")
a.axhspan(130, 150, color="0.85", lw=0); a.text(1, 140, "~140 USD/t", fontsize=5, color="0.3", va="center")
a.set_xlim(0, 52); a.set_ylim(0, 300); a.set_title("(c) Break-even $p_c$ [USD/t]", fontsize=6.5)
a.legend(fontsize=5.5, loc="upper left")
for a in axs:
    a.set_xticks([0, 10, 20, 30, 40, 50]); a.set_xlabel("Investment subsidy rate [%]", fontsize=6)
    a.tick_params(labelsize=5.5, length=2)
fig.tight_layout(w_pad=1.0); fig.savefig(os.path.join(FIG, "fig_sweep.pdf")); plt.close(fig)

# Figure 3: social cost of supply vs carbon price
pc = np.linspace(0, 250, 251)
fig, axs = plt.subplots(2, 1, figsize=(3.5, 3.9), sharex=True)
for ax, (ef, lab) in zip(axs, [(EF_AVG, "(a) Average EF, 0.103 t/MWh"), (EF_GAS, "(b) Gas-fired, 0.37 t/MWh")]):
    for r, name, col, ls in [(opt, "NPV-maximizing, no subsidy", C_DN, "-"), (sub_opt, "Subsidized, re-optimized", C_UP, "--")]:
        ax.plot(pc, r["LSsoc"] - K * pc * ef * r["dW"] / WL, color=col, ls=ls, lw=1.1, label=name)
    ps = thr / ef; ys = opt["LSsoc"] - X(opt, ps, ef)
    ax.plot([ps], [ys], "o", ms=3, color="k"); ax.axvline(ps, lw=0.5, color="0.6", ls=":")
    ax.text(ps + 4, ys + 6, f"{ps:.0f} USD/t", fontsize=5.5); ax.axhline(0, lw=0.5, color="0.6")
    ax.set_title(lab, fontsize=6.5); ax.set_xlim(0, 250)
for ax in axs:
    ax.set_ylabel("Social cost of supply [USD/MWh]", fontsize=6); ax.tick_params(labelsize=5.5, length=2)
axs[1].set_xlabel("Carbon price $p_c$ [USD/tCO$_2$]", fontsize=6)
axs[0].legend(fontsize=5.5, loc="lower left")
fig.tight_layout(h_pad=0.8, pad=0.3); fig.savefig(os.path.join(FIG, "fig_carbon.pdf")); plt.close(fig)

# ---------------------------------------------------------------------------
# 9. Report
# ---------------------------------------------------------------------------
with open(os.path.join(HERE, "key_numbers.txt"), "w") as f:
    f.write(f"CRF = {CRF:.6f}, CRF_e = {CRFe:.6f}, CRF/CRF_e = {K:.5f}\n")
    f.write(f"{'Check':60s} {'Computed':>14s} {'Reported':>12s}  Result\n")
    for lab, v, rv, ok in CHECKS:
        f.write(f"{lab:60s} {v:14.4f} {rv:12.4f}  {'PASS' if ok else 'FAIL'}\n")
npass = sum(c[3] for c in CHECKS)
print(open(os.path.join(HERE, "key_numbers.txt")).read())
print(f"{npass}/{len(CHECKS)} reported numbers reproduced. Tables in tables/, figures in figures/.")
sys.exit(0 if npass == len(CHECKS) else 1)
