#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat March 20 2026
@authors: Paulo De Oliveira (PDEOLIV@GMAIL.COM) y Alejandro Salas Durán
Optimal Sizing + BESS + PV  ─ Case 8760h
MILP version: the bilinear constraints (11)-(12), P_c <= P_bess^inv * w and
P_d <= P_bess^inv * (1 - w), are replaced by linear Big-M constraints.
CASE = 'optimal' : NPV-maximizing sizing (J. Energy Storage 180 (2026) 124014)
CASE = 'ember'   : Ember configuration (6000 kWp, 17000 kWh, no exports, imports <= 12 % of the load)
CASE = 'ember_island' : Ember configuration in island (no purchases, no sales); unserved load allowed and
                  minimized (penalty VOLL); cyclic SOC constraint active. Reproduces the Ember report LCOE
                  when the asset cost is divided by the served energy.
"""
import numpy as np
import numpy_financial as npf
import gurobipy as gp
from gurobipy import GRB, quicksum
import re
import csv
def read_inc(path):
    values = {}
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                if not line.lstrip().startswith("t"):
                    continue
                parts = re.split(r"\s+", line.strip())
                if len(parts) >= 2:
                    hour  = int(parts[0][1:])           
                    value = float(parts[1])
                    values[hour] = value
    except FileNotFoundError:
        return {t: 0.1 for t in range(1, 8761)}
    return values 
paths = {
#     Generic Industries:
#     'Plu'    :  'Plu2.inc',  # Industry load profile (Spanish) 2
#     'Plu'     : 'Plu.inc',    # Industry load profile  (Korean) 1
      'Plu'     : 'PluDataCenter.inc',    # DataCenter load factor = 1
# ----------------------
#    Spanish Electricity Market    
       'lambda'  : 'lambda_spain_localtime.inc', # 8760h 2024 Spanish Spot prices Eur/MWh
       'psi'     : 'psi.inc', # Tariff of Use of the network Eur/MWh
       'Ppvu'    : 'PpvuMadridSarah20052023_localtime.inc',  # Perfil Solar Madrid España Sarah 2005-2023
       "periodo" : "periodo.inc",  # calendar to assingn type of contracted power 1..6 per each hour
# ----------------------
}
series = {name: read_inc(route) for name, route in paths.items()}
T = range(1, 8761)
data = {t: {'lambda': series['lambda'].get(t, 0.0), 'Plu': series['Plu'].get(t, 0.0),
            'psi': series['psi'].get(t, 0.0), 'Ppvu': series['Ppvu'].get(t, 0.0)} for t in T}
periodo = {t: int(series["periodo"].get(t, 6)) for t in T}
CASE = 'optimal'                    # 'optimal', 'ember' or 'ember_island'
Plinst  = 1000.0                    # kW peak/installed Load DataCenter
#Plinst  = 564.3                    # kW peak/installed load Industry 1
#Plinst  = 630.620646500000         # kW peak/installed load Industry 2
Rmax   = 1000                       # W/m2
Area   = 2.4                        # Area module 500W
eta    = 0.2094                     # PV module efficiency (new modules)
eff_c  = 0.9624                     # BESS charge inverter efficiency
eff_d  = 0.9624                     # BESS discharge inverter efficiency
DoD    = 0.90                       # Depth of Discharge
PmaxF   = Plinst                    # límit frontier export/import
er = 1.1 # (exchange rate 1 Eur = 1.1 USD, average 2024)
BoP = 0  # Balance of Plant, USD 
Sc=1.2 # Soft costs, %
OaMpv = 12.5  #USD/kWp
OaMbess = 5.9    #USD/kWp
CAPEX_pv = 388 #USD/kWp
CAPEX_BESS = 185 #USD/kWh
CAPEX_inverter = 48 #USD/kWp
i=7.7/100 #annual discount rate (Lazard)
n=20 #project life span
e=2.5/100 # escalation rate
ir=(i-e)/(1+e) #equivalent annual discount rate including escalation
crf=(i*(i+1)**n)/((i+1)**n-1) #capital recovery factor
crfe=(1+e)*(ir*(ir+1)**n)/((ir+1)**n-1) #capital recovery factor modified
# Big-M for constraints (11)-(12): an upper bound of the BESS inverter power (kW).
# M must not cut the optimum: the script checks after solving that PinverterBESS < M_bess.
M_bess = 5*Plinst
# ── Capacity charges in Spain (USD/kW-year)
kappa = [0, 28.79187*er, 15.07764*er, 6.55917*er, 5.17209*er, 1.93281*er, 0.91609*er] 
# ── Backup charge for self-generators (USD/kW-year), applied to PmaxF. 0 in Spain.
# Colombia (AGGE, Res. CREG 015/2018): 2800 COP/kW-month * 12 / 3915 COP/USD = 8.58 USD/kW-year
Cresp = 0.0
Cresp_year = Cresp*PmaxF            # USD/year
all_results = []    
# ─────────────────────────────────────────────────────────────────────────────
# The model
# ─────────────────────────────────────────────────────────────────────────────
m = gp.Model('SizingModel')
#m.setParam('OutputFlag', 0)  
# Variables
Savings = m.addVar(name='Savings', lb=-GRB.INFINITY); OPEX = m.addVar(name='OPEX', lb=0)
CapacityP = m.addVar(name='Capacity', lb=0); CapacityP0 = m.addVar(name='Capacity0', lb=0)
npv_var = m.addVar(name='npv', lb=-GRB.INFINITY); Benefit = m.addVar(name='Benefit', lb=-GRB.INFINITY)
OPEX0 = m.addVar(name='OPEX0', lb=-GRB.INFINITY); Eb = m.addVar(name='Eb', lb=0); Eb0 = m.addVar(name='Eb0', lb=0)
Es = m.addVar(name='Es', lb=0); wpvmx = m.addVar(name='wpvmx', lb=-GRB.INFINITY); wpv = m.addVar(name='wpv', lb=-GRB.INFINITY)
wclipping = m.addVar(name='wclipping', lb=-GRB.INFINITY); Wb = m.addVar(name='Wb', lb=-GRB.INFINITY)
Ws = m.addVar(name='Ws', lb=-GRB.INFINITY); Wl = m.addVar(name='Wl', lb=-GRB.INFINITY)
Wd = m.addVar(name='Wd', lb=-GRB.INFINITY); Wc = m.addVar(name='Wc', lb=-GRB.INFINITY)
CashFlow_var = m.addVar(name='CashFlow', lb=-GRB.INFINITY); Investment0 = m.addVar(name='Investment', lb=0)
PinverterBESS = m.addVar(name='PinverterBESS', lb=0); PinverterPV = m.addVar(name='PinverterPV', lb=0);Ppvinst = m.addVar(name='Ppvinst', lb=0)
SOC0 = m.addVar(name='SOC0', lb=0); nx = m.addVar(name='nx', lb=0); C = m.addVar(name='C', lb=0)
Pbmax = {p: m.addVar(lb=0, name=f'PbmaxP{p}') for p in range(1, 7)}
SOC = m.addVars(T, lb=0, name='SOC'); Ppv = m.addVars(T, lb=0, name='Ppv')
Ppvmx = m.addVars(T, lb=0, name='Ppvmx'); Pd = m.addVars(T, lb=0, name='Pd')
Pc = m.addVars(T, lb=0, name='Pc'); Pb = m.addVars(T, lb=0, name='Pb'); Ps = m.addVars(T, lb=0, name='Ps')
w1 = m.addVars(T, vtype=GRB.BINARY, name='w1'); w3 = m.addVars(T, vtype=GRB.BINARY, name='w3')
U = m.addVars(T, lb=0, name='U'); Wu = m.addVar(name='Wu', lb=0)   # unserved load (ember_island case only)
# Constraints
m.addConstrs((Pbmax[periodo[t]] >= Pb[t] for t in T), "res_Pbmax")
m.addConstr(Pbmax[6] == PmaxF)
for p in range(1, 6): m.addConstr(Pbmax[p+1] >= Pbmax[p])
for t in T:
    m.addConstr(Pd[t] + Pb[t] + Ppv[t] + U[t] == Pc[t] + Ps[t] + Plinst * data[t]['Plu'])
    m.addConstr(Ppvmx[t] == Ppvinst  * data[t]['Ppvu']) 
    if t == 1:
        m.addConstr(SOC[t] == SOC0 + Pc[t]*eff_c - Pd[t]/eff_d)
    else:
        m.addConstr(SOC[t] == SOC[t-1] + Pc[t]*eff_c - Pd[t]/eff_d)
    m.addConstr(Ppv[t] <= Ppvmx[t])
    m.addConstr(Pc[t] <= PinverterBESS)          # (11a) inverter limit, charge
    m.addConstr(Pc[t] <= M_bess * w1[t])         # (11b) Big-M: charge only if w1 = 1
    m.addConstr(Pd[t] <= PinverterBESS)          # (12a) inverter limit, discharge
    m.addConstr(Pd[t] <= M_bess * (1-w1[t]))     # (12b) Big-M: discharge only if w1 = 0
    m.addConstr(Pb[t] <= PmaxF * w3[t])
    m.addConstr(Ps[t] <= PmaxF * (1-w3[t]))
    m.addConstr(SOC[t] <= ((1-DoD)/2 + DoD)*C)
    m.addConstr(SOC[t] >= ((1-DoD)/2)*C)
m.addConstr(SOC0 <= ((1-DoD)/2 + DoD)*C)
m.addConstr(SOC0 >= ((1-DoD)/2)*C)
# Optional cyclic SOC constraint: the end-of-year SOC must recover the initial SOC, so the initial charge
# is not counted as free energy every year. Commented out to reproduce the published results
# (J. Energy Storage 180 (2026) 124014). Enabling it changes the Madrid NPV by -0.2%.
#m.addConstr(SOC[T[-1]] >= SOC0, name='soc_cyclic')
# FINANCIAL EQUATIONS
m.addConstr(Es == er*(quicksum(data[t]['lambda'] * Ps[t] for t in T)))
m.addConstr(Eb == er*(quicksum((data[t]['lambda'] + data[t]['psi']) * Pb[t] for t in T)))
m.addConstr(Eb0 == er*(quicksum((data[t]['lambda'] + data[t]['psi']) * Plinst * data[t]['Plu'] for t in T)))
m.addConstr(CapacityP == quicksum(kappa[p]*Pbmax[p] for p in range(1, 7)))
m.addConstr(CapacityP0 == sum(kappa[p]*PmaxF for p in range(1, 7)))
m.addConstr(OPEX0 == Eb0 + CapacityP0)
m.addConstr(OPEX == CapacityP + Eb + OaMpv*Ppvinst + OaMbess*C + Cresp_year)
m.addConstr(Savings == OPEX0 - OPEX)
m.addConstr(Benefit == Es - OPEX)
m.addConstr(CashFlow_var == Es + OPEX0 - OPEX)
m.addConstr(Investment0 == BoP + Sc * (CAPEX_pv*Ppvinst + CAPEX_BESS*C + CAPEX_inverter*(PinverterBESS+PinverterPV)))
m.addConstr(npv_var == 1*CashFlow_var/(crfe) - 1*Investment0)
# ENERGY DISPATCH
m.addConstr(wpv == quicksum(Ppv[t] for t in T))
m.addConstr(wpvmx == quicksum(Ppvmx[t] for t in T))
m.addConstr(Wb == quicksum(Pb[t] for t in T))
m.addConstr(Ws == quicksum(Ps[t] for t in T))
m.addConstr(Wc == quicksum(Pc[t] for t in T))
m.addConstr(Wd == quicksum(Pd[t] for t in T))
m.addConstr(Wl == sum(Plinst * data[t]['Plu'] for t in T))
m.addConstr(wclipping == wpvmx - wpv)
m.addConstr(PinverterBESS <= C*2.0)
m.addConstr(PinverterBESS >= C*0.1)
m.addConstr(nx == 1000*Ppvinst/(Rmax*eta*Area))
m.addConstr(Wu == quicksum(U[t] for t in T))
VOLL = 100.0 if CASE == 'ember_island' else 0.0  # USD/kWh penalty on unserved load
m.setObjective(npv_var - VOLL*Wu, GRB.MAXIMIZE)
#m.addConstr(Wb <= 0.12*Wl) 
#m.addConstr(Ws == 0.00*Wl)
# m.addConstr(C == 17000, name='r10')
# m.addConstr(PinverterBESS == 2, name='r11')
#m.addConstr(Ppvinst    == 6000, name='r12')
# m.addConstr(Pbmax[1]   == PmaxF, name='r13')
# m.addConstr(Pbmax[2]   == PmaxF, name='r14')
# m.addConstr(Pbmax[3]  == PmaxF, name='r15')
# m.addConstr(Pbmax[4]  == PmaxF, name='r16')
# m.addConstr(Pbmax[5]   == PmaxF, name='r17')
# m.addConstr(Pbmax[6]  == PmaxF, name='r18')
m.addConstr(PinverterPV == Plinst+PmaxF+PinverterBESS)   
if CASE == 'ember_island':          # Ember in island: fixed sizes, no grid trading, unserved load allowed
    m.addConstr(C == 17000, name='ember_C'); m.addConstr(Ppvinst == 6000, name='ember_PV')
    m.addConstr(Wb == 0, name='island_nobuy'); m.addConstr(Ws == 0, name='island_nosell')
    m.addConstr(SOC[T[-1]] >= SOC0, name='soc_cyclic')
else:
    m.addConstrs((U[t] == 0 for t in T), name='no_unserved')
if CASE == 'ember':                 # Ember benchmark: fixed sizes, no exports, deficit <= 12 % of the load
    m.addConstr(C == 17000, name='ember_C')
    m.addConstr(Ppvinst == 6000, name='ember_PV')
    m.addConstr(Ws == 0.00*Wl, name='ember_noexport')
    m.addConstr(Wb <= 0.12*Wl, name='ember_deficit')
m.optimize()
if m.Status == GRB.OPTIMAL and PinverterBESS.X >= M_bess - 1e-6:
    print(f"WARNING: PinverterBESS = {PinverterBESS.X:,.2f} kW reached M_bess = {M_bess:,.2f} kW. Increase M_bess.")
Crate = PinverterBESS.X/C.X #resulting BESS Crate
BESSinverterCost = CAPEX_inverter*PinverterBESS.X
PVinverterCost = CAPEX_inverter*PinverterPV.X
BESSbatteryCost = CAPEX_BESS*C.X
PVsystemCost =CAPEX_pv*Ppvinst.X+CAPEX_inverter*PinverterPV.X
ILR= Ppvinst.X/PinverterPV.X  #MWdc/MWac, DC-to-AC Ratio or Inverter Loading Ratio (ILR)
if m.Status == GRB.OPTIMAL:
    Io = Investment0.X
    CF = CashFlow_var.X
    TIR = npf.irr([-Io] + [CF]*n) * 100 if Io > 0 else 0
    TIRe = npf.irr([-Io] + [CF*(1+e)**t for t in range(n)]) * 100 if Io > 0 else 0 # IRR with cash flow escalated at e (consistent with NPV)
    BCratio = npf.pv(rate=i, nper=n, pmt=-CF, fv=0) / Io if Io > 0 else 0 # constant cash flow
    BCratioe = (CF/crfe) / Io if Io > 0 else 0 # escalated cash flow (= 1 + NPV/Io)
    # ── Levelized costs (USD/MWh) ──────────────────────────────────────────
    OPEXgross = OaMpv*Ppvinst.X + OaMbess*C.X                       # O&M of the assets only
    Egrid = Eb.X + CapacityP.X + Cresp_year                          # grid payments: energy, capacity, backup
    # (1) Asset LCOE (stand-alone assets: CAPEX + O&M), three energy bases
    LCOEgross  = 1000*(Io + OPEXgross/crfe)/(Wl.X/crf)    if Wl.X    > 0 else 0  # load energy W_L (eq. 8 of the paper)
    LCOEgross2 = 1000*(Io + OPEXgross/crfe)/(wpv.X/crf)   if wpv.X   > 0 else 0  # PV energy used, after clipping (depends on dispatch)
    LCOEgross3 = 1000*(Io + OPEXgross/crfe)/(wpvmx.X/crf) if wpvmx.X > 0 else 0  # PV energy generated, before clipping
    # (2) Supply LCOE: asset LCOE + levelized net grid cost (purchases + capacity + backup - sales), load energy
    LCOEsupply = LCOEgross + 1000*crf*(Egrid - Es.X)/(crfe*Wl.X) if Wl.X > 0 else 0
    # (0) Grid-only LCOE: cost of buying all the load from the grid (no project)
    LCOE0 = 1000*OPEX0.X*crf/(crfe*Wl.X) if Wl.X > 0 else 0
    # Net LCOEs (load energy). Net 1 = -NPV*crf/W_L = LCOEsupply - LCOE0; Net 2 = LCOEsupply (same expression)
    LCOEnet1 = 1000*(Io + (OPEX.X - OPEX0.X - Es.X)/crfe)/(Wl.X/crf) if Wl.X > 0 else 0
    LCOEnet2 = 1000*(Io + (OPEX.X           - Es.X)/crfe)/(Wl.X/crf) if Wl.X > 0 else 0
    NPERaprox = Io/CF if CF > 0 else 0
    try:
        NPER = np.log((CF - i * 0) / (CF + i * (-Io))) / np.log(1 + i) if (CF + i * (-Io)) > 0 else 0 # discounted payback, constant cash flow
    except: NPER = 0
    try:
        NPERe = np.log(1 - Io*(i-e)/CF) / np.log((1+e)/(1+i)) if (CF > 0 and 1 - Io*(i-e)/CF > 0) else 0 # discounted payback, escalated cash flow
    except: NPERe = 0
  
# ─────────────────────────────────────────────────────────────────────────────
if m.Status == GRB.OPTIMAL:
    print(f"Results:")
    print(f"------ System sizing  ------------------")
    print(f"BESS capacity (WbessInst): {C.X:,.2f} kWh")
    print(f"BESS inverter capacity (PbessInst): {PinverterBESS.X:,.2f} kW")
    print(f"BESS C-rate:  {Crate:,.2f} 1/h,  {1/Crate:,.2f}  hours")
    print(f"PV System capacity (PpvInst): {Ppvinst.X:,.2f} kWdc")
    print(f"PV System inverter (Ppv_inverter): {PinverterPV.X:,.2f} kWpv")
    print(f"PV inverter load ratio (ILR, MWdc/MWp): {ILR:,.2f} ")
    print(f"Npan: {nx.X:,.2f} modules of 500W")
    print(f"Contracted power per period")
    print(f"PbmaxP1: {Pbmax[1].X:,.2f} kW")
    print(f"PbmaxP2: {Pbmax[2].X:,.2f} kW") 
    print(f"PbmaxP3: {Pbmax[3].X:,.2f} kW") 
    print(f"PbmaxP4: {Pbmax[4].X:,.2f} kW") 
    print(f"PbmaxP5: {Pbmax[5].X:,.2f} kW") 
    print(f"PbmaxP6: {Pbmax[6].X:,.2f} kW")       
    print(f"------Energy Dispatch--------------------")
    print(f"Energy load consumption (Wl): {Wl.X:,.2f} kWh/year")
    print(f"Energy bought from the market (Wb): {Wb.X:,.2f} kWh/year")
    print(f"Energy sold to the market (Ws): {Ws.X:,.2f} kWh/year")
    print(f"BESS Energy charged (Wc): {Wc.X:,.2f} kWh/year")
    print(f"BESS Energy discharged (Wd): {Wd.X:,.2f} kWh/year")
    print(f"PV energy generated (wpvmx): {wpvmx.X:,.2f} kWh/year")
    print(f"PV energy used after clipping (wpv): {wpv.X:,.2f} kWh/year")
    print(f"PV energy curtailed (Clipping) (Wclipping): {wclipping.X:,.2f} kWh/year")
    print(f"Initial SOC:  {SOC0.X:,.2f} kWh")
    print(f"Final SOC:    {SOC[T[-1]].X:,.2f} kWh")
    print(f"------Financial results--------------------")   
    print(f"Operational Benefit: {Benefit.X:,.2f} USD/year")
    print(f"Operational Expenditure with the project (OPEX): {OPEX.X:,.2f} USD/year")
    print(f"Operational Expenditure without the project (OPEX_0): {OPEX0.X:,.2f} USD/year")
    print(f"Operational Savings (Savings): {Savings.X:,.2f} USD/year")
    print(f"Energy expenses with the project (Eb): {Eb.X:,.2f} USD/year")
    print(f"Energy expenses without the project (Eb0): {Eb0.X:,.2f} USD/year")
    print(f"Energy earnings of the project (Es): {Es.X:,.2f} USD/year")
    print(f"Capacity Charges (CP): {CapacityP.X:,.2f} USD/year")
    print(f"Capacity Charges without the project (CP0): {CapacityP0.X:,.2f} USD/year")
    print(f"Capital Expenditure (CAPEX): {Investment0.X:,.2f} USD")
    print(f"BESS battery cost: {BESSbatteryCost:,.2f} USD")
    print(f"BESS inverter cost: {BESSinverterCost:,.2f} USD")
    print(f"PV inverter cost: {PVinverterCost:,.2f} USD")
    print(f"PV modules cost: {CAPEX_pv*Ppvinst.X:,.2f} USD")
    print(f"Soft Costs: {Investment0.X*(Sc-1)/Sc:,.2f} USD")   # CAPEX = Sc*(direct costs), so soft costs = CAPEX*(Sc-1)/Sc
    print(f"PV system cost (modules + PV inverter): {PVsystemCost:,.2f} USD")
    print(f"Net Present Value: {npv_var.X:,.2f} USD  (Objective)")
    print(f"Project Cash Flow: {CF:,.2f} USD/year")
    print(f"Internal Rate of Return (constant cash flow):      {TIR:,.2f} %")
    print(f"Internal Rate of Return (escalated cash flow):     {TIRe:,.2f} %")
    print(f"Simple Pay Back Time:    {NPERaprox:,.2f} years")
    print(f"Discounted Pay Back Time (constant cash flow):     {NPER:,.2f} years")
    print(f"Discounted Pay Back Time (escalated cash flow):    {NPERe:,.2f} years")
    print(f"------Levelized costs (USD/MWh)------------")
    print(f"Backup charge (Cresp): {Cresp_year:,.2f} USD/year")
    print(f"Asset LCOE 1 (CAPEX and O&M, load energy Wl):              {LCOEgross:,.2f} USD/MWh")
    print(f"Asset LCOE 2 (CAPEX and O&M, PV energy used wpv):          {LCOEgross2:,.2f} USD/MWh")
    print(f"Asset LCOE 3 (CAPEX and O&M, PV energy generated wpvmx):   {LCOEgross3:,.2f} USD/MWh")
    print(f"Supply LCOE (asset LCOE 1 + net grid cost, Wl):            {LCOEsupply:,.2f} USD/MWh")
    print(f"Grid-only LCOE 0 (no project, Wl):                         {LCOE0:,.2f} USD/MWh")
    print(f"Net LCOE 1 (= -NPV*crf/Wl = Supply - LCOE 0):              {LCOEnet1:,.2f} USD/MWh")
    print(f"Net LCOE 2 (= Supply LCOE):                                {LCOEnet2:,.2f} USD/MWh")
    print(f"Unserved load (Wu): {Wu.X:,.2f} kWh/year, served {100*(1-Wu.X/Wl.X):,.2f} % of the load energy")
    print(f"Asset LCOE over served energy (Ember definition):          {1000*(Io + OPEXgross/crfe)*crf/(Wl.X-Wu.X):,.2f} USD/MWh")
    print(f"Benefit-Cost Ratio (constant cash flow) {BCratio:,.2f}")
    print(f"Benefit-Cost Ratio (escalated cash flow) {BCratioe:,.2f}")  
    # Summary of all results, for tables
    summary = [
        ("Case", CASE, "-"),
        ("BESS capacity C", C.X, "kWh"), ("BESS inverter", PinverterBESS.X, "kW"),
        ("PV capacity", Ppvinst.X, "kWp"), ("PV inverter", PinverterPV.X, "kW"), ("ILR", ILR, "-"),
        ("Load energy Wl", Wl.X, "kWh/yr"), ("Energy bought Wb", Wb.X, "kWh/yr"), ("Energy sold Ws", Ws.X, "kWh/yr"),
        ("BESS charged Wc", Wc.X, "kWh/yr"), ("BESS discharged Wd", Wd.X, "kWh/yr"),
        ("PV generated wpvmx", wpvmx.X, "kWh/yr"), ("PV used wpv", wpv.X, "kWh/yr"), ("PV clipping", wclipping.X, "kWh/yr"),
        ("CAPEX I0", Io, "USD"), ("BESS battery cost", BESSbatteryCost, "USD"), ("BESS inverter cost", BESSinverterCost, "USD"),
        ("PV inverter cost", PVinverterCost, "USD"), ("PV modules cost", CAPEX_pv*Ppvinst.X, "USD"),
        ("Soft costs", Io*(Sc-1)/Sc, "USD"),
        ("O&M of assets", OPEXgross, "USD/yr"), ("Energy purchases Eb", Eb.X, "USD/yr"), ("Energy sales Es", Es.X, "USD/yr"),
        ("Capacity charges CP", CapacityP.X, "USD/yr"), ("Backup charge Cresp", Cresp_year, "USD/yr"),
        ("OPEX", OPEX.X, "USD/yr"), ("Purchases without project Eb0", Eb0.X, "USD/yr"),
        ("Capacity charges without project CP0", CapacityP0.X, "USD/yr"), ("OPEX0", OPEX0.X, "USD/yr"),
        ("Savings", Savings.X, "USD/yr"), ("Cash flow CF", CF, "USD/yr"),
        ("NPV", npv_var.X, "USD"), ("IRR constant CF", TIR, "%"), ("IRR escalated CF", TIRe, "%"),
        ("Simple payback", NPERaprox, "yr"), ("Discounted payback constant CF", NPER, "yr"),
        ("Discounted payback escalated CF", NPERe, "yr"),
        ("Benefit-cost ratio constant CF", BCratio, "-"), ("Benefit-cost ratio escalated CF", BCratioe, "-"),
        ("Asset LCOE 1 (Wl)", LCOEgross, "USD/MWh"), ("Asset LCOE 2 (wpv)", LCOEgross2, "USD/MWh"),
        ("Asset LCOE 3 (wpvmx)", LCOEgross3, "USD/MWh"), ("Supply LCOE (Wl)", LCOEsupply, "USD/MWh"),
        ("Grid-only LCOE 0 (Wl)", LCOE0, "USD/MWh"), ("Unserved load Wu", Wu.X, "kWh/yr"),
        ("Asset LCOE over served energy", 1000*(Io + OPEXgross/crfe)*crf/(Wl.X-Wu.X), "USD/MWh"), ("Net LCOE 1 (Wl)", LCOEnet1, "USD/MWh"),
        ("Net LCOE 2 (Wl)", LCOEnet2, "USD/MWh"),
        ("crf", crf, "-"), ("crfe", crfe, "-"),
    ]
    with open(f'results_summary_{CASE}.csv', 'w', newline='') as fsum:
        wsum = csv.writer(fsum); wsum.writerow(["Variable", "Value", "Unit"]); wsum.writerows(summary)
    print(f"Summary written to results_summary_{CASE}.csv")
else:
    print("N optimal solution, state code:", m.Status)
    
if m.status == gp.GRB.OPTIMAL:
    with open('solution_sizing_model.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write the header
        writer.writerow(["Hour", "Pb", "Ppv", "Pd", "Ppvmax", "Ps", "Pc", "SOC",	
"SpotPrice","npv" ])
        # Write each hour's data
        for t in T:
            writer.writerow([
                t,
                Pb[t].X,
                Ppv[t].X,
                Pd[t].X,
                Ppvmx[t].X,
                Ps[t].X,
                Pc[t].X,
                SOC[t].X,
                data[t]['lambda'],
                npv_var.X
            ])
    print("Solution written to solution_sizing_model.csv") 
    
    with open('solution_daily_summary.csv', 'w', newline='') as csvfile_daily:
        writer_daily = csv.writer(csvfile_daily)
        # Encabezados (Suma de energía diaria para potencias, Promedio para SOC)
        writer_daily.writerow(["Day", "Pb_Sum", "Ppv_Sum", "Pd_Sum", "Ppvmax_Sum", "Ps_Sum", "Pc_Sum", "SOC_Avg", "Benefit", "npv"])
        
        # Iterar sobre los 365 días
        for d in range(365):
            day_hours = range(d * 24 + 1, (d + 1) * 24 + 1)
            
            # Calculamos las sumas para las variables de flujo/potencia
            pb_sum = sum(Pb[t].X for t in day_hours)
            ppv_sum = sum(Ppv[t].X for t in day_hours)
            pd_sum = sum(Pd[t].X for t in day_hours)
            ppvmax_sum = sum(Ppvmx[t].X for t in day_hours)
            ps_sum = sum(Ps[t].X for t in day_hours)
            pc_sum = sum(Pc[t].X for t in day_hours)
            
            # Calculamos el promedio para el Estado de Carga (SOC)
            soc_avg = sum(SOC[t].X for t in day_hours) / 24
            
            writer_daily.writerow([
                d + 1,      # Día
                pb_sum,     # kWh/día comprados
                ppv_sum,    # kWh/día PV usados
                pd_sum,     # kWh/día descargados BESS
                ppvmax_sum, # kWh/día PV disponibles
                ps_sum,     # kWh/día vendidos
                pc_sum,     # kWh/día cargados BESS
                soc_avg,    # % promedio o kWh promedio en BESS
                Benefit.X,  # Valor constante del proyecto
                npv_var.X   # Valor constante del proyecto
            ])
    print("Daily summary written to solution_daily_summary.csv")