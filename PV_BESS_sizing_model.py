#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat March 20 2026
@authors: Paulo De Oliveira (PDEOLIV@GMAIL.COM) y Alejandro Salas Durán
Optimal Sizing + BESS + PV  ─ Case 8760h
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
       'lambda'  : 'lambda.inc', # 8760h 2024 Spanish Spot prices Eur/MWh
       'psi'     : 'psi.inc', # Tariff of Use of the network Eur/MWh
       'Ppvu'    : 'PpvuMadridSarah20052023.inc',  # Perfil Solar Madrid España Sarah 2005-2023
       "periodo" : "periodo.inc",  # calendar to assingn type of contracted power 1..6 per each hour
# ----------------------
}
series = {name: read_inc(route) for name, route in paths.items()}
T = range(1, 8761)
data = {t: {'lambda': series['lambda'].get(t, 0.0), 'Plu': series['Plu'].get(t, 0.0),
            'psi': series['psi'].get(t, 0.0), 'Ppvu': series['Ppvu'].get(t, 0.0)} for t in T}
periodo = {t: int(series["periodo"].get(t, 6)) for t in T}
Plinst  = 1000.0                    # kW peak/installed Load DataCenter
#Plinst  = 564.3                    # kW peak/installed load Industry 1
#Plinst  = 630.620646500000         # kW peak/installed load Industry 2
Rmax   = 1000                       # W/m2
Area   = 2.4                        # Area module 500W
eta    = 0.2094                     # PV module efficiency
eff_c  = 0.9624                     # BESS charge inverter efficiency
eff_d  = 0.9624                     # BESS discharge inverter efficiency
eff_pv = 0.9624                     # PV inverter efficiency
DoD    = 0.90                       # Depth of Discharge
PmaxF   = Plinst                    # límit frontier export/import
er = 1.1 # (exchange rate 1 Eur = 1.1 USD, average 2024)
BoP = 0  # Balance of Plant, USD 
Sc=1.2 # Soft costs, %
OaMpv = 12.5  #USD/kWp
OaMbess = 5.9    #USD/kWp
CAPEX_pv = 436 #USD/kWp
CAPEX_BESS = 185 #USD/kWh
CAPEX_BESS_inverter = 48 #USD/kWp
i=7.7/100 #annual discount rate (Lazard)
n=20 #project life span
e=2.5/100 # escalation rate
ir=(i-e)/(1+e) #equivalent annual discount rate including escalation
crf=(i*(i+1)**n)/((i+1)**n-1) #capital recovery factor
crfe=(1+e)*(ir*(ir+1)**n)/((ir+1)**n-1) #capital recovery factor modified
# ── Capacity charges in Spain (USD/kW-year)
kappa = [0, 28.79187*er, 15.07764*er, 6.55917*er, 5.17209*er, 1.93281*er, 0.91609*er] 
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
wcurtail = m.addVar(name='wcurtail', lb=-GRB.INFINITY); Wb = m.addVar(name='Wb', lb=-GRB.INFINITY)
Ws = m.addVar(name='Ws', lb=-GRB.INFINITY); Wl = m.addVar(name='Wl', lb=-GRB.INFINITY)
Wd = m.addVar(name='Wd', lb=-GRB.INFINITY); Wc = m.addVar(name='Wc', lb=-GRB.INFINITY)
CashFlow_var = m.addVar(name='CashFlow', lb=-GRB.INFINITY); Investment0 = m.addVar(name='Investment', lb=0)
PinverterBESS = m.addVar(name='PinverterBESS', lb=0); Ppvinst = m.addVar(name='Ppvinst', lb=0)
SOC0 = m.addVar(name='SOC0', lb=0); nx = m.addVar(name='nx', lb=0); C = m.addVar(name='C', lb=0)
Pbmax = {p: m.addVar(lb=0, name=f'PbmaxP{p}') for p in range(1, 7)}
SOC = m.addVars(T, lb=0, name='SOC'); Ppv = m.addVars(T, lb=0, name='Ppv')
Ppvmx = m.addVars(T, lb=0, name='Ppvmx'); Pd = m.addVars(T, lb=0, name='Pd')
Pc = m.addVars(T, lb=0, name='Pc'); Pb = m.addVars(T, lb=0, name='Pb'); Ps = m.addVars(T, lb=0, name='Ps')
w1 = m.addVars(T, vtype=GRB.BINARY, name='w1'); w3 = m.addVars(T, vtype=GRB.BINARY, name='w3')
# Constraints
m.addConstrs((Pbmax[periodo[t]] >= Pb[t] for t in T), "res_Pbmax")
m.addConstr(Pbmax[6] == PmaxF)
for p in range(1, 6): m.addConstr(Pbmax[p+1] >= Pbmax[p])
for t in T:
    m.addConstr(Pd[t] + Pb[t] + Ppv[t] == Pc[t] + Ps[t] + Plinst * data[t]['Plu'])
    m.addConstr(Ppvmx[t] == Ppvinst * eff_pv * data[t]['Ppvu']) 
    if t == 1:
        m.addConstr(SOC[t] == SOC0 + Pc[t]*eff_c - Pd[t]/eff_d)
    else:
        m.addConstr(SOC[t] == SOC[t-1] + Pc[t]*eff_c - Pd[t]/eff_d)
    m.addConstr(Ppv[t] <= Ppvmx[t])
    m.addConstr(Pc[t] <= PinverterBESS * w1[t])
    m.addConstr(Pd[t] <= PinverterBESS * (1-w1[t]))
    m.addConstr(Pb[t] <= PmaxF * w3[t])
    m.addConstr(Ps[t] <= PmaxF * (1-w3[t]))
    m.addConstr(SOC[t] <= ((1-DoD)/2 + DoD)*C)
    m.addConstr(SOC[t] >= ((1-DoD)/2)*C)
# FINANCIAL EQUATIONS
m.addConstr(Es == er*(quicksum(data[t]['lambda'] * Ps[t] for t in T)))
m.addConstr(Eb == er*(quicksum((data[t]['lambda'] + data[t]['psi']) * Pb[t] for t in T)))
m.addConstr(Eb0 == er*(quicksum((data[t]['lambda'] + data[t]['psi']) * Plinst * data[t]['Plu'] for t in T)))
m.addConstr(CapacityP == quicksum(kappa[p]*Pbmax[p] for p in range(1, 7)))
m.addConstr(CapacityP0 == sum(kappa[p]*PmaxF for p in range(1, 7)))
m.addConstr(OPEX0 == Eb0 + CapacityP0)
m.addConstr(OPEX == CapacityP + Eb + OaMpv*Ppvinst + OaMbess*C)
m.addConstr(Savings == OPEX0 - OPEX)
m.addConstr(Benefit == Es - OPEX)
m.addConstr(CashFlow_var == Es + OPEX0 - OPEX)
m.addConstr(Investment0 == BoP + Sc * (CAPEX_pv*Ppvinst + CAPEX_BESS*C + CAPEX_BESS_inverter*PinverterBESS))
m.addConstr(npv_var == CashFlow_var/(crfe) - Investment0)
# ENERGY DISPATCH
m.addConstr(wpv == quicksum(Ppv[t] for t in T))
m.addConstr(wpvmx == (1/eff_pv)*quicksum(Ppvmx[t] for t in T))
m.addConstr(Wb == quicksum(Pb[t] for t in T))
m.addConstr(Ws == quicksum(Ps[t] for t in T))
m.addConstr(Wc == quicksum(Pc[t] for t in T))
m.addConstr(Wd == quicksum(Pd[t] for t in T))
m.addConstr(Wl == sum(Plinst * data[t]['Plu'] for t in T))
m.addConstr(wcurtail == wpvmx/eff_pv - wpv)
m.addConstr(PinverterBESS <= C*2.0)
m.addConstr(PinverterBESS >= C*0.1)
m.addConstr(nx == 1000*Ppvinst/(Rmax*eta*Area))
m.setObjective(npv_var, GRB.MAXIMIZE)
#m.addConstr(Wb <= 0.12*Wl) 
#m.addConstr(Ws == 0.00*Wl)
# m.addConstr(C == 17000, name='r10')
# m.addConstr(PinverterBESS == 2, name='r11')
# m.addConstr(Ppvinst    == 6000, name='r12')
# m.addConstr(PbmaxP1   == PmaxF, name='r13')
# m.addConstr(PbmaxP2   == PmaxF, name='r14')
# m.addConstr(PbmaxP3   == PmaxF, name='r15')
# m.addConstr(PbmaxP4   == PmaxF, name='r16')
# m.addConstr(PbmaxP5   == PmaxF, name='r17')
# m.addConstr(PbmaxP6   == PmaxF, name='r18')  
m.optimize()
Crate = PinverterBESS.X/C.X #resulting BESS Crate
BESSinverterCost = CAPEX_BESS_inverter*PinverterBESS.X
BESSbatteryCost = CAPEX_BESS*C.X
PVsystemCost =CAPEX_pv*Ppvinst.X
if m.Status == GRB.OPTIMAL:
    Io = Investment0.X
    CF = CashFlow_var.X
    TIR = npf.irr([-Io] + [CF]*n) * 100 if Io > 0 else 0
    BCratio = npf.pv(rate=i, nper=n, pmt=-CF, fv=0) / Io if Io > 0 else 0
    OPEXgross = OaMpv*Ppvinst.X + OaMbess*C.X
    LCOEgross = 1000*(Io + OPEXgross/crfe)/(Wl.X/crf) if Wl.X > 0 else 0
    LCOEnet = 1000*(Io + (OPEX.X - OPEX0.X - Es.X)/crfe)/(Wl.X/crf) if Wl.X > 0 else 0
    NPERaprox = Io/CF if CF > 0 else 0
    try:
        NPER = np.log((CF - i * 0) / (CF + i * (-Io))) / np.log(1 + i) if (CF + i * (-Io)) > 0 else 0
    except: NPER = 0
  
# ─────────────────────────────────────────────────────────────────────────────
if m.Status == GRB.OPTIMAL:
    print(f"Results:")
    print(f"BESS capacity (WbessInst): {C.X:,.2f} kWh")
    print(f"BESS inverter capacity (PbessInst): {PinverterBESS.X:,.2f} kW")
    print(f"PV System capacity (PpvInst): {Ppvinst.X:,.2f} kW")
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
    print(f"PV energy injected (wpv): {wpv.X:,.2f} kWh/year")
    print(f"PV energy curtailed (Wcurtail): {wcurtail.X:,.2f} kWh/year")
    print(f"Initial SOC:  {SOC0.X:,.2f} MWh")
    print(f"BESS C-rate:  {Crate:,.2f} 1/s")
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
    print(f"Soft Costs: {Investment0.X*(Sc-1):,.2f} USD")
    print(f"PV system cost: {PVsystemCost:,.2f} USD")
    print(f"Net Present Value: {npv_var.X:,.2f} USD")
    print(f"Project Cash Flow: {CF:,.2f} USD/year")
    print(f"Internal Rate of Return:      {TIR:,.2f} %")
    print(f"Pay Back Time:     {NPER:,.2f} years")
    print(f"Simple Pay Back Time:    {NPERaprox:,.2f} years")
    print(f"Net LCOE (include earnings and savings):     {LCOEnet:,.2f} USD/MWh")
    print(f"Gross LCOE: (only CAPEX and OPEX)    {LCOEgross:,.2f} USD/MWh")
    print(f"Benefit-Cost Ratio {BCratio:,.2f}")  
else:
    print("N optimal solution, state code:", m.Status)
    
if m.status == gp.GRB.OPTIMAL:
    with open('solution_sizing_model.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write the header
        writer.writerow(["Hour", "Pb", "Ps", "Pc", "Pd", "SOC",	
"Benefit","npv" ])
        # Write each hour's data
        for t in T:
            writer.writerow([
                t,
                Pb[t].X,
                Ps[t].X,
                Pc[t].X,
                Pd[t].X,
                SOC[t].X,
                Benefit.X,
                npv_var.X
            ])
    print("Solution written to solution_sizing_model.csv")    