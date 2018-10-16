"""
Model for HENG
"""

import pulp
import pandas as pd

# Time-series constants
SBG = list(input_df['SBG(kWh)'])
D = list(input_df['NG_demand(m^3)'])
HOEP = list(input_df['HOEP'])

# Fixed constants
nu_electrolyzer = var['value']['electrolyzer_eff']
E_HHV_H2 = var['value']['E_hhv_h2']
E_electrolyzer_min = var['value']['min_E_cap']
E_electrolyzer_max = var['value']['max_E_cap']
beta = var['value']['beta']
C_0 = var['value']['C_0']
mu = var['value']['mu']
gamma = var['value']['gamma']
k = var['value']['k']
TC = var['value']['TC']
C_H2O = var['value']['C_H2O']
WCR = var['value']['water_cons_rate']

# MILP model
LP_2 = pulp.LpProblem('LP', pulp.LpMinimize)
N_electrolyzer_2 = pulp.LpVariable('N_electrolyzer_2',
                                 lowBound=0,
                                 cat='Continuous')
alpha_2 = pulp.LpVariable.dicts('alpha_2',
                          [str(i) for i in range(1, 31)],
                          cat='Binary')
E_2 = pulp.LpVariable.dicts('E_2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
H2_2 = pulp.LpVariable.dicts('H2_2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
CO2 = pulp.LpVariable.dicts('CO2_2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
NG_2 = pulp.LpVariable.dicts('NG_2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
CAPEX_2 = pulp.LpVariable('CAPEX_2', cat='Continuous')
OPEX_2 = pulp.LpVariable('OPEX_2', cat='Continuous')

for i, h in enumerate([str(i) for i in input_df.index]):
    # Energy and flow constraints
    LP_2 += H2_2[h] == nu_electrolyzer * E_2[h] / E_HHV_H2

    # Demand constraint
    LP_2 += NG_2[h] + H2_2[h] == D[i]

    # NG and Electrolyzer constraints
    LP_2 += 0.95 * H2_2[h] <= 0.05 * NG_2[h]
    LP_2 += N_electrolyzer_2 * E_electrolyzer_min <= E_2[h]
    LP_2 += N_electrolyzer_2 * E_electrolyzer_max >= E_2[h]
    LP_2 += E_2[h] <= SBG[i]
    if h == '0':
        LP_2 += pulp.lpSum(n * alpha_2[str(n)] for n in range(1, 31)) == N_electrolyzer_2
        LP_2 += pulp.lpSum(alpha_2) <= 1

# CAPEX
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, 31)]
LP_2 += pulp.lpSum(alpha_2[str(n)] * C_electrolyzer[n - 1] for n in range(1, 31)) == CAPEX_2

# OPEX
LP_2 += pulp.lpSum(E_2[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
      pulp.lpSum(H2_2[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_2

# Objective
LP_2 += CAPEX_2 + OPEX_2, 'Cost_2'

LP_2.solve()
