"""
Model for RNG/NG
"""

import pulp
import pandas as pd

input_df = pd.DataFrame() # import from Julie's

# Time-series constants
SBG = list(input_df['SBG'])
D = list(input_df['NG_demand'])
HOEP = list(input_df['HOEP'])

# Fixed constants
nu_electrolyzer = 3.55 / 4.63
E_HHV_H2 = 3.55
nu_reactor = 213.6 / 4 / 68.7
HHV_H2 = 68.7
HHV_NG = 10.55
CO2_available = 1300.26
E_electrolyzer_min = 0
E_electrolyzer_max = 1000
tau = 0.50
beta = 1.35
C_0 = 1324.3
mu = 0.707
gamma = 1714.8
k = 2e6
C_upgrading = 1172.3
C_CO2 = 0.172
TC = 0.008804
C_H2O = 0.00314
WCR = 0.4
OPEX_upgrading = 146.5

# MILP model
LP = pulp.LpProblem('LP', pulp.LpMinimize)
RNG_max = pulp.LpVariable('RNG_max',
                          lowBound=0,
                          cat='Continuous')
N_electrolyzer = pulp.LpVariable('N_electrolyzer',
                                 lowBound=0
                                 cat='Continuous')
alpha = pulp.LpVariable.dicts('alpha',
                          list(range(1, 31)),
                          cat='Binary')
E = pulp.LpVariable.dicts('E_h',
                          list(input_df.index),
                          lowBound=0,
                          cat='Continuous')
H2 = pulp.LpVariable.dicts('H2_h',
                          list(input_df.index),
                          lowBound=0,
                          cat='Continuous')
RNG = pulp.LpVariable.dicts('RNG_h',
                          list(input_df.index),
                          lowBound=0,
                          cat='Continuous')
CO2 = pulp.LpVariable.dicts('CO2_h',
                          list(input_df.index),
                          lowBound=0,
                          cat='Continuous')
NG = pulp.LpVariable.dicts('NG_h',
                          list(input_df.index),
                          lowBound=0,
                          cat='Continuous')
CAPEX = pulp.LpVariable('CAPEX', cat='Continuous')

for h in input_df.index:
    # Energy and flow constraints
    if h == 0:
        LP += RNG_max <= CO2_available
    LP += H2[h] == nu_electrolyzer * E[h] / E_HHV_H2
    LP += RNG[h] == nu_reactor * H2[h] * HHV_H2 / HHV_NG
    LP += CO2[h] == RNG[h]

    # Demand constraint
    LP += NG[h] == D[h] - RNG[h]

    # Electrolyzer and reactor constraints
    LP += N_electrolyzer * E_electrolyzer_min <= E[h]
    LP += N_electrolyzer * E_electrolyzer_max >= E[h]
    LP += E[h] <= SBG[h]
    LP += RNG[h] <= RNG_max
    if h > 0:
        LP += -RNG_max * tau <= RNG[h] - RNG[h - 1]
        LP += RNG_max * tau >= RNG[h] - RNG[h - 1]
    if h == 0:
        LP += pulp.LpAffineExpression(dict(zip(alpha, range(1, 31)))) == N_electrolyzer
        LP += pulp.lpSum(alpha) <= 1

# CAPEX
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, 31)]
LP += pulp.LpAffineExpression(dict(zip(alpha, C_electrolyzer))) + \
      gamma * RNG_max + k + C_upgrading * RNG_max == CAPEX

# OPEX
LP += pulp.LpAffineExpression(dict(zip(CO2, [C_CO2] * len(input_df.index)))) + \
      pulp.LpAffineExpression(dict(zip(E, [x + TC for x in HOEP]))) + \
      pulp.LpAffineExpression(dict(zip(H2, [C_H2O * WCR] * len(input_df.index)))) + \
      OPEX_upgrading * RNG_max == OPEX

# Objective
LP += CAPEX + OPEX, 'Z'

LP.solve()
print([x.varValue for x in LP.variables()])
print(pulp.value(LP.objective))
