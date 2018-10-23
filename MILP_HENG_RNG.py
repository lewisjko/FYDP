"""
Model for RNG/NG + HENG
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
nu_reactor = var['value']['meth_reactor_eff']
HHV_H2 = var['value']['HHV_H2']
HHV_NG = var['value']['HHV_NG']
CO2_available = var['value']['CO2_available']
E_electrolyzer_min = var['value']['min_E_cap']
E_electrolyzer_max = var['value']['max_E_cap']
tau = 0.50

EMF_NG = var['value']['EMF_NG']
EMF_comb = var['value']['EMF_combRNG']
EMF_nuc = var['value']['EMF_nuclear']
EMF_bio = var['value']['EMF_bioCO2']
EMF_electrolyzer = var['value']['EMF_electrolyzer']
EMF_reactor = var['value']['EMF_reactor']

beta = var['value']['beta']
C_0 = var['value']['C_0']
mu = var['value']['mu']
gamma = var['value']['gamma']
k = var['value']['k']
C_upgrading = var['value']['C_upgrading']
C_CO2 = var['value']['C_CO2']
TC = var['value']['TC']
C_H2O = var['value']['C_H2O']
WCR = var['value']['water_cons_rate']
OPEX_upgrading = var['value']['OPEX_upgrading']
TVM = var['value']['TVM']

# RNG + HENG model
LP = pulp.LpProblem('LP', pulp.LpMinimize)

# RNG Variables
RNG_max = pulp.LpVariable('RNG_max',
                          lowBound=0,
                          cat='Continuous')
N_electrolyzer_1 = pulp.LpVariable('N_electrolyzer_1',
                                 lowBound=0,
                                 cat='Integer')
alpha_1 = pulp.LpVariable.dicts('alpha_1',
                          [str(i) for i in range(1, 31)],
                          cat='Binary')

E_1 = pulp.LpVariable.dicts('E_1',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
H2_1 = pulp.LpVariable.dicts('H2_1',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
RNG = pulp.LpVariable.dicts('RNG',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
CO2 = pulp.LpVariable.dicts('CO2_1',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

CAPEX_1 = pulp.LpVariable('CAPEX_1', lowBound=0, cat='Continuous')
OPEX_1 = pulp.LpVariable('OPEX_1', lowBound=0, cat='Continuous')

# HENG Variables
N_electrolyzer_2 = pulp.LpVariable('N_electrolyzer_2',
                                 lowBound=0,
                                 cat='Integer')
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
CAPEX_2 = pulp.LpVariable('CAPEX_2', cat='Continuous')
OPEX_2 = pulp.LpVariable('OPEX_2', cat='Continuous')

# Shared Variables
NG = pulp.LpVariable.dicts('NG',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
em_offset_max = pulp.LpVariable('em_offset_max',
                          lowBound=0,
                          cat='Continuous')
em_offset = pulp.LpVariable('em_offset',
                          lowBound=0,
                          cat='Continuous')
em_rng = pulp.LpVariable('em_rng',
                          lowBound=0,
                          cat='Continuous')
em_heng = pulp.LpVariable('em_heng',
                          lowBound=0,
                          cat='Continuous')
em_ng = pulp.LpVariable('em_ng',
                          lowBound=0,
                          cat='Continuous')

for i, h in enumerate([str(i) for i in input_df.index]):
    # Energy and flow constraints
    LP += H2_1[h] + H2_2[h] == nu_electrolyzer * (E_1[h] + E_2[h]) / E_HHV_H2

    if h == '0':
        LP += RNG_max <= CO2_available
    LP += RNG[h] == nu_reactor * H2_1[h] * HHV_H2 / HHV_NG
    LP += CO2[h] == RNG[h]

    # Demand constraint
    LP += D[i] == RNG[h] + NG[h] + H2_2[h]

    # Supply constraint
    LP += E_1[h] + E_2[h] <= SBG[i]

    # Electrolyzer constraints
    LP += N_electrolyzer_1 * E_electrolyzer_min <= E_1[h]
    LP += N_electrolyzer_1 * E_electrolyzer_max >= E_1[h]
    LP += N_electrolyzer_2 * E_electrolyzer_min <= E_2[h]
    LP += N_electrolyzer_2 * E_electrolyzer_max >= E_2[h]

    # Reactor constraints
    LP += RNG[h] <= RNG_max
    if h != '0':
        LP += -RNG_max * tau <= RNG[h] - RNG[str(i - 1)]
        LP += RNG_max * tau >= RNG[h] - RNG[str(i - 1)]

    # HENG-specific constraints
    LP += 0.95 * H2_2[h] <= 0.05 * NG[h]

    # Integer constraints
    if h == '0':
        LP += pulp.lpSum(n * alpha_1[str(n)] for n in range(1, 31)) == N_electrolyzer_1
        LP += pulp.lpSum(alpha_1) <= 1
        LP += pulp.lpSum(n * alpha_2[str(n)] for n in range(1, 31)) == N_electrolyzer_2
        LP += pulp.lpSum(alpha_2) <= 1

# Emission constraints
LP += pulp.lpSum(EMF_comb * RNG[h] + EMF_nuc * E_1[h] + EMF_bio * CO2[h] + \
                 EMF_electrolyzer * H2_1[h] + EMF_reactor * RNG[h] \
                 for h in [str(x) for x in input_df.index]) == em_rng
LP += pulp.lpSum(EMF_NG * NG[h] + EMF_nuc * E_2[h] + EMF_electrolyzer * H2_2[h] \
                   for h in [str(x) for x in input_df.index]) == em_heng
LP += pulp.lpSum(EMF_NG * D[h] for h in input_df.index) == em_ng
LP += em_offset == em_ng - em_rng - em_heng

# Electrolyzer cost
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, 31)]

# CAPEX RNG
LP += pulp.lpSum(alpha_1[str(n)] * C_electrolyzer[n - 1] for n in range(1, 31)) + \
      gamma * RNG_max + k + C_upgrading * RNG_max == CAPEX_1

# CAPEX HENG
LP += pulp.lpSum(alpha_2[str(n)] * C_electrolyzer[n - 1] for n in range(1, 31)) == CAPEX_2

# OPEX RNG
LP += pulp.lpSum(CO2[str(n)] * C_CO2 for n in input_df.index) + \
      pulp.lpSum(E_1[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
      pulp.lpSum(H2_1[str(n)] * C_H2O * WCR for n in input_df.index) + \
      OPEX_upgrading * RNG_max == OPEX_1

# OPEX HENG
LP += pulp.lpSum(E_2[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
      pulp.lpSum(H2_2[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_2
