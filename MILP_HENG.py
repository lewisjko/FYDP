"""
Model for HENG
"""

import pulp

# Time-series constants
SBG = list(input_df['SBG(kWh)'])
D = list(input_df['NG_demand(m^3)'])
HOEP = list(input_df['HOEP'])
EMF = list(input_df['EMF(tonne/kWh)'])

# Fixed constants
N_max = 30000
nu_electrolyzer = var['value']['electrolyzer_eff']
E_HHV_H2 = var['value']['E_hhv_h2']
E_electrolyzer_min = var['value']['min_E_cap']
E_electrolyzer_max = var['value']['max_E_cap']

EMF_NG = var['value']['EMF_NG']
EMF_nuc = var['value']['EMF_nuclear']
EMF_electrolyzer = var['value']['EMF_electrolyzer']

beta = var['value']['beta']
C_0 = var['value']['C_0']
mu = var['value']['mu']
gamma = var['value']['gamma']
k = var['value']['k']
TC = var['value']['TC']
C_H2O = var['value']['C_H2O']
WCR = var['value']['water_cons_rate']
TVM = var['value']['TVM']

# Tank and compressor constants
I_max = var['value']['Imax'] # kmol
I_min= var['value']['Imin'] # kmol
F_max_prestorage =var['value']['Fmax_prestorage'] # kmol

CAPEX_prestorage = var['value']['CAPEX_prestorage'] # $
CAPEX_tank = var['value']['CAPEX_tank'] # $

ECF_prestorage = var['value']['ECF_prestorage'] # kWh/kmol H2

# HENG model
LP_cost = pulp.LpProblem('LP_cost', pulp.LpMinimize)
LP_eps = pulp.LpProblem('LP_eps', pulp.LpMaximize)
N_electrolyzer_2 = pulp.LpVariable('N_electrolyzer_2',
                          lowBound=0,
                          cat='Continuous')
alpha_2 = pulp.LpVariable.dicts('alpha_2',
                          [str(i) for i in range(1, N_max)],
                          cat='Binary')
E_2 = pulp.LpVariable.dicts('E_2',
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

em_offset_2 = pulp.LpVariable('em_offset_max_1',
                          lowBound=0,
                          cat='Continuous')
em_heng = pulp.LpVariable('em_heng',
                          lowBound=0,
                          cat='Continuous')
em_ng = pulp.LpVariable('em_ng',
                          lowBound=0,
                          cat='Continuous')

# Number of prestorage compressors
N_prestorage = pulp.LpVariable('N_prestorage',
                          lowBound=0,
                          cat='Integer')
# Number of tanks
N_tank = pulp.LpVariable('N_tank',
                          lowBound=0,
                          cat='Integer')
# H2 flow directly from electrolyzers (m^3/h)
H2_direct = pulp.LpVariable.dicts('H2_direct',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
# H2 flow going into tanks from electrolyzers (m^3/h)
H2_tank_in = pulp.LpVariable.dicts('H2_tank_in',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
# H2 flow coming out from tanks (m^3/h)
H2_tank_out = pulp.LpVariable.dicts('H2_tank_out',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
# Inventory level of H2 (m^3)
I_H2 = pulp.LpVariable.dicts('I_H2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

CAPEX_2 = pulp.LpVariable('CAPEX_2', cat='Continuous')
OPEX_2 = pulp.LpVariable('OPEX_2', cat='Continuous')

for LP in [LP_eps, LP_cost]:
    for i, h in enumerate([str(i) for i in input_df.index]):
        # Energy and flow constraints
        LP += H2_direct[h] + H2_tank_in[h] == nu_electrolyzer * E_2[h] * E_HHV_H2 ** -1
        
        # Hydrogen storage tank constraint
        if h == '0':
            LP += I_H2[h] == I_min * N_tank + H2_tank_in[h] - H2_tank_out[h]
        else:
            LP += I_H2[h] == I_H2[str(i - 1)] + H2_tank_in[h] - H2_tank_out[h]
        LP += I_H2[h] <= I_max * N_tank
        LP += I_H2[h] >= I_min * N_tank        

        # Demand constraint
        LP += NG_2[h] + H2_direct[h] + H2_tank_out[h] == D[i]

        # NG and Electrolyzer constraints
        LP += 0.95 * (H2_direct[h] + H2_tank_out[h]) <= 0.05 * NG_2[h]
        LP += N_electrolyzer_2 * E_electrolyzer_min <= E_2[h]
        LP += N_electrolyzer_2 * E_electrolyzer_max >= E_2[h]
        LP += E_2[h] <= SBG[i]
        if h == '0':
            LP += pulp.lpSum(n * alpha_2[str(n)] for n in range(1, N_max)) == N_electrolyzer_2
            LP += pulp.lpSum(alpha_2) <= 1

    # Emission constraints
    LP += pulp.lpSum(EMF_NG * NG_2[h] + EMF[int(h)] * E_2[h] + EMF_electrolyzer * (H2_direct[h] + H2_tank_in[h]) \
                     for h in [str(x) for x in input_df.index]) + \
          pulp.lpSum(EMF[int(h)] * ECF_prestorage * H2_tank_in[h] for h in [str(x) for x in input_df.index]) == em_heng
    LP += pulp.lpSum(EMF_NG * D[h] for h in input_df.index) == em_ng

# Eps Objective
LP_eps += em_ng - em_heng, 'max offset'
LP_eps.solve()
offset_max = LP_eps.objective.value()
print(LP_eps.status)

# CAPEX
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, N_max)]
LP_cost += pulp.lpSum(alpha_2[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_max)) + \
           (N_tank * CAPEX_tank + N_prestorage * CAPEX_prestorage) * 20 == CAPEX_2

# OPEX
LP_cost += pulp.lpSum(E_2[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
           pulp.lpSum((H2_direct[str(n)] + H2_tank_in[str(n)]) * C_H2O * WCR for n in input_df.index) + \
           pulp.lpSum(ECF_prestorage * H2_tank_in[str(n)] for n in input_df.index) == OPEX_2

# Objectives
phi = 0.80
LP_cost += em_ng - em_heng == em_offset_2
LP_cost += em_offset_2 >= phi * offset_max
LP_cost += CAPEX_2 + OPEX_2, 'Cost_2'

LP_cost.solve()
print(LP_cost.status)
