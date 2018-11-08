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

# HENG model
LP_2 = pulp.LpProblem('LP', pulp.LpMinimize)
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

em_offset_2 = pulp.LpVariable('em_offset_max_1',
                          lowBound=0,
                          cat='Continuous')
em_heng = pulp.LpVariable('em_heng',
                          lowBound=0,
                          cat='Continuous')
em_ng = pulp.LpVariable('em_ng',
                          lowBound=0,
                          cat='Continuous')

CAPEX_2 = pulp.LpVariable('CAPEX_2', cat='Continuous')
OPEX_2 = pulp.LpVariable('OPEX_2', cat='Continuous')

for LP in [LP_2, LP_eps]:
    for i, h in enumerate([str(i) for i in input_df.index]):
        # Energy and flow constraints
        LP += H2_2[h] == nu_electrolyzer * E_2[h] / E_HHV_H2

        # Demand constraint
        LP += NG_2[h] + H2_2[h] == D[i]

        # NG and Electrolyzer constraints
        LP += 0.95 * H2_2[h] <= 0.05 * NG_2[h]
        LP += N_electrolyzer_2 * E_electrolyzer_min <= E_2[h]
        LP += N_electrolyzer_2 * E_electrolyzer_max >= E_2[h]
        LP += E_2[h] <= SBG[i]
        if h == '0':
            LP += pulp.lpSum(n * alpha_2[str(n)] for n in range(1, N_max)) == N_electrolyzer_2
            LP += pulp.lpSum(alpha_2) <= 1

    # Emission constraints
    LP += pulp.lpSum(EMF_NG * NG_2[h] + EMF_nuc * E_2[h] + EMF_electrolyzer * H2_2[h] \
                       for h in [str(x) for x in input_df.index]) == em_heng
    LP += pulp.lpSum(EMF_NG * D[h] for h in input_df.index) == em_ng

# Eps Objective
LP_eps += em_ng - em_heng, 'max offset'
LP_eps.solve()
offset_max = LP_eps.objective.value()

# CAPEX
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, N_max)]
LP_2 += pulp.lpSum(alpha_2[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_max)) == CAPEX_2

# OPEX
LP_2 += pulp.lpSum(E_2[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
      pulp.lpSum(H2_2[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_2

# Cost Objective
phi = 0.80
LP_2 += em_ng - em_heng == em_offset_2
LP_2 += em_offset_2 >= phi * offset_max
LP_2 += CAPEX_2 + OPEX_2, 'Cost_2'

LP_2.solve()
