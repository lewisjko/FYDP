"""
Combined model (RNG+HENG+Mobility+Industry)
"""

import pulp
from numpy import count_nonzero

# Time-series constants
SBG = list(input_df['SBG(kWh)'])
NG_demand = list(input_df['NG_demand(m^3)'])
mobility_demand = list(input_df['mobility_demand(m^3)'])
industry_demand = list(input_df['industry_demand(m^3)']) # 5% of actual demand
HOEP = list(input_df['HOEP'])
EMF = list(input_df['EMF(tonne/kWh)'])

# Electrolyzer and flow constants
N_max = 30000
N_max += 1
nu_electrolyzer = var['value']['electrolyzer_eff']
E_HHV_H2 = var['value']['E_hhv_h2'] # kwh/m^3
nu_reactor = var['value']['meth_reactor_eff']
HHV_H2 = var['value']['HHV_H2'] # MMBtu/kmol
HHV_NG = var['value']['HHV_NG'] # MMBtu/kmol
CO2_available_total = var['value']['CO2_total'] # m^3/year
CO2_available = float(CO2_available_total) / count_nonzero(SBG)
E_electrolyzer_min = var['value']['min_E_cap'] # kwh
E_electrolyzer_max = var['value']['max_E_cap'] # kwh
tau = 0.50

# Emission constants
EMF_NG = var['value']['EMF_NG'] # tonne CO2/m^3 H2
EMF_comb = var['value']['EMF_combRNG'] # tonne CO2/m^3 RNG
EMF_bio = var['value']['EMF_bioCO2'] # tonne CO2/kWh
EMF_electrolyzer = var['value']['EMF_electrolyzer'] # tonne CO2/m^3 H2
EMF_reactor = var['value']['EMF_reactor'] # tonne CO2/m^3 RNG
EMF_vehicle = var['value']['emission_gasoline_v'] # tonne CO2/car/year
num_vehicle = var['value']['N_gasoline_v']
FCV_penetration = var['value']['FCV_penetration'] # FCV market penetration (estimated)
EMF_SMR = var['value']['EMF_SMR'] # kg CO2/kmol of H2

# Electrolyzer and flow cost
beta = var['value']['beta']
C_0 = var['value']['C_0'] # $/kW
mu = var['value']['mu']
gamma = var['value']['gamma']
k = var['value']['k'] # $
C_upgrading = var['value']['C_upgrading'] # $/m^3 reactor capacity
C_CO2 = var['value']['C_CO2'] # $/m^3 CO2
TC = var['value']['TC'] # $/kWh
C_H2O = var['value']['C_H2O'] # $/L
WCR = var['value']['water_cons_rate'] # L H2O/m^3 H2
OPEX_upgrading = var['value']['OPEX_upgrading'] # $/m^3 reactor capacity
TVM = var['value']['TVM']

# Tank and compressor constants
I_max = var['value']['Imax'] # kmol
I_min= var['value']['Imin'] # kmol
F_max_booster = var['value']['Fmax_booster'] # kmol
F_max_prestorage =var['value']['Fmax_prestorage'] # kmol

CAPEX_booster = var['value']['CAPEX_booster'] # $
CAPEX_prestorage = var['value']['CAPEX_prestorage'] # $
CAPEX_tank = var['value']['CAPEX_tank'] # $

ECF_prestorage = var['value']['ECF_prestorage'] # kWh/kmol H2

z_booster = var['value']['z_booster'] # compressibility factor for booster compressor
R = var['value']['R'] # kJ/kmolK
T = var['value']['T'] # K
comp_efficiency = var['value']['comp_efficiency']  # isentropic compressor efficiency
heat_cap_ratio = var['value']['heat_cap_ratio'] # heat capacaity ratio of hydrogen
P_in_booster = var['value']['P_in_booster'] # inlet pressure of booster compressor
P_out_booster = var['value']['P_out_booster'] # outlet pressure of booster compressor
N_stage_booster = var['value']['N_stage_booster']

ECF_booster = z_booster * R * T * N_stage_booster / comp_efficiency * \
              heat_cap_ratio / (heat_cap_ratio - 1) * \
              (((P_out_booster / P_in_booster) ** ((heat_cap_ratio - 1) / N_stage_booster / heat_cap_ratio ))-1) \
              / 3600 # converting kJ to kWh ECF booster in kWh/kmol

# Converting the tank and compressor constants to m^3
MW_H2 = var['value']['MW_H2'] # kg/kmol H2
density_H2 = var['value']['density_H2'] # kg/m^3

I_max = I_max * MW_H2 / density_H2 # m^3
I_min = I_min * MW_H2 / density_H2 # m^3
F_max_booster = F_max_booster * MW_H2 / density_H2 # m^3
F_max_prestorage = F_max_prestorage * MW_H2 / density_H2 # m^3

ECF_booster = ECF_booster / MW_H2 * density_H2 # kWh/m^3
ECF_prestorage = ECF_prestorage / MW_H2 * density_H2 # kWh/m^3
EMF_SMR = EMF_SMR * 0.001 / MW_H2 * density_H2 # tonne CO2/m^3 of H2

"""LP"""
# Maximize emission offset
LP_eps = pulp.LpProblem('LP_eps', pulp.LpMaximize)

# Minimize cost
LP_cost = pulp.LpProblem('LP_cost', pulp.LpMinimize)

"""RNG"""
# Max RNG flow (m^3/h)
RNG_max = pulp.LpVariable('RNG_max',
                          lowBound=0,
                          cat='Continuous')
# H2 flow to RNG (m^3/h)
H2_1 = pulp.LpVariable.dicts('H2_1',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
# RNG flow (m^3/h)
RNG = pulp.LpVariable.dicts('RNG',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
# CO2 flow (m^3/h)
CO2 = pulp.LpVariable.dicts('CO2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
# 0 if RNG_max = 0, 1 otherwise
alpha_RNG = pulp.LpVariable('alpha_RNG',
                          cat='Binary')

"""HENG"""
# H2 flow to HENG (m^3/h)
H2_2 = pulp.LpVariable.dicts('H2_2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

"""Mobility"""
# H2 flow to mobility fuel (m^3/h)
H2_3 = pulp.LpVariable.dicts('H2_3',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

"""Industry"""
# H2 flow to industry feed (m^3/h)
H2_4 = pulp.LpVariable.dicts('H2_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

"""Shared"""
# Number of booster compressors
N_booster = pulp.LpVariable('N_booster',
                          lowBound=0,
                          cat='Integer')
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
# Total SBG usage (kWh)
E = pulp.LpVariable.dicts('E',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
# Inventory level of H2 (m^3)
I_H2 = pulp.LpVariable.dicts('I_H2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
# Number of electrolyzers
N_electrolyzer = pulp.LpVariable('N_electrolyzer',
                          lowBound=0,
                          cat='Integer')
# alpha_n = 1 if N_electrolyzer = n, 0 otherwise
alpha = pulp.LpVariable.dicts('alpha',
                          [str(i) for i in range(1, N_max)],
                          cat='Binary')
# NG flow (m^3/h)
NG = pulp.LpVariable.dicts('NG',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

# Emission variables (tonnes CO2/yr)
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
em_electrolyzer = pulp.LpVariable('em_electrolyzer',
                          lowBound=0,
                          cat='Continuous')
em_booster_comp = pulp.LpVariable('em_booster_comp',
                          lowBound=0,
                          cat='Continuous')
em_pre_comp = pulp.LpVariable('em_pre_comp',
                          lowBound=0,
                          cat='Continuous')
em_smr = pulp.LpVariable('em_smr',
                          lowBound=0,
                          cat='Continuous')
em_sbg = pulp.LpVariable('em_sbg',
                          lowBound=0,
                          cat='Continuous')

# Electrolyzer CAPEX ($)
CAPEX_electrolyzer = pulp.LpVariable('CAPEX_Electrolyzer', lowBound=0, cat='Continuous')
# Reactor CAPEX ($)
CAPEX_reactor = pulp.LpVariable('CAPEX_reactor', lowBound=0, cat='Continuous')
# Tank + prestorage comp CAPEX ($)
CAPEX_storage = pulp.LpVariable('CAPEX_storage', lowBound=0, cat='Continuous')
# Booster compressor CAPEX ($)
CAPEX_booster_comp = pulp.LpVariable('CAPEX_booster_comp', lowBound=0, cat='Continuous')

# SBG treatment OPEX (electrolyzer + prestorage + tank) ($/yr)
OPEX_SBG = pulp.LpVariable('OPEX_SBG', lowBound=0, cat='Continuous')
# Reactor OPEX ($/yr)
OPEX_reactor = pulp.LpVariable('OPEX_reactor', lowBound=0, cat='Continuous')
# Booster compressor OPEX ($/yr)
OPEX_booster_comp = pulp.LpVariable('OPEX_booster_comp', lowBound=0, cat='Continuous')

# Total cost
total_cost = pulp.LpVariable('total_cost', lowBound=0, cat='Continuous')


for LP in [LP_eps, LP_cost]:
    for i, h in enumerate([str(i) for i in input_df.index]):
        # Energy and flow constraints
        LP += H2_direct[h] + H2_tank_in[h] == nu_electrolyzer * E[h] * E_HHV_H2 ** (-1)
        LP += H2_direct[h] + H2_tank_out[h] == H2_1[h] + H2_2[h] + H2_3[h] + H2_4[h]

        # Hydrogen storage tank constraint
        if h == '0':
            LP += I_H2[h] == I_min * N_tank + H2_tank_in[h] - H2_tank_out[h]
        else:
            LP += I_H2[h] == I_H2[str(i - 1)] + H2_tank_in[h] - H2_tank_out[h]
        LP += I_H2[h] <= I_max * N_tank
        LP += I_H2[h] >= I_min * N_tank

        # Prestorage compressor capacity constraint
        LP += H2_tank_in[h] <= N_prestorage * F_max_prestorage

        # Reactor constraints
        if h == '0':
            LP += RNG_max <= CO2_available
        LP += RNG[h] == nu_reactor * H2_1[h] * HHV_H2 / HHV_NG
        LP += CO2[h] == RNG[h]

        # NG demand constraint
        LP += NG_demand[i] == RNG[h] + NG[h] + H2_2[h]

        # Mobility demand constraint
        LP += H2_3[h] == mobility_demand[i]

        # Industry demand constraint
        LP += H2_4[h] == industry_demand[i]

        # Supply constraint
        LP += E[h] <= SBG[i]

        # Booster compressor capcity constraint (mobility)
        LP += H2_3[h] <= N_booster * F_max_booster

        # Electrolyzer constraints
        LP += N_electrolyzer * E_electrolyzer_min <= E[h]
        LP += N_electrolyzer * E_electrolyzer_max >= E[h]

        # Reactor constraints
        LP += RNG[h] <= RNG_max
        if h != '0':
            LP += -RNG_max * tau <= RNG[h] - RNG[str(i - 1)]
            LP += RNG_max * tau >= RNG[h] - RNG[str(i - 1)]

        # NG component constraints
        LP += 0.90 * RNG[h] <= 0.10 * (NG[h] + H2_2[h])
        LP += 0.95 * H2_2[h] <= 0.05 * (NG[h] + RNG[h])

    # Integer constraints
    LP += pulp.lpSum(n * alpha[str(n)] for n in range(1, N_max)) == N_electrolyzer
    LP += pulp.lpSum(alpha) <= 1

    # Emission constraints
    LP += pulp.lpSum(EMF_comb * RNG[h] + EMF_bio * CO2[h] + EMF_reactor * RNG[h] \
                     for h in [str(x) for x in input_df.index]) == em_rng
    LP += pulp.lpSum(EMF_NG * NG[h] for h in [str(x) for x in input_df.index]) == em_heng
    LP += pulp.lpSum(EMF_NG * NG_demand[h] for h in input_df.index) == em_ng
    em_gas_vehicle = 100000 * EMF_vehicle
    LP += pulp.lpSum(EMF_SMR * industry_demand[h] for h in input_df.index) == em_smr
    LP += pulp.lpSum(EMF[int(h)] * (E[h]) for h in [str(x) for x in input_df.index]) == em_sbg
    LP += pulp.lpSum(EMF_electrolyzer * (H2_direct[h] + H2_tank_in[h]) \
                     for h in [str(x) for x in input_df.index]) == em_electrolyzer
    LP += pulp.lpSum(EMF[n] * ECF_booster * H2_3[h] for n in input_df.index) == em_booster_comp
    LP += pulp.lpSum(EMF[n] * ECF_prestorage * H2_tank_in[str(n)] for n in input_df.index) == em_pre_comp

"""
Emission objective model
"""

LP_eps += em_ng + em_gas_vehicle + em_smr - \
          (em_rng + em_heng  + em_sbg + em_electrolyzer + em_booster_comp + em_pre_comp), 'Offset'
LP_eps.solve()
print(LP_eps.status)
offset_max = LP_eps.objective.value()

"""
Cost objective model
"""

# CAPEX electrolyzer
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, N_max)]
LP_cost += pulp.lpSum(alpha[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_max)) \
           == CAPEX_electrolyzer

# CAPEX reactor (RNG)
LP_cost += gamma * RNG_max + k * alpha_RNG + C_upgrading * RNG_max == CAPEX_reactor
LP_cost += alpha_RNG <= RNG_max * 10E10
LP_cost += alpha_RNG >= RNG_max * 10E-10

# CAPEX booster compressor (Mobility)
LP_cost += N_booster * CAPEX_booster * 20 == CAPEX_booster_comp

# CAPEX storage
LP_cost += (N_tank * CAPEX_tank + N_prestorage * CAPEX_prestorage) * 20 == CAPEX_storage

# OPEX SBG (electrolyzer + prestorage + tank)
LP_cost += pulp.lpSum(E[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
           pulp.lpSum((H2_direct[str(n)] + H2_tank_in[str(n)]) * C_H2O * WCR for n in input_df.index) + \
           pulp.lpSum(ECF_prestorage * H2_tank_in[str(n)] for n in input_df.index) == OPEX_SBG

# OPEX reactor
LP_cost += pulp.lpSum(CO2[str(n)] * C_CO2 for n in input_df.index) + \
           OPEX_upgrading * RNG_max == OPEX_reactor

# OPEX booster comp
LP_cost += pulp.lpSum(ECF_booster * H2_3[str(n)] * (HOEP[n] + TC) for n in input_df.index) \
           == OPEX_booster_comp

# Cost LP Objective
phi = 0.80
LP_cost += em_ng + em_gas_vehicle + em_smr - \
           (em_rng + em_heng  + em_sbg + em_electrolyzer + em_booster_comp + em_pre_comp) \
           == em_offset
LP_cost += em_offset >= phi * offset_max
LP_cost += (CAPEX_electrolyzer + CAPEX_reactor + CAPEX_booster_comp + CAPEX_storage) + \
           (OPEX_reactor + OPEX_booster_comp + OPEX_SBG) * TVM == total_cost
LP_cost += total_cost, 'Cost'
LP_cost.solve()
print(LP_cost.status)
