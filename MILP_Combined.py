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

# Fixed constants
N_max = 30000
N_max += 1
nu_electrolyzer = var['value']['electrolyzer_eff']
E_HHV_H2 = var['value']['E_hhv_h2'] #kwh/m^3
nu_reactor = var['value']['meth_reactor_eff']
HHV_H2 = var['value']['HHV_H2'] #MMBtu/kmol
HHV_NG = var['value']['HHV_NG'] #MMBtu/kmol
CO2_available_total = var['value']['CO2_total'] #m^3/year
CO2_available = float(CO2_available_total) / count_nonzero(SBG)
E_electrolyzer_min = var['value']['min_E_cap'] #kwh
E_electrolyzer_max = var['value']['max_E_cap'] #kwh
tau = 0.50

EMF_NG = var['value']['EMF_NG'] #tonne CO2/m^3 H2
EMF_comb = var['value']['EMF_combRNG'] #tonne CO2/m^3 RNG
EMF_bio = var['value']['EMF_bioCO2'] #tonne CO2/kWh
EMF_electrolyzer = var['value']['EMF_electrolyzer'] #tonne CO2/m^3 H2
EMF_reactor = var['value']['EMF_reactor'] #tonne CO2/m^3 RNG
EMF_vehicle = var['value']['emission_gasoline_v'] #tonne CO2/car/year
num_vehicle = var['value']['N_gasoline_v']
FCV_penetration = var['value']['FCV_penetration'] #FCV market penetration (estimated)
EMF_SMR = var['value']['EMF_SMR'] #kg CO2/kmol of H2

beta = var['value']['beta']
C_0 = var['value']['C_0'] #$/kW
mu = var['value']['mu']
gamma = var['value']['gamma']
k = var['value']['k'] #$
C_upgrading = var['value']['C_upgrading'] #$/m^3 reactor capacity
C_CO2 = var['value']['C_CO2'] #$/m^3 CO2
TC = var['value']['TC'] #$/kWh
C_H2O = var['value']['C_H2O'] #$/L
WCR = var['value']['water_cons_rate'] #L H2O/m^3 H2
OPEX_upgrading = var['value']['OPEX_upgrading'] #$/m^3 reactor capacity
TVM = var['value']['TVM']

# Fixed constants for transportation + industry models
I_max = var['value']['I_max'] #kmol
I_min= var['value']['I_min'] #kmol
F_max_booster = var['value']['F_max_booster'] #kmol
F_max_prestorage =var['value']['F_max_prestorage'] #kmol

CAPEX_booster = var['value']['CAPEX_booster'] #$
CAPEX_prestorage = var['value']['CAPEX_prestorage'] #$
CAPEX_tank = var['value']['CAPEX_tank'] #$

ECF_prestorage = var['value']['ECF_prestorage'] #kWh/kmol H2

z_booster = var['value']['z_booster'] #compressibility factor for booster compressor
R = var['value']['R'] #kJ/kmolK
T = var['value']['T'] #K
comp_efficiency = var['value']['comp_efficiency']  #isentropic compressor efficiency
heat_cap_ratio = var['value']['heat_cap_ratio'] #heat capacaity ratio of hydrogen
P_iN_booster_3 = var['value']['P_iN_booster_3'] #inlet pressure of booster compressor
P_out_booster = var['value']['P_out_booster'] #outlet pressure of booster compressor
N_stage_booster = var['value']['N_stage_booster']

ECF_booster = z_booster * R * T * N_stage_booster / comp_efficiency * \
              heat_cap_ratio / (heat_cap_ratio - 1) * \
              (((P_out_booster / P_iN_booster_3) ** ((heat_cap_ratio - 1) / N_stage_booster / heat_cap_ratio ))-1) \
              / 3600 #converting kJ to kWh ECF booster in kWh/kmol

# Converting the transportation constants to m^3
MW_H2 = var['value']['MW_H2'] #kg/kmol H2
density_H2 = var['value']['density_H2'] #kg/m^3

I_max = I_max * MW_H2 / density_H2 # m^3
I_min = I_min * MW_H2 / density_H2 # m^3
F_max_booster = F_max_booster * MW_H2 / density_H2 # m^3
F_max_prestorage = F_max_prestorage * MW_H2 / density_H2 # m^3

ECF_booster = ECF_booster / MW_H2 * density_H2 #kWh/m^3
ECF_prestorage = ECF_prestorage / MW_H2 * density_H2 #kWh/m^3
EMF_SMR = EMF_SMR * 0.001 / MW_H2 * density_H2 # tonne CO2/m^3 of H2

# RNG + HENG model
LP_eps = pulp.LpProblem('LP_eps', pulp.LpMaximize)
LP_cost = pulp.LpProblem('LP_cost', pulp.LpMinimize)

# RNG Variables
RNG_max = pulp.LpVariable('RNG_max',
                          lowBound=0,
                          cat='Continuous')
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
alpha_RNG = pulp.LpVariable('alpha_RNG',
                          cat='Binary')

CAPEX_reactor = pulp.LpVariable('CAPEX_reactor', lowBound=0, cat='Continuous')
OPEX_1 = pulp.LpVariable('OPEX_1', lowBound=0, cat='Continuous')

# HENG Variables
E_2 = pulp.LpVariable.dicts('E_2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
H2_2 = pulp.LpVariable.dicts('H2_2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

OPEX_2 = pulp.LpVariable('OPEX_2', lowBound=0, cat='Continuous')

# Transportation Variables
N_booster_3 = pulp.LpVariable('N_booster_3',
                          lowBound=0,
                          cat='Integer')
N_prestorage = pulp.LpVariable('N_prestorage',
                          lowBound=0,
                          cat='Integer')
N_tank = pulp.LpVariable('N_tank',
                          lowBound=0,
                          cat='Integer')
E_3 = pulp.LpVariable.dicts('E_3',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
H2_3 = pulp.LpVariable.dicts('H2_3',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
H2_direct = pulp.LpVariable.dicts('H2_direct',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
H2_tank_in = pulp.LpVariable.dicts('H2_tank_in',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
H2_tank_out = pulp.LpVariable.dicts('H2_tank_out',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
I_H2 = pulp.LpVariable.dicts('I_H2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

CAPEX_mobility = pulp.LpVariable('CAPEX_mobility', lowBound=0, cat='Continuous')
OPEX_3 = pulp.LpVariable('OPEX_3', lowBound=0, cat='Continuous')

# Industry Variables
N_booster_4 = pulp.LpVariable('N_booster_4',
                          lowBound=0,
                          cat='Integer')
E_4 = pulp.LpVariable.dicts('E_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
H2_4 = pulp.LpVariable.dicts('H2_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

CAPEX_industry = pulp.LpVariable('CAPEX_industry', cat='Continuous')
OPEX_4 = pulp.LpVariable('OPEX_4', cat='Continuous')

# Shared Variables
N_electrolyzer = pulp.LpVariable('N_electrolyzer_2',
                          lowBound=0,
                          cat='Integer')
alpha = pulp.LpVariable.dicts('alpha_2',
                          [str(i) for i in range(1, N_max)],
                          cat='Binary')
CO2 = pulp.LpVariable.dicts('CO2_2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
NG = pulp.LpVariable.dicts('NG',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
CAPEX_electrolyzer = pulp.LpVariable('CAPEX_Electrolyzer',
                          lowBound=0,
                          cat='Continuous')

# Emission Variables
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

for LP in [LP_eps, LP_cost]:
    for i, h in enumerate([str(i) for i in input_df.index]):
        # Energy and flow constraints
        LP += H2_1[h] == nu_electrolyzer * E_1[h] / E_HHV_H2
        LP += H2_2[h] == nu_electrolyzer * E_2[h] / E_HHV_H2
        LP += H2_3[h] == nu_electrolyzer * E_3[h] / E_HHV_H2
        LP += H2_4[h] == nu_electrolyzer * E_4[h] / E_HHV_H2

        if h == '0':
            LP += RNG_max <= CO2_available
        LP += RNG[h] == nu_reactor * H2_1[h] * HHV_H2 / HHV_NG
        LP += CO2[h] == RNG[h]
        LP += H2_3[h] == H2_tank_in[h] + H2_direct[h]

        # NG demand constraint
        LP += NG_demand[i] == RNG[h] + NG[h] + H2_2[h]

        # Mobility demand constraint
        LP += H2_tank_out[h] + H2_direct[h] == mobility_demand[i]

        # Industry demand constraint
        LP += H2_4[h] == industry_demand[i]

        # Supply constraint
        LP += E_1[h] + E_2[h] + E_3[h] + E_4[h] <= SBG[i]

        # Hydrogen storage inventory constraint
        if h == '0': # at hour zero, no accumulation exists
            LP += I_H2[h] == I_min * N_tank + H2_tank_in[h] - H2_tank_out[h]
            LP += I_H2[h] <= I_max * N_tank + H2_tank_in[h] - H2_tank_out[h]
        else: # else, accumulation exists
            LP += I_H2[h] == I_H2[str(i - 1)] + H2_tank_in[h] - H2_tank_out[h]
        LP += I_H2[h] <= I_max * N_tank
        LP += I_H2[h] >= I_min * N_tank

        # Compressor capacity constraint
        LP += H2_tank_in[h] <= N_prestorage * F_max_prestorage
        LP += H2_tank_out[h] + H2_direct[h] <= N_booster_3 * F_max_booster
        LP += H2_4[h] <= N_booster_4 * F_max_booster

        # Electrolyzer constraints
        LP += N_electrolyzer * E_electrolyzer_min <= E_1[h] + E_2[h] + E_3[h] + E_4[h]
        LP += N_electrolyzer * E_electrolyzer_max >= E_1[h] + E_2[h] + E_3[h] + E_4[h]

        # Reactor constraints
        LP += RNG[h] <= RNG_max
        if h != '0':
            LP += -RNG_max * tau <= RNG[h] - RNG[str(i - 1)]
            LP += RNG_max * tau >= RNG[h] - RNG[str(i - 1)]

        # HENG-specific constraints
        LP += 0.95 * H2_2[h] <= 0.05 * (NG[h] + RNG[h])

    # Integer constraints
    LP += pulp.lpSum(n * alpha[str(n)] for n in range(1, N_max)) == N_electrolyzer
    LP += pulp.lpSum(alpha) <= 1

    # Emission constraints
    LP += pulp.lpSum(EMF_comb * RNG[h] + EMF_bio * CO2[h] + EMF_reactor * RNG[h] \
                     for h in [str(x) for x in input_df.index]) == em_rng
    LP += pulp.lpSum(EMF_NG * NG[h] for h in [str(x) for x in input_df.index]) == em_heng
    em_gas_vehicle = num_vehicle * FCV_penetration * EMF_vehicle
    LP += pulp.lpSum(EMF_NG * NG_demand[h] for h in input_df.index) == em_ng
    LP += pulp.lpSum(EMF_SMR * industry_demand[h] for h in input_df.index) == em_smr
    LP += pulp.lpSum(EMF[int(h)] * (E_1[h] + E_2[h] + E_3[h] + E_4[h]) \
                     for h in [str(x) for x in input_df.index]) == em_sbg
    LP += pulp.lpSum(EMF_electrolyzer * (H2_1[h] + H2_2[h] + H2_3[h] + H2_4[h]) \
                     for h in [str(x) for x in input_df.index]) == em_electrolyzer
    LP += pulp.lpSum(EMF[n] * (ECF_booster * (H2_tank_out[str(n)] + H2_direct[str(n)] + H2_4[h])) \
                     for n in input_df.index) == em_booster_comp
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

# CAPEX Electrolyzer
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, N_max)]
LP_cost += pulp.lpSum(alpha[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_max)) \
           == CAPEX_electrolyzer

# CAPEX Reactor (RNG)
LP_cost += gamma * RNG_max + k * alpha_RNG + C_upgrading * RNG_max == CAPEX_reactor
LP_cost += alpha_RNG <= RNG_max * 10E10
LP_cost += alpha_RNG >= RNG_max / 10E10

# CAPEX Compressors + Tank (Mobility)
LP_cost += (N_booster_3 * CAPEX_booster + N_prestorage * CAPEX_prestorage + N_tank * CAPEX_tank) * TVM == CAPEX_mobility

# CAPEX Compressor (Industry)
LP_cost += N_booster_4 * CAPEX_booster * TVM == CAPEX_industry

# OPEX RNG
LP_cost += pulp.lpSum(CO2[str(n)] * C_CO2 for n in input_df.index) + \
           pulp.lpSum(E_1[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
           pulp.lpSum(H2_1[str(n)] * C_H2O * WCR for n in input_df.index) + \
           OPEX_upgrading * RNG_max == OPEX_1

# OPEX HENG
LP_cost += pulp.lpSum(E_2[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
           pulp.lpSum(H2_2[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_2

# OPEX Mobility
LP_cost += pulp.lpSum((E_3[str(n)] + ECF_booster * (H2_tank_out[str(n)] + H2_direct[str(n)]) + \
                       ECF_prestorage * H2_tank_in[str(n)]) * (HOEP[n] + TC) for n in input_df.index) + \
           pulp.lpSum(H2_3[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_3

# OPEX Industry
LP_cost += pulp.lpSum((E_4[str(n)] + ECF_booster * (H2_4[str(n)])) * (HOEP[n] + TC) for n in input_df.index) + \
           pulp.lpSum(H2_4[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_4


# Cost LP Objective
phi = 0.80
LP_cost += em_ng - em_rng - em_heng + em_gas_vehicle - em_electrolyzer == em_offset
LP_cost += em_offset >= phi * offset_max
LP_cost += (CAPEX_electrolyzer + CAPEX_reactor + CAPEX_mobility + CAPEX_industry) + \
           (OPEX_1 + OPEX_2 + OPEX_3 + OPEX_4) * TVM, 'Cost'
LP_cost.solve()
print(LP_cost.status)
