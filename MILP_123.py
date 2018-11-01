"""
Combined model
"""

import pulp

input_df = None
var = None

# Time-series constants
SBG = list(input_df['SBG(kWh)'])
NG_demand = list(input_df['NG_demand(m^3)'])
mobility_demand = list(input_df['mobility_demand(m^3)']) #m^3
HOEP = list(input_df['HOEP'])
EMF = list(input_df['EMF(tonne/kWh)'])

# Fixed constants
N_max = 1000
N_max += 1
nu_electrolyzer = var['value']['electrolyzer_eff']
E_HHV_H2 = var['value']['E_hhv_h2'] #kwh/m^3
nu_reactor = var['value']['meth_reactor_eff']
HHV_H2 = var['value']['HHV_H2'] #MMBtu/kmol
HHV_NG = var['value']['HHV_NG'] #MMBtu/kmol
CO2_available = var['value']['CO2_available'] #m^3/h
E_electrolyzer_min = var['value']['min_E_cap'] #kwh
E_electrolyzer_max = var['value']['max_E_cap'] #kwh
tau = 0.50

EMF_NG = var['value']['EMF_NG'] #tonne CO2/m^3 H2
EMF_comb = var['value']['EMF_combRNG'] #tonne CO2/m^3 RNG
EMF_bio = var['value']['EMF_bioCO2'] #tonne CO2/kWh
EMF_electrolyzer = var['value']['EMF_electrolyzer'] #tonne CO2/m^3 H2
EMF_reactor = var['value']['EMF_reactor'] #tonne CO2/m^3 RNG
EMF_vehicle = var['Value']['emission_gasoline_v'] #tonne CO2/car/year
num_vehicle = var['Value']['N_gasoline_v']
FCV_penetration = var['Value']['FCV_penetration'] #FCV market penetration (estimated)

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

# Fixed constants for transportation models
Imax = var['value']['Imax'] #kmol
Imin= var['value']['Imin'] #kmol
Fmax_booster = var['value']['Fmax_booster'] #kmol
Fmax_prestorage =var['value']['Fmax_prestorage'] #kmol

CAPEX_booster = var['value']['CAPEX_booster'] #$
CAPEX_prestorage = var['value']['CAPEX_prestorage'] #$
CAPEX_tank = var['value']['CAPEX_tank'] #$

ECF_prestorage = var['value']['ECF_prestorage'] #kWh/kmol H2

z_booster = var['value']['z_booster'] #compressibility factor for booster compressor
R = var['value']['R'] #kJ/kmolK
T = var['value']['T'] #K
comp_efficiency = var['value']['comp_efficiency']  #isentropic compressor efficiency
heat_cap_ratio = var['value']['heat_cap_ratio'] #heat capacaity ratio of hydrogen
P_in_booster = var['value']['P_in_booster'] #inlet pressure of booster compressor
P_out_booster = var['value']['P_out_booster'] #outlet pressure of booster compressor
N_stage_booster = var['value']['N_stage_booster']

ECF_booster = z_booster * R * T * N_stage_booster / comp_efficiency * \
              heat_cap_ratio / (heat_cap_ratio - 1) * \
              (((P_out_booster / P_in_booster) ** ((heat_cap_ratio - 1) / N_stage_booster / heat_cap_ratio ))-1) \
              / 3600 #converting kJ to kWh ECF booster in kWh/kmol

# Converting the transportation constants to m^3
MW_H2 = var['value']['MW_H2'] #kg/kmol H2
density_H2 = var['value']['density_H2'] #kg/m^3

Imax = Imax * MW_H2 / density_H2 # m^3
Imin = Imin * MW_H2 / density_H2 # m^3
Fmax_booster = Fmax_booster * MW_H2 / density_H2 # m^3
Fmax_prestorage = Fmax_prestorage * MW_H2 / density_H2 # m^3

ECF_booster = ECF_booster / MW_H2 * density_H2 #kWh/m^3
ECF_prestorage = ECF_prestorage / MW_H2 * density_H2 #kWh/m^3

# RNG + HENG model
LP_eps = pulp.LpProblem('LP_eps', pulp.LpMaximize)
LP_cost = pulp.LpProblem('LP_cost', pulp.LpMinimize)

# RNG Variables
RNG_max = pulp.LpVariable('RNG_max',
                          lowBound=0,
                          cat='Continuous')
N_electrolyzer_1 = pulp.LpVariable('N_electrolyzer_1',
                                 lowBound=0,
                                 cat='Integer')
alpha_1 = pulp.LpVariable.dicts('alpha_1',
                          [str(i) for i in range(1, N_max)],
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
CAPEX_2 = pulp.LpVariable('CAPEX_2', lowBound=0, cat='Continuous')
OPEX_2 = pulp.LpVariable('OPEX_2', lowBound=0, cat='Continuous')

# Transportation Variables
N_electrolyzer_3 = pulp.LpVariable('N_electrolyzer_3',
                          lowBound=0,
                          cat='Integer')

N_booster = pulp.LpVariable('N_booster',
                          lowBound=0,
                          cat='Integer')

N_prestorage = pulp.LpVariable('N_prestorage',
                          lowBound=0,
                          cat='Integer')

N_tank = pulp.LpVariable('N_tank',
                          lowBound=0,
                          cat='Integer')


alpha_3 = pulp.LpVariable.dicts('alpha_3',
                          [str(i) for i in range(1, N_max+1)],
                          cat='Binary')


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

CAPEX_3 = pulp.LpVariable('CAPEX_3', lowBound=0, cat='Continuous')
OPEX_3 = pulp.LpVariable('OPEX_3', lowBound=0, cat='Continuous')

# Shared Variables
CO2 = pulp.LpVariable.dicts('CO2_2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
NG = pulp.LpVariable.dicts('NG',
                          [str(i) for i in input_df.index],
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
em_offset_fcv = pulp.LpVariable('em_offset_fcv',
                          lowBound=0,
                          cat='Continuous')

for LP in [LP_eps, LP_cost]:
    for i, h in enumerate([str(i) for i in input_df.index]):
        # Energy and flow constraints
        LP += H2_1[h] == nu_electrolyzer * E_1[h] / E_HHV_H2
        LP += H2_2[h] == nu_electrolyzer * E_2[h] / E_HHV_H2
        LP += H2_3[h] == nu_electrolyzer * E_3[h] / E_HHV_H2

        if h == '0':
            LP += RNG_max <= CO2_available
        LP += RNG[h] == nu_reactor * H2_1[h] * HHV_H2 / HHV_NG
        LP += CO2[h] == RNG[h]
        LP += H2_3[h] == H2_tank_in[h] + H2_direct[h]

        # NG demand constraint
        LP += NG_demand[i] == RNG[h] + NG[h] + H2_2[h]

        # Mobility demand constraint
        LP += H2_tank_out[h] + H2_direct[h] == mobility_demand[i]

        # Supply constraint
        LP += E_1[h] + E_2[h] + E_3[h] <= SBG[i]

        # Hydrogen storage inventory constraint
        if h == '0': # at hour zero, no accumulation exists
            LP += I_H2[h] >= Imin * N_tank + H2_tank_in[h] - H2_tank_out[h]
            LP += I_H2[h] <= Imax * N_tank + H2_tank_in[h] - H2_tank_out[h]
        else: # else, accumulation exists
            LP += I_H2[h] == I_H2[str(i - 1)] + H2_tank_in[h] - H2_tank_out[h]
        LP += I_H2[h] <= Imax * N_tank
        LP += I_H2[h] >= Imin * N_tank

        # Compressor capacity constraint
        LP += H2_tank_in[h] <= N_prestorage * Fmax_prestorage
        LP += H2_tank_out[h] + H2_direct[h] <= N_booster * Fmax_booster

        # Electrolyzer constraints
        LP += N_electrolyzer_1 * E_electrolyzer_min <= E_1[h]
        LP += N_electrolyzer_1 * E_electrolyzer_max >= E_1[h]
        LP += N_electrolyzer_2 * E_electrolyzer_min <= E_2[h]
        LP += N_electrolyzer_2 * E_electrolyzer_max >= E_2[h]
        LP += N_electrolyzer_3 * E_electrolyzer_min <= E_3[h]
        LP += N_electrolyzer_3 * E_electrolyzer_max >= E_3[h]

        # Reactor constraints
        LP += RNG[h] <= RNG_max
        if h != '0':
            LP += -RNG_max * tau <= RNG[h] - RNG[str(i - 1)]
            LP += RNG_max * tau >= RNG[h] - RNG[str(i - 1)]

        # HENG-specific constraints
        LP += 0.95 * H2_2[h] <= 0.05 * (NG[h] + RNG[h])

    # Integer constraints
    LP += pulp.lpSum(n * alpha_1[str(n)] for n in range(1, N_max)) == N_electrolyzer_1
    LP += pulp.lpSum(alpha_1) <= 1
    LP += pulp.lpSum(n * alpha_2[str(n)] for n in range(1, N_max)) == N_electrolyzer_2
    LP += pulp.lpSum(alpha_2) <= 1
    LP += pulp.lpSum(n * alpha_3[str(n)] for n in range(1, N_max + 1)) == N_electrolyzer_3
    LP += pulp.lpSum(alpha_3) <= 1

    # Emission constraints
    LP += pulp.lpSum(EMF_comb * RNG[h] + EMF[h] * E_1[h] + EMF_bio * CO2[h] + \
                     EMF_electrolyzer * H2_1[h] + EMF_reactor * RNG[h] \
                     for h in [str(x) for x in input_df.index]) == em_rng
    LP += pulp.lpSum(EMF_NG * NG[h] + EMF[h] * E_2[h] + EMF_electrolyzer * H2_2[h] \
                     for h in [str(x) for x in input_df.index]) == em_heng
    LP += num_vehicle * FCV_penetration * EMF_vehicle == em_offset_fcv
    LP += pulp.lpSum(EMF_NG * NG_demand[h] for h in input_df.index) == em_ng

# Epsilon LP Objective
LP_eps += em_ng - em_rng - em_heng + em_offset_fcv, 'Offset'
LP_eps.solve()
print(LP_eps.status)
offset_max = LP_eps.objective.value()

# Electrolyzer cost
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, N_max)]

# CAPEX RNG
LP_cost += pulp.lpSum(alpha_1[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_max)) + \
           gamma * RNG_max + k + C_upgrading * RNG_max == CAPEX_1

# CAPEX HENG
LP_cost += pulp.lpSum(alpha_2[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_max)) == CAPEX_2

# CAPEX Mobility
LP_cost += pulp.lpSum(alpha_3[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_max + 1)) + \
           (N_booster * CAPEX_booster + N_prestorage * CAPEX_prestorage + N_tank * CAPEX_tank) * TVM == CAPEX_3

# OPEX RNG
LP_cost += pulp.lpSum(CO2[str(n)] * C_CO2 for n in input_df.index) + \
           pulp.lpSum(E_1[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
           pulp.lpSum(H2_1[str(n)] * C_H2O * WCR for n in input_df.index) + \
           OPEX_upgrading * RNG_max == OPEX_1

# OPEX HENG
LP_cost += pulp.lpSum(E_2[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
           pulp.lpSum(H2_2[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_2

# OPEX Mobility
LP_cost += pulp.lpSum((E_3[str(n)] + \
                       ECF_booster * (H2_tank_out[str(n)] + H2_direct[str(n)]) + \
                       ECF_prestorage * H2_tank_in[str(n)]) * (HOEP[n] + TC) for n in input_df.index) + \
           pulp.lpSum(H2_3[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_3

# Cost LP Objective
phi = 0.80
LP_cost += em_ng - em_rng - em_heng + em_offset_fcv == em_offset
LP_cost += em_offset >= phi * offset_max
LP_cost += (CAPEX_1 + CAPEX_2 + CAPEX_3) + (OPEX_1 + OPEX_2 + OPEX_3) * TVM, 'Cost'
LP_cost.solve()
print(LP_cost.status)
