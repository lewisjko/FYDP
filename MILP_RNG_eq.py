"""
Model for RNG/NG
"""


'''
Functions for creating variable df and exporting as a csv file
'''
def create_var_df(LP_object):
    variable_tuple = [(v.name,v.varValue) for v in LP_object.variables()]

    variable_df = pd.DataFrame(data=variable_tuple,columns=['variable','value'])

    return variable_df

def export_to_csv(df,filename):
    df.to_csv(filename+'.csv')



import pulp
import pandas as pd
import time
from numpy import count_nonzero

#Estimating the time taken to solve this optimzation problem
#start time
start_time = time.time()
print(start_time)


# Time-series constants
SBG = list(input_df['SBG(kWh)'])
D = list(input_df['NG_demand(m^3)'])
HOEP = list(input_df['HOEP'])
EMF = list(input_df['EMF(tonne/kWh)'])

# Fixed constants
N_max = 3510
N_max += 1
nu_electrolyzer = var['value']['electrolyzer_eff']
E_HHV_H2 = var['value']['E_hhv_h2']
nu_reactor = var['value']['meth_reactor_eff']
HHV_H2 = var['value']['HHV_H2']
HHV_NG = var['value']['HHV_NG']

CO2_available_total = var['value']['CO2_total'] # m^3
CO2_available = float(CO2_available_total) / count_nonzero(SBG)

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

# Tank and compressor constants
Imax = var['value']['Imax'] # kmol
Imin= var['value']['Imin'] # kmol
Fmax_booster = var['value']['Fmax_booster'] # kmol
Fmax_prestorage =var['value']['Fmax_prestorage'] # kmol

CAPEX_prestorage = var['value']['CAPEX_prestorage'] # $
CAPEX_tank = var['value']['CAPEX_tank'] # $

ECF_prestorage = var['value']['ECF_prestorage'] # kWh/kmol H2


# Unit Conversions
#converting the transportation constants to m^3
MW_H2 = var['value']['MW_H2'] #kg/kmol H2
density_H2 = var['value']['density_H2'] #kg/m^3
Imax = Imax * MW_H2 / density_H2 # m^3
Imin = Imin * MW_H2 / density_H2 # m^3
Fmax_prestorage = Fmax_prestorage * MW_H2 / density_H2 # m^3
ECF_prestorage = ECF_prestorage / MW_H2 * density_H2 #kWh/m^3


# RNG model
LP_eps = pulp.LpProblem('LP_eps', pulp.LpMaximize)
LP_cost = pulp.LpProblem('LP_cost', pulp.LpMinimize)

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
RNG = pulp.LpVariable.dicts('RNG',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
CO2 = pulp.LpVariable.dicts('CO2_1',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
NG_1 = pulp.LpVariable.dicts('NG_1',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

em_offset_1 = pulp.LpVariable('em_offset_1',
                          lowBound=0,
                          cat='Continuous')
em_rng = pulp.LpVariable('em_rng',
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

CAPEX_1 = pulp.LpVariable('CAPEX_1', lowBound=0, cat='Continuous')
OPEX_1 = pulp.LpVariable('OPEX_1', lowBound=0, cat='Continuous')

for LP in [LP_eps, LP_cost]:
    for i, h in enumerate([str(i) for i in input_df.index]):
        # Energy and flow constraints
#        if h == '0':
#            LP += RNG_max <= CO2_available
        LP += H2_direct[h] + H2_tank_in[h] == nu_electrolyzer * E_1[h] * E_HHV_H2 ** (-1)
        LP += RNG[h] == nu_reactor * (H2_direct[h] + H2_tank_out[h]) * HHV_H2 / HHV_NG
        LP += CO2[h] == RNG[h]

        # Hydrogen storage tank constraint
        if h == '0':
            LP += I_H2[h] == I_min * N_tank + H2_tank_in[h] - H2_tank_out[h]
        else:
            LP += I_H2[h] == I_H2[str(i - 1)] + H2_tank_in[h] - H2_tank_out[h]
        LP += I_H2[h] <= I_max * N_tank
        LP += I_H2[h] >= I_min * N_tank

        # Demand constraint
        LP += NG_1[h] == D[i] - RNG[h]

        # Electrolyzer and reactor constraints
        LP += N_electrolyzer_1 * E_electrolyzer_min <= E_1[h]
        LP += N_electrolyzer_1 * E_electrolyzer_max >= E_1[h]
        LP += E_1[h] <= SBG[i]
        LP += 0.90 * (H2_direct[h] + H2_tank_out[h]) == 0.10 * NG_1[h]
        LP += RNG[h] <= RNG_max
        if h != '0':
            LP += -RNG_max * tau <= RNG[h] - RNG[str(i - 1)]
            LP += RNG_max * tau >= RNG[h] - RNG[str(i - 1)]
        if h == '0':
            LP += pulp.lpSum(n * alpha_1[str(n)] for n in range(1, N_max)) == N_electrolyzer_1
            LP += pulp.lpSum(alpha_1) <= 1

    # Emission constraints
    LP += pulp.lpSum(EMF_NG * NG_1[h] + EMF_comb * RNG[h] + EMF[int(h)] * E_1[h] + \
                       EMF_bio * CO2[h] + EMF_electrolyzer * (H2_direct[h] + H2_tank_in[h]) + EMF_reactor * RNG[h] + \
                       EMF[int(h)] * ECF_prestorage * H2_tank_in[h] for h in [str(x) for x in input_df.index]) \
        == em_rng
    LP += pulp.lpSum(EMF_NG * D[h] for h in input_df.index) == em_ng

# Eps Objective
LP_eps += em_ng - em_rng, 'Offset_1'
print('eps start')
LP_eps.solve()
offset_max_1 = LP_eps.objective.value()
print(LP_eps.status)

# CAPEX
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, N_max)]
LP_cost += pulp.lpSum(alpha_1[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_max)) + \
           gamma * RNG_max + k + C_upgrading * RNG_max + (N_tank * CAPEX_tank + N_prestorage * CAPEX_prestorage) * 20 \
    == CAPEX_1

# OPEX
LP_cost += pulp.lpSum(CO2[str(n)] * C_CO2 for n in input_df.index) + \
           pulp.lpSum(E_1[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
           pulp.lpSum((H2_direct[str(n)] + H2_tank_in[str(n)]) * C_H2O * WCR for n in input_df.index) + \
           pulp.lpSum((ECF_prestorage * H2_tank_in[str(n)])* (HOEP[n] + TC) for n in input_df.index) + \
           OPEX_upgrading * RNG_max == OPEX_1

# Objectives
##################
phi = 0.5
#################

LP_cost += em_ng - em_rng == em_offset_1
LP_cost += em_offset_1 >= phi * offset_max_1
LP_cost += CAPEX_1 + OPEX_1 * TVM, 'Cost_1'


print('cost start')
LP_cost.solve()



end_time = time.time()

#time difference
time_difference = end_time - start_time


print(time_difference)
print(LP_cost.status)
print(phi)


my_result = create_var_df(LP_cost)
my_result = my_result.append({'variable' : 'LP_cost_status', 'value' : LP_cost.status} , ignore_index=True)
my_result = my_result.append({'variable' : 'solving_time', 'value' : time_difference} , ignore_index=True)
my_result = my_result.append({'variable' : 'offset_max', 'value' : offset_max_1} , ignore_index=True)
my_result = my_result.append({'variable' : 'phi', 'value' : phi} , ignore_index=True)
filename = 'RNG_result_eq_' + str(phi)
export_to_csv(my_result, filename)

