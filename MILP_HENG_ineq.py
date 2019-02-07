"""
Model for HENG
"""

import pulp
import pandas as pd
import time

'''
Functions for creating variable df and exporting as a csv file
'''
def create_var_df(LP_object):
    variable_tuple = [(v.name,v.varValue) for v in LP_object.variables()]

    variable_df = pd.DataFrame(data=variable_tuple,columns=['variable','value'])

    return variable_df

def export_to_csv(df,filename):
    df.to_csv(filename+'.csv')



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
Imax = var['value']['Imax'] # kmol
Imin= var['value']['Imin'] # kmol
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



# HENG model
LP_cost = pulp.LpProblem('LP_cost', pulp.LpMinimize)
LP_eps = pulp.LpProblem('LP_eps', pulp.LpMaximize)
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
CO2 = pulp.LpVariable.dicts('CO2_2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
NG_2 = pulp.LpVariable.dicts('NG_2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

em_offset_2 = pulp.LpVariable('em_offset_2',
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

CAPEX_2 = pulp.LpVariable('CAPEX_2', lowBound=0, cat='Continuous')
OPEX_2 = pulp.LpVariable('OPEX_2', lowBound=0, cat='Continuous')

# Total cost
total_cost = pulp.LpVariable('total_cost', lowBound=0, cat='Continuous')


for LP in [LP_eps, LP_cost]:
    for i, h in enumerate([str(i) for i in input_df.index]):
        # Energy and flow constraints
        LP += H2_direct[h] + H2_tank_in[h] == nu_electrolyzer * E_2[h] * E_HHV_H2 ** (-1)

        # Hydrogen storage tank constraint
        if h == '0':
            LP += I_H2[h] == Imin * N_tank + H2_tank_in[h] - H2_tank_out[h]
        else:
            LP += I_H2[h] == I_H2[str(i - 1)] + H2_tank_in[h] - H2_tank_out[h]
        LP += I_H2[h] <= Imax * N_tank
        LP += I_H2[h] >= Imin * N_tank

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
LP_eps += em_ng - em_heng, 'Offset_2'
print('eps start')
LP_eps.solve()
offset_max_2 = LP_eps.objective.value()
print(LP_eps.status)

# CAPEX
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, N_max)]
LP_cost += pulp.lpSum(alpha_2[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_max)) + \
           (N_tank * CAPEX_tank + N_prestorage * CAPEX_prestorage) * 20 == CAPEX_2

# OPEX
LP_cost += pulp.lpSum(E_2[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
           pulp.lpSum((H2_direct[str(n)] + H2_tank_in[str(n)]) * C_H2O * WCR for n in input_df.index) + \
           pulp.lpSum((ECF_prestorage * H2_tank_in[str(n)])* (HOEP[n] + TC) for n in input_df.index) == OPEX_2

# Objectives

######################
phi = 0.5
######################

LP_cost += em_ng - em_heng == em_offset_2
LP_cost += em_offset_2 >= phi * offset_max_2

LP_cost += CAPEX_2 + OPEX_2 * TVM == total_cost
LP_cost += total_cost,'Cost_2'



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
my_result = my_result.append({'variable' : 'offset_max', 'value' : offset_max_2} , ignore_index=True)
my_result = my_result.append({'variable' : 'phi', 'value' : phi} , ignore_index=True)
filename = 'HENG_result_ineq_' + str(phi)
export_to_csv(my_result, filename)
