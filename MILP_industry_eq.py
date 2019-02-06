'''
Model for Industry
'''
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
SBG = list(input_df['SBG(kWh)']) #kwh
industry_demand = list(input_df['industry_demand(m^3)']) #m^3
HOEP = list(input_df['HOEP']) #$/kwh
EMF = list(input_df['EMF(tonne/kWh)'])

# Fixed constants for other models
nu_electrolyzer = var['value']['electrolyzer_eff'] #dimensionless
E_HHV_H2 = var['value']['E_hhv_h2'] #kwh/m^3
nu_reactor = var['value']['meth_reactor_eff'] #dimensionless
HHV_H2 = var['value']['HHV_H2'] #MMBtu/kmol
HHV_NG = var['value']['HHV_NG'] #MMBtu/kmol
# CO2_available = var['value']['CO2_available'] #m^3/h
E_electrolyzer_min = var['value']['min_E_cap'] #kwh
E_electrolyzer_max = var['value']['max_E_cap'] #kwh
tau = 0.50


EMF_NG = var['value']['EMF_NG'] #tonne CO2/m^3 H2
EMF_comb = var['value']['EMF_combRNG'] #tonne CO2 /m^3 H2
EMF_nuc = var['value']['EMF_nuclear'] #tonne CO2 / kWh
EMF_bio = var['value']['EMF_bioCO2'] #tonne CO2/m^3 bio CO2
EMF_electrolyzer = var['value']['EMF_electrolyzer'] #tonne CO2 /m^3 H2
EMF_reactor = var['value']['EMF_reactor'] #tonne CO2 /m^3 RNG produced
EMF_SMR=var['value']['EMF_SMR'] #kg CO2/kmol of H2

beta = var['value']['beta']
C_0 = var['value']['C_0'] #$/kW
mu = var['value']['mu']
gamma = var['value']['gamma']
k = var['value']['k'] #$
C_upgrading = var['value']['C_upgrading'] #$/m^3 reactor capacity
C_CO2 = var['value']['C_CO2']#$/m^3 CO2
TC = var['value']['TC'] #$/kWh
C_H2O = var['value']['C_H2O'] #$/L
WCR = var['value']['water_cons_rate'] #L H2o/m^3 H2
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

#converting the transportation constants to m^3
MW_H2 = var['value']['MW_H2'] #kg/kmol H2
density_H2 = var['value']['density_H2'] #kg/m^3

Imax = Imax * MW_H2 / density_H2 # m^3
Imin = Imin * MW_H2 / density_H2 # m^3

Fmax_prestorage = Fmax_prestorage * MW_H2 / density_H2 # m^3
ECF_prestorage = ECF_prestorage / MW_H2 * density_H2 #kWh/m^3

EMF_SMR=EMF_SMR*0.001/MW_H2*density_H2 #tonne CO2/m^3 of H2

#number of electrolyzer max
N_electrolyzer_max = int(3510)



# Transportation model
LP_eps_4 = pulp.LpProblem('LP_eps_4', pulp.LpMaximize)
LP_cost_4 = pulp.LpProblem('LP_cost_4', pulp.LpMinimize)


N_electrolyzer_4 = pulp.LpVariable('N_electrolyzer_4',
                          lowBound=0,
                          cat='Integer')



N_prestorage_4 = pulp.LpVariable('N_prestorage_4',
                          lowBound=0,
                          cat='Integer')

N_tank_4 = pulp.LpVariable('N_tank_4',
                          lowBound=0,
                          cat='Integer')


alpha_4 = pulp.LpVariable.dicts('alpha_4',
                          [str(i) for i in range(1, N_electrolyzer_max+1)],
                          cat='Binary')


E_4 = pulp.LpVariable.dicts('E_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

H2_4 = pulp.LpVariable.dicts('H2_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

#hydrogen delivered directly to the booster compressor
H2_direct_4 = pulp.LpVariable.dicts('H2_direct_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

#hydrogen going into the tank
H2_tank_in_4 = pulp.LpVariable.dicts('H2_tank_in_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
#hydrogen leaving the tank
H2_tank_out_4 = pulp.LpVariable.dicts('H2_tank_out_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

#hydrogen inventory in the tank
I_H2_4 = pulp.LpVariable.dicts('I_H2_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

em_offset_4 = pulp.LpVariable('em_offset_4',
                          lowBound=0,
                          cat='Continuous')
em_compressor_4 = pulp.LpVariable('em_compressor_4',
                          lowBound=0,
                          cat='Continuous')
em_electrolyzer_4 = pulp.LpVariable('em_electrolyzer_4',
                          lowBound=0,
                          cat='Continuous')
em_before = pulp.LpVariable('em_before',
                          lowBound=0,
                          cat='Continuous')
em_sbg_4 = pulp.LpVariable('em_sbg_4',
                          lowBound=0,
                          cat='Continuous')

CAPEX_4 = pulp.LpVariable('CAPEX_4', lowBound=0, cat='Continuous')
OPEX_4 = pulp.LpVariable('OPEX_4', lowBound=0, cat='Continuous')

for LP_4 in [LP_eps_4, LP_cost_4]:
    for i, h in enumerate([str(i) for i in input_df.index]):

        # Energy and flow constraints
        LP_4 += H2_4[h] == nu_electrolyzer * E_4[h] / E_HHV_H2
        LP_4 += H2_4[h] == H2_tank_in_4[h] + H2_direct_4[h]

        #hydrogen storage inventory constraint
        if h == '0':
            LP_4 += I_H2_4[h] == Imin * N_tank_4 + H2_tank_in_4[h] - H2_tank_out_4[h]
        else:
            LP_4 += I_H2_4[h] == I_H2_4[str(i-1)] + H2_tank_in_4[h] - H2_tank_out_4[h]

        # Demand constraint
        LP_4 += H2_tank_out_4[h] + H2_direct_4[h] == industry_demand[i]

        # Electrolyzer constraints
        LP_4 += N_electrolyzer_4 * E_electrolyzer_min <= E_4[h]
        LP_4 += N_electrolyzer_4 * E_electrolyzer_max >= E_4[h]
        LP_4 += E_4[h] <= SBG[i]

        #storage inventory constraint
        LP_4 += I_H2_4[h] <= Imax * N_tank_4
        LP_4 += I_H2_4[h] >= Imin * N_tank_4

        #compressor capacity constraint
        LP_4 += H2_tank_in_4[h] <= N_prestorage_4 * Fmax_prestorage

    #Number of eletrolyzer constraint
    LP_4 += pulp.lpSum(n * alpha_4[str(n)] for n in range(1, N_electrolyzer_max+1)) == N_electrolyzer_4
    LP_4 += pulp.lpSum(alpha_4) <= 1

    #Emission constraints
    LP_4 += pulp.lpSum(EMF_SMR * industry_demand[h] for h in input_df.index) == em_before

    #Emission calculation
    LP_4 += pulp.lpSum(EMF[n] * (ECF_prestorage * H2_tank_in_4[str(n)]) for n in input_df.index)  == em_compressor_4

    LP_4 += pulp.lpSum(EMF_electrolyzer * H2_4[h] for h in [str(x) for x in input_df.index]) == em_electrolyzer_4
    LP_4 += pulp.lpSum(EMF[int(h)] * (E_4[h]) for h in [str(x) for x in input_df.index]) == em_sbg_4



# Objectives
LP_eps_4 += em_before - em_compressor_4 - em_electrolyzer_4 - em_sbg_4, 'Offset_4'
print('eps start')
LP_eps_4.solve()
print(LP_eps_4.status)
offset_max_4 = LP_eps_4.objective.value()

# Cost of electrolyzer list
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, N_electrolyzer_max+1)]

# CAPEX
LP_cost_4 += pulp.lpSum(alpha_4[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_electrolyzer_max+1)) + \
        (N_prestorage_4 * CAPEX_prestorage + \
        N_tank_4 * CAPEX_tank) * 20 == CAPEX_4


LP_cost_4 += pulp.lpSum((E_4[str(n)] + \
                    ECF_prestorage * H2_tank_in_4[str(n)]) * (HOEP[n] + TC) for n in input_df.index) + \
        pulp.lpSum(H2_4[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_4

#################
phi = 0.5
################


LP_cost_4 += em_before - em_compressor_4 - em_electrolyzer_4 - em_sbg_4 == em_offset_4
LP_cost_4 += em_offset_4 >= phi * offset_max_4
LP_cost_4 += CAPEX_4 + OPEX_4 * TVM, 'Cost_4'


print('cost start')
LP_cost_4.solve()


end_time = time.time()

#time difference
time_difference = end_time - start_time

print(time_difference)
print(LP_cost_4.status)
print(phi)

my_result = create_var_df(LP_cost_4)
my_result = my_result.append({'variable' : 'LP_cost_status', 'value' : LP_cost_4.status} , ignore_index=True)
my_result = my_result.append({'variable' : 'solving_time', 'value' : time_difference} , ignore_index=True)
my_result = my_result.append({'variable' : 'offset_max', 'value' : offset_max_4} , ignore_index=True)
my_result = my_result.append({'variable' : 'phi', 'value' : phi} , ignore_index=True)
filename = 'industry_result_eq_' + str(phi)
export_to_csv(my_result,filename)
