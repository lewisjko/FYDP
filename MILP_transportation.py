"""
Model for Transportation
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

# Time-series constants
SBG = list(input_df['SBG(kWh)']) #kwh
mobility_demand = list(input_df['mobility_demand(m^3)']) #m^3
HOEP = list(input_df['HOEP']) #$/kwh
EMF = list(input_df['EMF(tonne/kWh)'])

#Electrolyzer efficiency, capacity and and conversion factors
nu_electrolyzer = var['value']['electrolyzer_eff'] #dimensionless
E_HHV_H2 = var['value']['E_hhv_h2'] #kwh/m^3
HHV_H2 = var['value']['HHV_H2'] #MMBtu/kmol
E_electrolyzer_min = var['value']['min_E_cap'] #kwh
E_electrolyzer_max = var['value']['max_E_cap'] #kwh

#Emission factor of electrolyzer and gasoline vehicle 
EMF_electrolyzer = var['value']['EMF_electrolyzer']
EMF_vehicle = var['value']['emission_gasoline_v'] #tonne CO2/car/year

#Cost realated constants 
beta = var['value']['beta']
C_0 = var['value']['C_0'] #$/kW
mu = var['value']['mu']
gamma = var['value']['gamma']
TC = var['value']['TC'] #$/kWh
C_H2O = var['value']['C_H2O'] #$/L 
WCR = var['value']['water_cons_rate'] #L H2o/m^3 H2
TVM = var['value']['TVM'] #time money value
CAPEX_booster = var['value']['CAPEX_booster'] #$/year
CAPEX_prestorage = var['value']['CAPEX_prestorage'] #$/year
CAPEX_tank = var['value']['CAPEX_tank'] #$/year


#Constants for tank and compressor capacity 
Imax = var['value']['Imax'] #kmol
Imin= var['value']['Imin'] #kmol
Fmax_booster = var['value']['Fmax_booster'] #kmol
Fmax_prestorage =var['value']['Fmax_prestorage'] #kmol

#Electricity consumption rate for prestorage 
ECF_prestorage = var['value']['ECF_prestorage'] #kWh/kmol H2 

#Electricity consumption rate for booster compressor - calculated using power equation
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
Fmax_booster = Fmax_booster * MW_H2 / density_H2 # m^3
Fmax_prestorage = Fmax_prestorage * MW_H2 / density_H2 # m^3


ECF_booster = ECF_booster / MW_H2 * density_H2 #kWh/m^3
ECF_prestorage = ECF_prestorage / MW_H2 * density_H2 #kWh/m^3


#number of electrolyzer max
N_electrolyzer_max = int(3510)


# LP objective variable for transportation model
LP_eps_3 = pulp.LpProblem('LP_eps_3', pulp.LpMaximize)
LP_cost_3 = pulp.LpProblem('LP_cost_3', pulp.LpMinimize)


N_electrolyzer_3 = pulp.LpVariable('N_electrolyzer_3',
                          lowBound=0,
                          cat='Integer')

N_booster_3 = pulp.LpVariable('N_booster_3',
                          lowBound=0,
                          cat='Integer')

N_prestorage_3 = pulp.LpVariable('N_prestorage_3',
                          lowBound=0,
                          cat='Integer')

N_tank_3 = pulp.LpVariable('N_tank_3',
                          lowBound=0,
                          cat='Integer')


alpha_3 = pulp.LpVariable.dicts('alpha_3',
                          [str(i) for i in range(1, N_electrolyzer_max+1)],
                          cat='Binary')


E_3 = pulp.LpVariable.dicts('E_3',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

H2_3 = pulp.LpVariable.dicts('H2_3',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

#hydrogen delivered directly to the booster compressor 
H2_direct_3 = pulp.LpVariable.dicts('H2_direct_3',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

#hydrogen going into the tank 
H2_tank_in_3 = pulp.LpVariable.dicts('H2_tank_in_3',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
#hydrogen leaving the tank
H2_tank_out_3 = pulp.LpVariable.dicts('H2_tank_out_3',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

#hydrogen inventory in the tank
I_H2_3 = pulp.LpVariable.dicts('I_H2_3',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

#compressor emission 
em_compressor_3 = pulp.LpVariable('em_compressor_3',
                          lowBound=0,
                          cat='Continuous')

#electrolyzer emission
em_electrolyzer_3 = pulp.LpVariable('em_electrolyzer_3',
                          lowBound=0,
                          cat='Continuous')

#emission offset
em_offset_3 = pulp.LpVariable('em_offset_3',
                          lowBound=0,
                          cat='Continuous')

#emission for sbg
em_sbg_3 = pulp.LpVariable('em_sbg_3',
                          lowBound=0,
                          cat='Continuous')



#CAPEX and OPEX
CAPEX_3 = pulp.LpVariable('CAPEX_3', lowBound=0, cat='Continuous')
OPEX_3 = pulp.LpVariable('OPEX_3', lowBound=0, cat='Continuous')

# This for loop creates two stage - the first stage is minimizing the emission offset
# The second stage is minimizing the total cost 
for LP in [LP_eps_3,LP_cost_3]:
    for i, h in enumerate([str(i) for i in input_df.index]):

        # Energy and flow constraints
        LP += H2_3[h] == nu_electrolyzer * E_3[h] / E_HHV_H2
        LP += H2_3[h] == H2_tank_in_3[h] + H2_direct_3[h] 

        #hydrogen storage inventory constraint 
        if h == '0': #at hour zero, accumulation assumed to be Imin*Ntank
            LP += I_H2_3[h] == Imin * N_tank_3 + H2_tank_in_3[h] - H2_tank_out_3[h]
        else: #at hour non zero accumulation exists from the previous hour
            LP += I_H2_3[h] == I_H2_3[str(i-1)] + H2_tank_in_3[h] - H2_tank_out_3[h]

        #Demand constraint
        LP += H2_tank_out_3[h] + H2_direct_3[h] == mobility_demand[i] 

        #Electrolyzer constraints 
        LP += N_electrolyzer_3 * E_electrolyzer_min <= E_3[h]
        LP += N_electrolyzer_3 * E_electrolyzer_max >= E_3[h]
        LP += E_3[h] <= SBG[i]
        
        #storage inventory constraint 
        LP += I_H2_3[h] <= Imax * N_tank_3
        LP += I_H2_3[h] >= Imin * N_tank_3
        
        #compressor capacity constraint
        LP += H2_tank_in_3[h] <= N_prestorage_3 * Fmax_prestorage
        LP += H2_tank_out_3[h] + H2_direct_3[h] <= N_booster_3 * Fmax_booster
        
    #Number of eletrolyzer constraint 
    LP += pulp.lpSum(n * alpha_3[str(n)] for n in range(1, N_electrolyzer_max+1)) == N_electrolyzer_3
    LP += pulp.lpSum(alpha_3) <= 1
     
    #Emission calculation 
    LP += pulp.lpSum(EMF[n] * (ECF_booster * (H2_tank_out_3[str(n)] + H2_direct_3[str(n)]) + \
                    ECF_prestorage * H2_tank_in_3[str(n)]) for n in input_df.index)  == em_compressor_3
    

    LP += pulp.lpSum(EMF_electrolyzer * H2_3[h] for h in [str(x) for x in input_df.index]) == em_electrolyzer_3
    LP += pulp.lpSum(EMF[int(h)] * (E_3[h]) for h in [str(x) for x in input_df.index]) == em_sbg_3
    

    
#emission offset by FCV is emission offset due to replacing gasoline vehicle
em_offset_fcv = 100000 * EMF_vehicle

# Epsilon LP Objective
LP_eps_3 += em_offset_fcv - em_compressor_3 - em_electrolyzer_3 - em_sbg_3, 'Offset_3'
LP_eps_3.solve()
print(LP_eps_3.status)
offset_max_3 = LP_eps_3.objective.value()

# Cost of electrolyzer list
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, N_electrolyzer_max+1)]

# CAPEX
LP_cost_3 += pulp.lpSum(alpha_3[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_electrolyzer_max+1)) + \
        (N_booster_3 * CAPEX_booster + \
        N_prestorage_3 * CAPEX_prestorage + \
        N_tank_3 * CAPEX_tank) * 20 == CAPEX_3

# OPEX 
LP_cost_3 += pulp.lpSum((E_3[str(n)] + \
                    ECF_booster * (H2_tank_out_3[str(n)] + H2_direct_3[str(n)]) + \
                    ECF_prestorage * H2_tank_in_3[str(n)]) * (HOEP[n] + TC) for n in input_df.index) + \
        pulp.lpSum(H2_3[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_3

# percentage of maximum emission offset
# This is going to be the contraint in the cost minimization
phi = 0.5

# Cost LP Objective
LP_cost_3 += em_offset_fcv - em_compressor_3 - em_electrolyzer_3 -em_sbg_3 == em_offset_3
LP_cost_3 += em_offset_3 >= phi * offset_max_3
LP_cost_3 += CAPEX_3 + OPEX_3 * TVM, 'Cost_3'


#Estimating the time taken to solve this optimzation problem 
#start time 
start_time_cost = time.time()

print(start_time_cost)

LP_cost_3.solve()


print(LP_cost_3.status)

end_time_cost = time.time()

#time difference 
time_difference_cost = end_time_cost - start_time_cost

print(time_difference_cost)



my_result = create_var_df(LP_cost_3)
my_result = my_result.append({'variable' : 'LP_cost_status', 'value' : LP_cost_3.status} , ignore_index=True)
my_result = my_result.append({'variable' : 'LP_cost_time', 'value' : time_difference_cost} , ignore_index=True)
my_result = my_result.append({'variable' : 'offset_max', 'value' : offset_max_3} , ignore_index=True)
my_result = my_result.append({'variable' : 'phi', 'value' : phi} , ignore_index=True)
filename = 'transportation_result'
export_to_csv(my_result,filename)