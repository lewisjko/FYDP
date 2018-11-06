"""
Model for Transportation
"""

import pulp
import pandas as pd

# Time-series constants
SBG = list(input_df['SBG(kWh)']) #kwh
mobility_demand = list(input_df['mobility_demand(m^3)']) #m^3
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
EMF_vehicle = var['value']['emission_gasoline_v'] #tonne CO2/car/year
num_vehicle = var['value']['N_gasoline_v']
FCV_penetration = var['value']['FCV_penetration'] #FCV market penetration (estimated)


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
Fmax_booster = Fmax_booster * MW_H2 / density_H2 # m^3
Fmax_prestorage = Fmax_prestorage * MW_H2 / density_H2 # m^3


ECF_booster = ECF_booster / MW_H2 * density_H2 #kWh/m^3
ECF_prestorage = ECF_prestorage / MW_H2 * density_H2 #kWh/m^3


#number of electrolyzer max
N_electrolyzer_max = int(10000)



# Transportation model
LP_3 = pulp.LpProblem('LP', pulp.LpMinimize)


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
H2_direct = pulp.LpVariable.dicts('H2_direct',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

#hydrogen going into the tank 
H2_tank_in = pulp.LpVariable.dicts('H2_tank_in',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
#hydrogen leaving the tank
H2_tank_out = pulp.LpVariable.dicts('H2_tank_out',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

#hydrogen inventory in the tank
I_H2 = pulp.LpVariable.dicts('I_H2',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')

em_compressor = pulp.LpVariable('em_compressor',
                          lowBound=0,
                          cat='Continuous')


CAPEX_3 = pulp.LpVariable('CAPEX_3', lowBound=0, cat='Continuous')
OPEX_3 = pulp.LpVariable('OPEX_3', lowBound=0, cat='Continuous')

for i, h in enumerate([str(i) for i in input_df.index]):
    
    # Energy and flow constraints
    LP_3 += H2_3[h] == nu_electrolyzer * E_3[h] / E_HHV_H2
    LP_3 += H2_3[h] == H2_tank_in[h] + H2_direct[h] 

    #hydrogen storage inventory constraint 
    if h == '0': #at hour zero, accumulation assumed to be Imin*Ntank
        LP_3 += I_H2[h] == Imin * N_tank + H2_tank_in[h] - H2_tank_out[h]
    else: #at hour non zero accumulation exists 
        LP_3 += I_H2[h] == I_H2[str(i-1)] + H2_tank_in[h] - H2_tank_out[h]
    

    # Demand constraint
    LP_3 += H2_tank_out[h] + H2_direct[h] == mobility_demand[i] 

    # Electrolyzer constraints 
    LP_3 += N_electrolyzer_3 * E_electrolyzer_min <= E_3[h]
    LP_3 += N_electrolyzer_3 * E_electrolyzer_max >= E_3[h]
    LP_3 += E_3[h] <= SBG[i]
    
    
    if h == '0':
        #Number of eletrolyzer constraint 
        LP_3 += pulp.lpSum(n * alpha_3[str(n)] for n in range(1, N_electrolyzer_max+1)) == N_electrolyzer_3
        LP_3 += pulp.lpSum(alpha_3) <= 1
         
    #storage inventory constraint 
    LP_3 += I_H2[h] <= Imax * N_tank  
    LP_3 += I_H2[h] >= Imin * N_tank

    #compressor capacity constraint
    LP_3 += H2_tank_in[h] <= N_prestorage * Fmax_prestorage
    LP_3 += H2_tank_out[h] + H2_direct[h] <= N_booster * Fmax_booster

    
# Cost of electrolyzer list
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, N_electrolyzer_max+1)]

# CAPEX
LP_3 += pulp.lpSum(alpha_3[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_electrolyzer_max+1)) + \
        (N_booster * CAPEX_booster + \
        N_prestorage * CAPEX_prestorage + \
        N_tank * CAPEX_tank) * 20 == CAPEX_3


LP_3 += pulp.lpSum((E_3[str(n)] + \
                    ECF_booster * (H2_tank_out[str(n)] + H2_direct[str(n)]) + \
                    ECF_prestorage * H2_tank_in[str(n)]) * (HOEP[n] + TC) for n in input_df.index) + \
        pulp.lpSum(H2_3[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_3


em_offset_fcv = num_vehicle * FCV_penetration * EMF_vehicle

LP_3 += pulp.lpSum(EMF[n] * (ECF_booster * (H2_tank_out[str(n)] + H2_direct[str(n)]) + \
                    ECF_prestorage * H2_tank_in[str(n)]) for n in input_df.index)  == em_compressor


# Objective
LP_3 += CAPEX_3 + OPEX_3 * TVM, 'Cost_3'

LP_3.solve()
print([x.varValue for x in LP_3.variables()])
print(pulp.value(LP_3.objective))
print(LP_3.status)
print(N_electrolyzer_3.value())
print(N_booster.value())
print(N_prestorage.value())
print(N_tank.value())
print(CAPEX_3.value(), OPEX_3.value())

