import pulp
import pandas as pd

# Time-series constants
SBG = list(input_df['SBG(kWh)']) #kwh
D = list(input_df['industry_demand(m^3)']) #m^3
HOEP = list(input_df['HOEP']) #$/kwh

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


# Fixed constants for industry models 

Fmax_booster = var['value']['Fmax_booster'] #kmol


CAPEX_booster = var['value']['CAPEX_booster'] #$
CAPEX_prestorage = var['value']['CAPEX_prestorage'] #$


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


Fmax_booster = Fmax_booster * MW_H2 / density_H2 # m^3



ECF_booster = ECF_booster / MW_H2 * density_H2 #kWh/m^3

EMF_SMR=EMF_SMR*0.001/MW_H2*density_H2 #tonne CO2/m^3 of H2


#number of electrolyzer max
N_electrolyzer_max = int(5000)



# Transportation model
LP_4 = pulp.LpProblem('LP', pulp.LpMinimize)


N_electrolyzer_4 = pulp.LpVariable('N_electrolyzer_4',
                          lowBound=0,
                          cat='Integer')

N_booster_4 = pulp.LpVariable('N_booster_4',
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
em_offset_max_4 = pulp.LpVariable('em_offset_max_4',
                          lowBound=0,
                          cat='Continuous')
em_old = pulp.LpVariable('em_old',
                          lowBound=0,
                          cat='Continuous')
em_ind = pulp.LpVariable('em_ind',
                          lowBound=0,
                          cat='Continuous')

CAPEX_4 = pulp.LpVariable('CAPEX_4', cat='Continuous')
OPEX_4 = pulp.LpVariable('OPEX_4', cat='Continuous')

for i, h in enumerate([str(i) for i in input_df.index]):
    # Energy and flow constraints
    LP_4 += H2_4[h] == nu_electrolyzer * E_4[h] / E_HHV_H2

    # Demand constraint
    # Actual H2 Demand is 0.05 of total industrial demand
    LP_4 += H2_4[h] == D[i]

    # Electrolyzer constraints
    
    LP_4 += N_electrolyzer_4 * E_electrolyzer_min <= E_4[h]
    LP_4 += N_electrolyzer_4 * E_electrolyzer_max >= E_4[h]
    LP_4 += E_4[h] <= SBG[i]
    
    
    if h == '0':
        LP_4 += pulp.lpSum(n * alpha_4[str(n)] for n in range(1, N_electrolyzer_max+1)) == N_electrolyzer_4
        LP_4 += pulp.lpSum(alpha_4) <= 1

#Emission
LP_4 += pulp.lpSum(EMF_SMR*D[h] for h in input_df.index)==em_old
LP_4 +=pulp.lpSum(EMF_avg*E_4[h]+EMF_electrolyzer*H2_4[h] for h in [str(x) for x in input_df.index])==em_ind
LP_4 +=em_old-em_ind==em_offset_max_4
     
#Compressor Capacity Constraint   
LP_4 += H2_4[h] <= N_booster_4 * Fmax_booster
        
# CAPEX (electrolyzer+Booster pump)
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, N_electrolyzer_max+1)]
LP_4 += pulp.lpSum(alpha_4[str(n)] * C_electrolyzer[n - 1] for n in range(1, N_electrolyzer_max+1))+\
                   N_booster_4*CAPEX_booster*TVM == CAPEX_4

# OPEX (+electericity consumed by running electrolyzer and booster)
LP_4 += pulp.lpSum((E_4[str(n)]+ECF_booster* (H2_4[str(n)])) * (HOEP[n] + TC) for n in input_df.index) + \
      pulp.lpSum(H2_4[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_4

# Objective
LP_4 += CAPEX_4 + OPEX_4 * TVM, 'Cost_4'

LP_4.solve()
print([x.varValue for x in LP_4.variables()])
print(pulp.value(LP_4.objective))