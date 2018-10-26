# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 21:25:38 2018
4- industry model
@author: Wei Yu
"""

import pulp
import pandas as pd

# Time-series constants
SBG = list(input_df['SBG(kWh)'])
D = list(input_df['NG_demand(m^3)'])
HOEP = list(input_df['HOEP'])

# Fixed constants
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

# Industry model
LP_4 = pulp.LpProblem('LP', pulp.LpMinimize)
N_electrolyzer_4 = pulp.LpVariable('N_electrolyzer_4',
                                 lowBound=0,
                                 cat='Continuous')
alpha_4 = pulp.LpVariable.dicts('alpha_4',
                          [str(i) for i in range(1, 31)],
                          cat='Binary')
E_4 = pulp.LpVariable.dicts('E_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
H2_4 = pulp.LpVariable.dicts('H2_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
CO2 = pulp.LpVariable.dicts('CO2_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')
NG_4 = pulp.LpVariable.dicts('NG_4',
                          [str(i) for i in input_df.index],
                          lowBound=0,
                          cat='Continuous')



CAPEX_4 = pulp.LpVariable('CAPEX_4', cat='Continuous')
OPEX_4 = pulp.LpVariable('OPEX_4', cat='Continuous')

for i, h in enumerate([str(i) for i in input_df.index]):
    # Energy and flow constraints
    LP_4 += H2_4[h] == nu_electrolyzer * E_4[h] / E_HHV_H2

    # Demand constraint
    LP_4 += H2_4[h] == D[i]

    # Electrolyzer constraints
    
    LP_4 += N_electrolyzer_4 * E_electrolyzer_min <= E_4[h]
    LP_4 += N_electrolyzer_4 * E_electrolyzer_max >= E_4[h]
    LP_4 += E_4[h] <= SBG[i]
    
    
    if h == '0':
        LP_4 += pulp.lpSum(n * alpha_4[str(n)] for n in range(1, 31)) == N_electrolyzer_4
        LP_4 += pulp.lpSum(alpha_4) <= 1
        
#number of booster constraint
        LP_4 += N_booster_4 <=N_booster_max
        

# CAPEX
C_electrolyzer = [beta * C_0 * i ** mu for i in range(1, 31)]
LP_4 += pulp.lpSum(alpha_4[str(n)] * C_electrolyzer[n - 1] for n in range(1, 31)) == CAPEX_4

# OPEX
LP_4 += pulp.lpSum(E_4[str(n)] * (HOEP[n] + TC) for n in input_df.index) + \
      pulp.lpSum(H2_4[str(n)] * C_H2O * WCR for n in input_df.index) == OPEX_4

# Objective
LP_4 += CAPEX_4 + OPEX_4, 'Cost_4'

LP_4.solve()