# -*- coding: utf-8 -*-

import pandas as pd
var = pd.read_csv('Variables.csv')
var.set_index('variable_name', inplace=True)
