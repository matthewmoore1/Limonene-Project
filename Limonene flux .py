#Limonene flux analysis

#Import cobra and required tools
import cobra
import os
from os.path import join
import matplotlib.pyplot as plt
import pandas as pd
from __future__ import print_function
from cobra import Model, Reaction, Metabolite
import numpy as np

#import iJO1366 - Limo model - complex media
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_Lim_v2.xml"))
#set limonene flux as objective function
model.objective = "LIMS"

#removal of MEV pathway
model.reactions.HMGCOAS.upper_bound = 0
model.optimize()
model.summary()
model.metabolites.co2_c.summary()

#restore MEV pathway
model.reactions.HMGCOAS.upper_bound = 1000

#removal of DXP pathway
model.reactions.DXPS.upper_bound = 0

#optimized for limonene production without DXP pathway
model.optimize()
model.summary()
model.metabolites.co2_c.summary()

#restore DXP pathway
model.reactions.DXPS.upper_bound = 1000

#removal of individual metabolites on biomass
model.objective = "Ec_biomass_iJO1366_core_53p95M"
model.reactions.EX_#metabolite.lower_bound = 0
model.optimize()
model.reactions.EX_#metabolite.lower_bound = #restore flux to that in complex media
#repeat for every metabolite in complex media

#removal of individual metabolites on limonene
model.objective # "LIMS"
model.reactions.EX_#metabolite.lower_bound = 0
model.optimize()
model.reactions.EX_#metabolite.lower_bound = #restore flux in complex media
#repeat for every metabolite in complex media

#increasing flux of individual metabolites on biomass
model.objective = "Ec_biomass_iJO1366_core_53p95M"
model.reactions.EX_#metabolite.lower_bound = #triple flux of metabolite in complex media
model.optimize()
model.reactions.EX_#metabolite.lower_bound = #restore flux to that in complex media

#increasing flux of individual metabolites on limonene
model.objective = "Ec_biomass_iJO1366_core_53p95M"
model.reactions.EX_#metabolite.lower_bound = #triple flux of metabolite in complex media
model.optimize()
model.reactions.EX_#metabolite.lower_bound = #restore flux to that in complex media

#metabolite summary command used throughout
model.metabolites.#metabolite_c.summary()