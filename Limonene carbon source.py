#Changing carbon source

#import cobra and required tools
import cobra
import os
from os.path import join
import matplotlib.pyplot as plt
import pandas as pd
from __future__ import print_function
from cobra import Model, Reaction, Metabolite
import numpy as np

#import iJO1366 - limo model - complex media
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_Lim_v2.xml"))

#biomass with supplementary glycerol (4 g/L)
model.objective = "Ec_biomass_iJO1366_core_53p95M"
model.optimize().f

#biomass no carbohydrate supplement
model.reactions.EX_glc_e.lower_bound = -16.6 (carbohydrate from yeast extract)
model.reactions.EX_glyc_e.lower_bound = 0
model.optimize().f

#biomass with supplementary glucose (4 g/L + carbohydrate from yeast extract as set in complex media) 
model.reactions.EX_glc_e.lower_bound = -38.8
model.reactions.EX_glyc_e.lower_bound = 0
model.optimize().f

#biomass with supplementary malate (4 g/L)
model.reactions.EX_glc_e.lower_bound = -16.6
model.reactions.EX_glyc_e.lower_bound = 0
model.reactions.EX_mal__L_e.lower_bound = -30.3
model.optimize().f

#biomass with supplementary lactate (4 g/L)
model.reactions.EX_glc_e.lower_bound = -16.6
model.reactions.EX_glyc_e.lower_bound = 0
model.reactions.EX_mal__L_e.lower_bound = 0
model.reactions.EX_lac__L_e.lower_bound = -44.9
model.optimize().f

#biomass with supplementary pyruvate (4 g/L)
model.reactions.EX_glc_e.lower_bound = -16.6
model.reactions.EX_glyc_e.lower_bound = 0
model.reactions.EX_lac__L_e.lower_bound = 0
model.reactions.EX_pyr_e.lower_bound = -45.4
model.optimize().f

#limonene with supplementary fructose (4 g/L)
model.reactions.EX_glc_e.lower_bound = -16.6
model.reactions.EX_glyc_e.lower_bound = 0
model.reactions.EX_pyr_e.lower_bound = 0
model.reactions.EX_fru_e.lower_bound = -22.2
model.optimize().f

#limonene with supplementary glycerol (4 g/L)
model.objective = "LIMS"
model.optimize().f
model.metabolites.glyc_c.summary()

#limonene with no carbohydrate supplement
model.reactions.EX_glc_e.lower_bound = -16.6 (carbohydrate from yeast extract)
model.reactions.EX_glyc_e.lower_bound = 0
model.optimize().f

#limonene with supplementary glucose (4 g/L + carbohydrate from yeast extract as set in complex media)
model.reactions.EX_glc_e.lower_bound = -38.8
model.reactions.EX_glyc_e.lower_bound = 0
model.optimize().f

#limonene with supplementary malate (4 g/L)
model.reactions.EX_glc_e.lower_bound = -16.6
model.reactions.EX_glyc_e.lower_bound = 0
model.reactions.EX_mal__L_e.lower_bound = -30.3
model.optimize().f

#limonene with supplementary lactate (4 g/L)
model.reactions.EX_glc_e.lower_bound = -16.6
model.reactions.EX_glyc_e.lower_bound = 0
model.reactions.EX_mal__L_e.lower_bound = 0
model.reactions.EX_lac__L_e.lower_bound = -44.9
model.optimize().f

#limonene with supplementary pyruvate (4 g/L)
model.reactions.EX_glc_e.lower_bound = -16.6
model.reactions.EX_glyc_e.lower_bound = 0
model.reactions.EX_lac__L_e.lower_bound = 0
model.reactions.EX_pyr_e.lower_bound = -45.4
model.optimize().f

#limonene with supplementary fructose (4 g/L)
model.reactions.EX_glc_e.lower_bound = -16.6
model.reactions.EX_glyc_e.lower_bound = 0
model.reactions.EX_pyr_e.lower_bound = 0
model.reactions.EX_fru_e.lower_bound = -22.2
model.optimize().f




