#Setting minimal media

#import cobra and required tools
import cobra
import os
from os.path import join
import matplotlib.pyplot as plt
import pandas as pd
from __future__ import print_function
from cobra import Model, Reaction, Metabolite
import numpy as np
from cobra.medium import minimal_medium

#import iJO1366 - Limo model
data_dir="C:\Work\ecoli_mod"
model=cobra.io.read_sbml_model(join(data_dir,"eco_Lim_v1.xml"))
model.objective="Ec_biomass_iJO1366_core_53p95M"

##minimal medium biomass
cobra.medium.minimal_medium(model)

#minimal medium limonene production
model.objective = 'LIMS'
cobra.medium.minimal_medium(model)

#set minimal media
data_dir="C:\Work\ecoli_mod"
model=cobra.io.read_sbml_model(join(data_dir,"eco_Lim_v1.xml"))
model.objective="Ec_biomass_iJO1366_core_53p95M" 

model.reactions.EX_ca2_e.lower_bound = -0.000521
model.reactions.EX_cl_e.lower_bound = -0.000521
model.reactions.EX_cobalt2_e.lower_bound = -0.000003
model.reactions.EX_cu2_e.lower_bound = -0.000071
model.reactions.EX_fe2_e.lower_bound = -0.000825
model.reactions.EX_fe3_e.lower_bound = -0.000781
model.reactions.EX_glc_e.lower_bound = -1.465741
model.reactions.EX_k_e.lower_bound = -0.019519
model.reactions.EX_mg2_e.lower_bound = -0.000868
model.reactions.EX_mn2_e.lower_bound = -0.000069
model.reactions.EX_mobd_e.lower_bound = -0.000013
model.reactions.EX_nh4_e.lower_bound = -1.080082
model.reactions.EX_ni2_e.lower_bound = -0.000032
model.reactions.EX_o2_e.lower_bound = -2.039876
model.reactions.EX_pi_e.lower_bound = -0.096463
model.reactions.EX_so4_e.lower_bound = -0.025221
model.reactions.EX_zn2_e.lower_bound = -0.000034

#save model in minimal media
cobra.io.write_sbml_model(model,filename="C:/Work/ecoli_mod/eco_Lim_v7.xml")

#Check maximal production for limonene - Define the maximal production of the LIM for the constraints
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_Lim_v7.xml"))

model.objective = "LIMS"
LIMo = model.optimize().f 

#Application
data_dir= "C:/Work/ecoli_mod/"
model = cobra.io.read_sbml_model(join(data_dir, "eco_Lim_v7.xml"))
model.objective = "Ec_biomass_iJO1366_core_53p95M"
 
LIM = list(np.linspace(0, LIMo, 100))
bm = {}
flux = pd.DataFrame()

b= 0
for a in LIM:
    model.reactions.LIMS.upper_bound = a
    model.reactions.LIMS.lower_bound = a
    flux[b] =  model.optimize().fluxes
    bm[b] = model.optimize().f
    b=b+1

#BiomassExport = pd.DataFrame(bm)
Lim = pd.DataFrame(LIM)
bm=pd.Series(bm)
Lim['Biomass']=bm
bm.to_csv("C:/Work/ecoli_mod/LIM_biomass.csv")
Lim.to_csv("C:/Work/ecoli_mod/LIM_biomass.csv") #of Lim constraints in loop & Optimization of biomass

